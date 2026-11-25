"""多次独立曝光的整批关联与测光零点校正.

每次曝光独立求解:候选沿用原球面 haversine 浮点距离与原半径判定;曝光内
每条检测与每颗参考星至多参与一对,先最大化配对数,再最小化整数边代价之和
(毫角秒,floor(浮点角距离*3600000+0.5))。全部最优方案用最小费用流的零约化
费用残差图强连通分量刻画,不枚举方案。测光零点只取在所有几何最优方案中都
存在且参考星标为 calibrator 的配对(锚点),偏移/中位数/支持度/星等残差均按
十进制输入的精确值(Fraction)计算;零点有效时仅在几何最优方案内再按绝对
星等残差之和精化,测光引起的关联变化不回流参与定标。
"""
from math import sin,cos,atan2,sqrt,radians,degrees,floor
from heapq import heappush,heappop
from .checks import fields,array,ident,number

def _distance(ra1,dec1,ra2,dec2):
    a,b=radians(float(dec1)),radians(float(dec2))
    h=sin((a-b)/2)**2+cos(a)*cos(b)*sin(radians(float(ra1-ra2))/2)**2
    h=min(1.0,max(0.0,h))
    return degrees(2*atan2(sqrt(h),sqrt(1-h)))

def _decimal_str(fr):
    if fr<0:return '-'+_decimal_str(-fr)
    num,den=fr.numerator,fr.denominator
    twos=0
    while den%2==0:den//=2;twos+=1
    fives=0
    while den%5==0:den//=5;fives+=1
    if den!=1:raise ValueError('non-terminating decimal')
    k=max(twos,fives)
    scaled=num*2**(k-twos)*5**(k-fives)
    if not k:return str(scaled)
    digits=str(scaled).rjust(k+1,'0')
    return (digits[:-k]+'.'+digits[-k:]).rstrip('0').rstrip('.')

def _min_cost_max_flow(size,arcs,source,sink,zero):
    """逐次最短路最小费用最大流;费用为按字典序比较的等长元组.

    arcs: (u,v,cost) 列表,容量均为 1。返回 (graph,potential,located);
    graph[u]=[[to,rev,cap,cost],...](含反向弧),potential 为最优对偶势
    (全部残差弧约化费用非负),located[k]=(u,pos) 给出第 k 条输入弧位置。
    """
    graph=[[] for _ in range(size)]
    located=[]
    for u,v,cost in arcs:
        located.append((u,len(graph[u])))
        graph[u].append([v,len(graph[v]),1,cost])
        graph[v].append([u,len(graph[u])-1,0,tuple(-c for c in cost)])
    potential=[zero]*size
    while True:
        dist=[None]*size;prev=[None]*size
        dist[source]=zero;heap=[(zero,source)]
        while heap:
            d,u=heappop(heap)
            if d!=dist[u]:continue
            for i,e in enumerate(graph[u]):
                if e[2]<=0:continue
                w=e[0]
                nd=tuple(a+b+c for a,b,c in zip(d,e[3],potential[u]))
                nd=tuple(a-b for a,b in zip(nd,potential[w]))
                if dist[w] is None or nd<dist[w]:
                    dist[w]=nd;prev[w]=(u,i);heappush(heap,(nd,w))
        if dist[sink] is None:break
        for v in range(size):
            if dist[v] is not None:
                potential[v]=tuple(a+b for a,b in zip(potential[v],dist[v]))
        v=sink
        while v!=source:
            u,i=prev[v];e=graph[u][i]
            e[2]-=1;graph[v][e[1]][2]+=1;v=u
    return graph,potential,located

def _feasible_potential(graph,zero):
    """为残差图重算可行势位(虚拟超源最短路),使全部残差弧约化费用非负。

    逐次最短路的对偶势位在"不可达→可达"的残差弧上可能为负;最优点无负费用
    残差环,故最短路势位存在且唯一地给出零费用环流刻画所需的非负约化费用。
    """
    size=len(graph)
    dist=[zero]*size
    for _ in range(size-1):
        updated=False
        for u in range(size):
            du=dist[u]
            for e in graph[u]:
                if e[2]>0:
                    w=e[0]
                    nd=tuple(a+b for a,b in zip(du,e[3]))
                    if nd<dist[w]:dist[w]=nd;updated=True
        if not updated:break
    return dist

def _zero_components(graph,potential,zero):
    """零约化费用残差图的强连通分量(迭代 Tarjan)。"""
    size=len(graph)
    adj=[[] for _ in range(size)]
    for u in range(size):
        pu=potential[u]
        for e in graph[u]:
            if e[2]>0:
                rc=tuple(a+b-c for a,b,c in zip(e[3],pu,potential[e[0]]))
                if rc==zero:adj[u].append(e[0])
    index=[None]*size;low=[0]*size;comp=[None]*size
    stack=[];onstack=[False]*size
    counter=0;ncomp=0
    for root in range(size):
        if index[root] is not None:continue
        work=[(root,0)]
        while work:
            v,pi=work[-1]
            if pi==0:
                index[v]=low[v]=counter;counter+=1
                stack.append(v);onstack[v]=True
            descended=False
            while pi<len(adj[v]):
                w=adj[v][pi];pi+=1
                if index[w] is None:
                    work[-1]=(v,pi);work.append((w,0))
                    descended=True;break
                if onstack[w] and index[w]<low[v]:low[v]=index[w]
            if descended:continue
            work.pop()
            if low[v]==index[v]:
                while True:
                    w=stack.pop();onstack[w]=False;comp[w]=ncomp
                    if w==v:break
                ncomp+=1
            if work:
                p=work[-1][0]
                if low[v]<low[p]:low[p]=low[v]
    return comp

def _arc_states(graph,potential,comp,zero,u,pos):
    """单位容量弧 (u,pos) 的 (flow,can_be_one,can_be_zero).

    最优解集合恰为当前最优流叠加零约化费用残差图上的环流;某条弧能在某个
    最优解中改变取值,当且仅当它(流为 0 时正向弧,流为 1 时反向弧)落在零
    费用残差环上,即弧为零约化费用且两端点同属一个强连通分量。
    """
    e=graph[u][pos];v=e[0]
    if e[2]==0:
        b=graph[v][e[1]]
        rc=tuple(a+b-c for a,b,c in zip(b[3],potential[v],potential[u]))
        return 1,True,rc==zero and comp[v]==comp[u]
    rc=tuple(a+b-c for a,b,c in zip(e[3],potential[u],potential[v]))
    return 0,rc==zero and comp[u]==comp[v],True

def _validate(data):
    fields(data,('reference','threshold','anchor_tolerance','exposures'))
    threshold=number(data['threshold'])
    if threshold<=0:raise ValueError('positive threshold required')
    tolerance=number(data['anchor_tolerance'])
    if tolerance<0:raise ValueError('non-negative anchor_tolerance required')
    reference={}
    for x in array(data['reference'],1,1000):
        fields(x,('id','ra','dec','magnitude','calibrator'))
        k=ident(x['id'])
        if k in reference:raise ValueError('duplicate reference')
        if type(x['calibrator']) is not bool:raise ValueError('calibrator must be boolean')
        ra,dec,mag=number(x['ra']),number(x['dec']),number(x['magnitude'])
        if not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid sky coordinates')
        reference[k]=(ra,dec,mag,x['calibrator'])
    exposures=[];exp_ids=set();total=0
    for e in array(data['exposures'],1,2000):
        fields(e,('id','observations'))
        ek=ident(e['id'])
        if ek in exp_ids:raise ValueError('duplicate exposure')
        exp_ids.add(ek)
        observations=[];obs_ids=set()
        for x in array(e['observations'],1,2000):
            fields(x,('id','ra','dec','sigma','instrumental_magnitude'))
            k=ident(x['id'])
            if k in obs_ids:raise ValueError('duplicate observation')
            obs_ids.add(k)
            ra,dec=number(x['ra']),number(x['dec'])
            sigma,instr=number(x['sigma']),number(x['instrumental_magnitude'])
            if sigma<=0:raise ValueError('sigma must be positive')
            if not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid sky coordinates')
            observations.append({'id':k,'ra':ra,'dec':dec,'sigma':sigma,'instrumental':instr})
        total+=len(observations)
        exposures.append({'id':ek,'observations':observations})
    if not 1<=total<=2000:raise ValueError('total observations 1..2000 required')
    return threshold,tolerance,reference,exposures

def _process_exposure(exp,threshold,tolerance,reference,anchor_blocked=None):
    observations=exp['observations'];n=len(observations)
    blocked_anchors=set() if anchor_blocked is None else set(anchor_blocked)
    refs=sorted(reference.items());m=len(refs)
    candidates=[]
    for o in observations:
        radius=threshold*o['sigma'];entries=[]
        for rid,(ra,dec,mag,cal) in refs:
            dist=_distance(o['ra'],o['dec'],ra,dec)
            if dist<=radius:entries.append((dist,rid,floor(dist*3600000+0.5)))
        entries.sort(key=lambda z:(z[0],z[1]))
        candidates.append(entries)
    size=n+m+2;source=0;sink=size-1
    refnode={rid:1+n+j for j,(rid,_) in enumerate(refs)}
    def solve(zeropoint):
        zero=(0,) if zeropoint is None else (0,0)
        arcs=[];kind=[]
        for i in range(n):
            arcs.append((source,1+i,zero));kind.append(('s',i))
        for i in range(n):
            for dist,rid,cost in candidates[i]:
                if zeropoint is None:edge=(cost,)
                else:edge=(cost,abs(observations[i]['instrumental']+zeropoint-reference[rid][2]))
                arcs.append((1+i,refnode[rid],edge));kind.append(('c',i,rid))
        for j,(rid,_) in enumerate(refs):
            arcs.append((1+n+j,sink,zero));kind.append(('r',j))
        graph,potential,located=_min_cost_max_flow(size,arcs,source,sink,zero)
        potential=_feasible_potential(graph,zero)
        comp=_zero_components(graph,potential,zero)
        states=[_arc_states(graph,potential,comp,zero,u,pos) for u,pos in located]
        return states,kind
    states,kind=solve(None)
    anchors=[]
    for idx,(tag,*rest) in enumerate(kind):
        if tag=='c':
            i,rid=rest
            flow,can1,can0=states[idx]
            if flow and not can0 and reference[rid][3] and observations[i]['id'] not in blocked_anchors:
                anchors.append((observations[i]['id'],rid,reference[rid][2]-observations[i]['instrumental']))
    zeropoint=None
    if len(anchors)>=3:
        offs=sorted(a[2] for a in anchors);q=len(offs)
        median=offs[q//2] if q%2 else (offs[q//2-1]+offs[q//2])/2
        supporters=sum(1 for a in anchors if abs(a[2]-median)<=tolerance)
        if supporters>=3 and 2*supporters>len(anchors):
            zeropoint=median;status='calibrated'
        else:status='no_consensus'
    else:status='insufficient_anchors'
    if zeropoint is not None:states,kind=solve(zeropoint)
    src_idx=[None]*n;cand_idx=[[] for _ in range(n)]
    for idx,(tag,*rest) in enumerate(kind):
        if tag=='s':src_idx[rest[0]]=idx
        elif tag=='c':cand_idx[rest[0]].append(idx)
    detections=[]
    for i,o in enumerate(observations):
        may_unassigned=states[src_idx[i]][2]
        possible=sorted(rid for idx,(dist,rid,cost) in zip(cand_idx[i],candidates[i]) if states[idx][1])
        selected=possible[0] if len(possible)==1 and not may_unassigned else None
        if selected is not None:st='matched'
        elif not candidates[i]:st='unmatched'
        elif not possible:st='unassigned'
        else:st='ambiguous'
        corrected=None if zeropoint is None else _decimal_str(o['instrumental']+zeropoint)
        detections.append({'observation':o['id'],
            'candidates':[{'reference':rid,'distance':str(dist),'cost':cost} for dist,rid,cost in candidates[i]],
            'possible':possible,
            'may_be_unassigned':may_unassigned,
            'selected':selected,
            'status':st,
            'corrected_magnitude':corrected})
    detections.sort(key=lambda x:x['observation'])
    anchor_list=sorted(({'observation':oid,'reference':rid,'offset':_decimal_str(off)} for oid,rid,off in anchors),key=lambda x:(x['observation'],x['reference']))
    return {'exposure':exp['id'],
        'calibration':{'status':status,
            'zero_point':None if zeropoint is None else _decimal_str(zeropoint),
            'anchors':anchor_list},
        'detections':detections}

def process_exposures(data):
    """整批处理多次曝光,返回可 JSON 序列化的关联与校正测光结果。"""
    threshold,tolerance,reference,exposures=_validate(data)
    out=[_process_exposure(e,threshold,tolerance,reference) for e in exposures]
    out.sort(key=lambda x:x['exposure'])
    return {'exposures':out}

def _geometric_evidence(threshold,tolerance,reference,exposures):
    """只按几何最优方案签发身份,供测光发布件使用。

    与 ``_process_exposure`` 共用同一套候选半径、整数角距代价、最小费用流与
    零费用残差环流刻画,但完全不做星等破平局:一条参考源-检测边当且仅当在
    *所有*几何最优方案中都被选用(flow=1 且 can0=False)时,才把该检测签给
    该参考源。返回每帧的几何强制配对(含原始仪器星等)、每个参考源的原始
    几何候选检测 ID 集合,以及完全沿用原锚点规则得到的有效绝对零点(无效为
    None)。转移零点不回流改变锚点或身份。
    """
    out=[]
    for exp in exposures:
        observations=exp['observations'];n=len(observations)
        refs=sorted(reference.items());m=len(refs)
        candidates=[]
        for o in observations:
            radius=threshold*o['sigma'];entries=[]
            for rid,(ra,dec,mag,cal) in refs:
                dist=_distance(o['ra'],o['dec'],ra,dec)
                if dist<=radius:entries.append((dist,rid,floor(dist*3600000+0.5)))
            entries.sort(key=lambda z:(z[0],z[1]))
            candidates.append(entries)
        size=n+m+2;source=0;sink=size-1
        refnode={rid:1+n+j for j,(rid,_) in enumerate(refs)}
        zero=(0,);arcs=[];kind=[]
        for i in range(n):
            arcs.append((source,1+i,zero));kind.append(('s',i))
        for i in range(n):
            for dist,rid,cost in candidates[i]:
                arcs.append((1+i,refnode[rid],(cost,)));kind.append(('c',i,rid))
        for j,(rid,_) in enumerate(refs):
            arcs.append((1+n+j,sink,zero));kind.append(('r',j))
        graph,potential,located=_min_cost_max_flow(size,arcs,source,sink,zero)
        potential=_feasible_potential(graph,zero)
        comp=_zero_components(graph,potential,zero)
        states=[_arc_states(graph,potential,comp,zero,u,pos) for u,pos in located]
        forced={}
        for idx,(tag,*rest) in enumerate(kind):
            if tag=='c':
                i,rid=rest
                flow,can1,can0=states[idx]
                if flow and not can0:
                    forced[rid]=(observations[i]['id'],observations[i]['instrumental'])
        incidence={}
        for i,entries in enumerate(candidates):
            for dist,rid,cost in entries:
                incidence.setdefault(rid,set()).add(observations[i]['id'])
        anchors=[]
        for rid,(oid,instr) in forced.items():
            if reference[rid][3]:
                anchors.append((oid,rid,reference[rid][2]-instr))
        zeropoint=None
        if len(anchors)>=3:
            offs=sorted(a[2] for a in anchors);q=len(offs)
            median=offs[q//2] if q%2 else (offs[q//2-1]+offs[q//2])/2
            supporters=sum(1 for a in anchors if abs(a[2]-median)<=tolerance)
            if supporters>=3 and 2*supporters>len(anchors):
                zeropoint=median
        out.append({'exposure':exp['id'],'forced':forced,
            'incidence':{rid:sorted(ids) for rid,ids in incidence.items()},
            'zero_point':zeropoint})
    return out
