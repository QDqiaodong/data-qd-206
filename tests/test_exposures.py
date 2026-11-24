import copy,json,random,unittest
from fractions import Fraction
from math import floor
from star_crossmatch import process_exposures
from star_crossmatch.exposures import _distance,_decimal_str

def ref(rid,ra,dec,mag,cal=False):
    return {'id':rid,'ra':ra,'dec':dec,'magnitude':mag,'calibrator':cal}

def obs(oid,ra,dec,sigma,instr):
    return {'id':oid,'ra':ra,'dec':dec,'sigma':sigma,'instrumental_magnitude':instr}

def batch(reference,exposures,threshold='1',tolerance='0.5'):
    return {'reference':reference,'threshold':threshold,
        'anchor_tolerance':tolerance,'exposures':exposures}

def one_exposure(result,eid='e1'):
    assert len(result['exposures'])==1
    exp=result['exposures'][0]
    assert exp['exposure']==eid
    return exp

def detmap(exp):
    return {d['observation']:d for d in exp['detections']}

class Association(unittest.TestCase):
    def test_forced_pair_cross_ra_and_reuse(self):
        # o1 更靠近 A 但只能二选一配成两对:o1->B、o2->A;跨赤经 0 度的 o3->C。
        data=batch([
            ref('A','10','0','20'),ref('B','10.0003','0','20'),ref('C','359.95','0','20')],[
            {'id':'e1','observations':[
                obs('o1','10.0001','0','0.00025','10'),
                obs('o2','10','0','0.00025','10'),
                obs('o3','0.05','0','0.2','10')]},
            {'id':'e2','observations':[obs('o1','10','0','0.00025','10')]}])
        r=process_exposures(data)
        self.assertEqual([e['exposure'] for e in r['exposures']],['e1','e2'])
        d=detmap(r['exposures'][0])
        self.assertEqual(d['o1']['selected'],'B')
        self.assertEqual(d['o2']['selected'],'A')
        self.assertEqual(d['o3']['selected'],'C')
        self.assertEqual(d['o1']['possible'],['B'])
        self.assertEqual(d['o2']['possible'],['A'])
        self.assertFalse(d['o1']['may_be_unassigned'])
        self.assertEqual(d['o1']['status'],'matched')
        # 角距证据按原球面判定,跨赤经 0 度;整数代价毫角秒
        costs={c['reference']:c['cost'] for c in d['o1']['candidates']}
        self.assertEqual(costs,{'A':360,'B':720})
        dist={c['reference']:c['distance'] for c in d['o3']['candidates']}
        self.assertAlmostEqual(float(dist['C']),0.1,places=10)
        self.assertEqual(d['o3']['candidates'][0]['cost'],360000)
        # 另一曝光可再次使用同一参考星 A;观测 ID 只在本曝光内唯一
        d2=detmap(r['exposures'][1])
        self.assertEqual(d2['o1']['selected'],'A')
        # 无定标星:锚点不足,校正值为空
        self.assertEqual(r['exposures'][0]['calibration'],
            {'status':'insufficient_anchors','zero_point':None,'anchors':[]})
        self.assertIsNone(d['o1']['corrected_magnitude'])

    def test_geometry_dominates_zero_residual(self):
        # 三对孤立定标锚点(偏移 +2)给出零点;争夺组几何严格最优为 o1->A、o2->B,
        # 而测光更偏好几何更贵的 o1->B、o2->A(残差 2<4),几何必须优先。
        data=batch([
            ref('C1','50','0','15',True),ref('C2','60','0','15',True),ref('C3','70','0','15',True),
            ref('A','10','0','18'),ref('B','10.000002','0','15')],[
            {'id':'e1','observations':[
                obs('k1','50','0','0.001','13'),obs('k2','60','0','0.001','13'),obs('k3','70','0','0.001','13'),
                obs('o1','10','0','0.00001','13'),obs('o2','10.000002','0','0.00001','14')]}])
        exp=one_exposure(process_exposures(data))
        self.assertEqual(exp['calibration']['status'],'calibrated')
        self.assertEqual(exp['calibration']['zero_point'],'2')
        d=detmap(exp)
        self.assertEqual(d['o1']['selected'],'A')
        self.assertEqual(d['o2']['selected'],'B')
        self.assertEqual(d['o1']['possible'],['A'])
        self.assertEqual(d['o2']['possible'],['B'])
        self.assertEqual(d['o1']['corrected_magnitude'],'15')
        self.assertEqual(d['o2']['corrected_magnitude'],'16')
        costs={c['reference']:c['cost'] for c in d['o1']['candidates']}
        self.assertEqual(costs,{'A':0,'B':7})

class Calibration(unittest.TestCase):
    def base_data(self):
        return batch([
            ref('C1','50','0','15',True),ref('C2','60','0','15',True),ref('C3','70','0','15',True),
            ref('N1','80','0','15'),ref('N2','80','0','17'),ref('N3','90','0','16')],[
            {'id':'e1','observations':[
                obs('k1','50','0','0.001','13'),obs('k2','60','0','0.001','13'),obs('k3','70','0','0.001','13'),
                obs('q','80','0','0.001','13'),
                obs('z','85','0','0.001','13')]}])

    def test_zero_point_selects_bright_reference(self):
        exp=one_exposure(process_exposures(self.base_data()))
        cal=exp['calibration']
        self.assertEqual(cal['status'],'calibrated')
        self.assertEqual(cal['zero_point'],'2')
        self.assertEqual([(a['observation'],a['reference'],a['offset']) for a in cal['anchors']],
            [('k1','C1','2'),('k2','C2','2'),('k3','C3','2')])
        d=detmap(exp)
        # 校正星等 15:只配星等 15 的 N1,星等 17 的 N2 被测光精化排除
        self.assertEqual(d['q']['selected'],'N1')
        self.assertEqual(d['q']['possible'],['N1'])
        self.assertEqual(d['q']['corrected_magnitude'],'15')
        self.assertEqual({c['reference'] for c in d['q']['candidates']},{'N1','N2'})
        # 无候选检测同样交付校正值
        self.assertEqual(d['z']['status'],'unmatched')
        self.assertIsNone(d['z']['selected'])
        self.assertEqual(d['z']['corrected_magnitude'],'15')
        self.assertEqual(d['k1']['corrected_magnitude'],'15')

    def test_anchor_loss_clears_zero_point(self):
        # 只移动非定标星 N3 与锚点 C1 同位置,成为 k1 的等价替代:
        # 锚点 3->2 不足三个,零点置空,远处检测恢复歧义,无候选检测校正值变空。
        data=self.base_data()
        data['reference'][5]['ra']='50'
        exp=one_exposure(process_exposures(data))
        cal=exp['calibration']
        self.assertEqual(cal['status'],'insufficient_anchors')
        self.assertIsNone(cal['zero_point'])
        self.assertEqual([(a['observation'],a['reference']) for a in cal['anchors']],
            [('k2','C2'),('k3','C3')])
        d=detmap(exp)
        self.assertEqual(d['q']['possible'],['N1','N2'])
        self.assertFalse(d['q']['may_be_unassigned'])
        self.assertIsNone(d['q']['selected'])
        self.assertEqual(d['q']['status'],'ambiguous')
        self.assertIsNone(d['q']['corrected_magnitude'])
        self.assertIsNone(d['z']['corrected_magnitude'])
        self.assertEqual(d['k1']['possible'],['C1','N3'])
        self.assertEqual(d['k1']['status'],'ambiguous')

    def test_majority_required(self):
        # 9 个锚点分三簇 {2,2,2,5,5,5,9,9,9}:中位数 5 只有 3 个支持者,非严格多数。
        refs=[];ob=[]
        for i,mag in enumerate([15,15,15,18,18,18,22,22,22]):
            refs.append(ref(f'C{i}','%d'%(10*(i+1)),'0',str(mag),True))
            ob.append(obs(f'k{i}','%d'%(10*(i+1)),'0','0.001','13'))
        exp=one_exposure(process_exposures(batch(refs,[{'id':'e1','observations':ob}])))
        self.assertEqual(exp['calibration']['status'],'no_consensus')
        self.assertIsNone(exp['calibration']['zero_point'])
        self.assertEqual(len(exp['calibration']['anchors']),9)
        # 3 个锚点仅 2 个支持:同样一致性不足
        refs2=[ref('C1','10','0','15',True),ref('C2','20','0','15',True),ref('C3','30','0','22',True)]
        ob2=[obs('k1','10','0','0.001','13'),obs('k2','20','0','0.001','13'),obs('k3','30','0','0.001','13')]
        exp2=one_exposure(process_exposures(batch(refs2,[{'id':'e1','observations':ob2}])))
        self.assertEqual(exp2['calibration']['status'],'no_consensus')

    def test_even_median_mean(self):
        # 偶数锚点取中间两值均值:偏移 {2,2,3,3} -> 零点 2.5
        refs=[ref('C1','10','0','15',True),ref('C2','20','0','15',True),
            ref('C3','30','0','16',True),ref('C4','40','0','16',True)]
        ob=[obs(f'k{i}','%d'%(10*(i+1)),'0','0.001','13') for i in range(4)]
        exp=one_exposure(process_exposures(batch(refs,[{'id':'e1','observations':ob}])))
        self.assertEqual(exp['calibration']['status'],'calibrated')
        self.assertEqual(exp['calibration']['zero_point'],'2.5')
        d=detmap(exp)
        self.assertEqual(d['k1']['corrected_magnitude'],'15.5')

class Ties(unittest.TestCase):
    def test_tied_contention_both_may_be_unassigned(self):
        # 两条检测同价争同一颗星、无有效零点:都可能配到,也都可能留空。
        data=batch([ref('R','10','0','20')],[
            {'id':'e1','observations':[obs('p1','10','0','0.001','13'),obs('p2','10','0','0.001','13')]}])
        exp=one_exposure(process_exposures(data))
        d=detmap(exp)
        for oid in ('p1','p2'):
            self.assertEqual(d[oid]['possible'],['R'])
            self.assertTrue(d[oid]['may_be_unassigned'])
            self.assertIsNone(d[oid]['selected'])
            self.assertEqual(d[oid]['status'],'ambiguous')
            self.assertIsNone(d[oid]['corrected_magnitude'])

    def test_crowded_band(self):
        # 40 个连通小组(2^40 个最优方案,不可枚举):相邻组的星进候选,
        # 但最终每条检测只能配本组任一颗星且不能留空。
        refs=[];ob=[]
        for i in range(40):
            ra='0.%03d'%(9*i)
            refs.append(ref(f'g{i}a',ra,'0','20'))
            refs.append(ref(f'g{i}b',ra,'0','20'))
            ob.append(obs(f'd{i}x',ra,'0','0.011','10'))
            ob.append(obs(f'd{i}y',ra,'0','0.011','10'))
        exp=one_exposure(process_exposures(batch(refs,[{'id':'band','observations':ob}])),eid='band')
        self.assertEqual(exp['calibration'],
            {'status':'insufficient_anchors','zero_point':None,'anchors':[]})
        d=detmap(exp)
        self.assertEqual(len(d),80)
        for i in range(40):
            for s in ('x','y'):
                row=d[f'd{i}{s}']
                self.assertEqual(row['possible'],[f'g{i}a',f'g{i}b'])
                self.assertFalse(row['may_be_unassigned'])
                self.assertIsNone(row['selected'])
                self.assertEqual(row['status'],'ambiguous')
                got={c['reference']:c['cost'] for c in row['candidates']}
                want={f'g{i}a':0,f'g{i}b':0}
                if i>0:want.update({f'g{i-1}a':32400,f'g{i-1}b':32400})
                if i<39:want.update({f'g{i+1}a':32400,f'g{i+1}b':32400})
                self.assertEqual(got,want)

class Validation(unittest.TestCase):
    def minimal(self):
        return batch([ref('r1','10','0','15',True)],
            [{'id':'e1','observations':[obs('o1','10','0','0.1','13')]}])

    def expect_bad(self,mutate):
        data=self.minimal();mutate(data)
        with self.assertRaises(ValueError):process_exposures(data)

    def test_fields_and_types(self):
        self.expect_bad(lambda d:d.pop('anchor_tolerance'))
        self.expect_bad(lambda d:d.update(extra='1'))
        self.expect_bad(lambda d:d['reference'][0].pop('calibrator'))
        self.expect_bad(lambda d:d['reference'][0].update(calibrator='true'))
        self.expect_bad(lambda d:d['reference'][0].update(calibrator=1))
        self.expect_bad(lambda d:d['exposures'][0]['observations'][0].pop('instrumental_magnitude'))
        self.expect_bad(lambda d:d['exposures'][0]['observations'][0].update(color='red'))
        self.expect_bad(lambda d:d['exposures'][0].pop('id'))
        self.expect_bad(lambda d:d['reference'][0].update(magnitude=15.0))
        self.expect_bad(lambda d:d['exposures'][0]['observations'][0].update(id='bad id!'))

    def test_bounds_and_identities(self):
        self.expect_bad(lambda d:d.update(anchor_tolerance='-0.5'))
        self.expect_bad(lambda d:d.update(threshold='0'))
        self.expect_bad(lambda d:d['reference'][0].update(ra='360'))
        self.expect_bad(lambda d:d['reference'][0].update(dec='-90.000001'))
        self.expect_bad(lambda d:d['exposures'][0]['observations'][0].update(sigma='0'))
        self.expect_bad(lambda d:d['reference'].append(d['reference'][0].copy()))
        self.expect_bad(lambda d:d['exposures'].append(copy.deepcopy(d['exposures'][0])))
        self.expect_bad(lambda d:d['exposures'][0]['observations'].append(d['exposures'][0]['observations'][0].copy()))
        self.expect_bad(lambda d:d.update(exposures=[]))
        self.expect_bad(lambda d:d['exposures'][0].update(observations=[]))
        def too_many(d):
            d['exposures']=[{'id':f'e{k}','observations':[
                obs(f'o{i}','10','0','0.1','13') for i in range(1001)]} for k in range(2)]
        self.expect_bad(too_many)
        # 观测 ID 只在本曝光内唯一:跨曝光重复合法
        ok=self.minimal()
        ok['exposures'].append({'id':'e2','observations':[obs('o1','10','0','0.1','13')]})
        self.assertEqual(len(process_exposures(ok)['exposures']),2)

    def test_input_not_mutated_and_json(self):
        data=self.minimal()
        data['exposures'][0]['observations'].append(obs('o2','30','0','0.1','13'))
        snapshot=copy.deepcopy(data)
        result=process_exposures(data)
        self.assertEqual(data,snapshot)
        payload=json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False)
        self.assertEqual(json.loads(payload),json.loads(json.dumps(result)))

    def test_order_invariance(self):
        data=batch([
            ref('C1','50','0','15',True),ref('C2','60','0','15',True),ref('C3','70','0','15',True),
            ref('N1','80','0','15'),ref('N2','80','0','17')],[
            {'id':'e1','observations':[
                obs('k1','50','0','0.001','13'),obs('k2','60','0','0.001','13'),
                obs('k3','70','0','0.001','13'),obs('q','80','0','0.001','13')]},
            {'id':'e2','observations':[obs('w','10','0','0.001','9')]}])
        shuffled=copy.deepcopy(data)
        shuffled['reference']=list(reversed(shuffled['reference']))
        shuffled['exposures']=list(reversed(shuffled['exposures']))
        for e in shuffled['exposures']:e['observations']=list(reversed(e['observations']))
        self.assertEqual(process_exposures(data),process_exposures(shuffled))

class BruteForceCrosscheck(unittest.TestCase):
    """小规模随机用例对照暴力枚举全部匹配,验证最优方案刻画。"""
    def brute(self,data):
        threshold=Fraction(data['threshold']);tol=Fraction(data['anchor_tolerance'])
        refs={r['id']:(Fraction(r['ra']),Fraction(r['dec']),Fraction(r['magnitude']),r['calibrator']) for r in data['reference']}
        out={}
        for exp in data['exposures']:
            ob=[(o['id'],Fraction(o['ra']),Fraction(o['dec']),Fraction(o['sigma']),Fraction(o['instrumental_magnitude'])) for o in exp['observations']]
            cand=[]
            for oid,ra,dc,sigma,instr in ob:
                row=[]
                for rid,(rra,rdec,mag,cal) in refs.items():
                    dist=_distance(ra,dc,rra,rdec)
                    if dist<=threshold*sigma:row.append((rid,floor(dist*3600000+0.5)))
                cand.append(row)
            sols=[]
            def rec(i,used,assign,cost):
                if i==len(ob):sols.append((len(assign),cost,dict(assign)));return
                rec(i+1,used,assign,cost)
                for rid,c in cand[i]:
                    if rid not in used:
                        used.add(rid);assign[i]=rid
                        rec(i+1,used,assign,cost+c)
                        del assign[i];used.remove(rid)
            rec(0,set(),{},0)
            card=max(s[0] for s in sols)
            cost=min(s[1] for s in sols if s[0]==card)
            optima=[s[2] for s in sols if s[0]==card and s[1]==cost]
            forced=None
            for m in optima:
                items=set(m.items())
                forced=items if forced is None else forced&items
            anchors=sorted((i,rid) for i,rid in forced if refs[rid][3])
            offsets=[refs[rid][2]-ob[i][4] for i,rid in anchors]
            zero=None
            if len(offsets)>=3:
                so=sorted(offsets);q=len(so)
                med=so[q//2] if q%2 else (so[q//2-1]+so[q//2])/2
                sup=sum(1 for o in offsets if abs(o-med)<=tol)
                if sup>=3 and 2*sup>len(offsets):zero=med;status='calibrated'
                else:status='no_consensus'
            else:status='insufficient_anchors'
            final=optima
            if zero is not None:
                def photo(m):return sum(abs(ob[i][4]+zero-refs[rid][2]) for i,rid in m.items())
                best=min(photo(m) for m in optima)
                final=[m for m in optima if photo(m)==best]
            det={}
            for i,(oid,*_) in enumerate(ob):
                poss=tuple(sorted({m[i] for m in final if i in m}))
                det[oid]={'possible':poss,'may':any(i not in m for m in final),
                    'corrected':None if zero is None else _decimal_str(ob[i][4]+zero)}
            out[exp['id']]={'status':status,'zero':None if zero is None else _decimal_str(zero),
                'anchors':sorted((ob[i][0],rid,_decimal_str(refs[rid][2]-ob[i][4])) for i,rid in anchors),
                'det':det}
        return out

    def case(self,rng):
        def dec(x):return f'{x:.6f}'
        if rng.random()<0.5:
            # 孤立定标锚点对 + 争夺点:经常产生有效零点
            refs=[];ob=[]
            base_off=rng.choice(['2','2.5','-1'])
            for i in range(rng.randint(0,5)):
                mag=13+float(base_off)
                if rng.random()<0.25:mag+=rng.choice([5,-3,0.7])
                refs.append({'id':f'c{i}','ra':dec(10.0+i*30.0),'dec':'0.000000',
                    'magnitude':dec(mag),'calibrator':True})
                ob.append({'id':f'k{i}','ra':dec(10.0+i*30.0),'dec':'0.000000',
                    'sigma':'0.001000','instrumental_magnitude':'13.000000'})
            for i in range(rng.randint(0,3)):
                ra=5.0+i*7.0
                for j in range(rng.randint(1,3)):
                    refs.append({'id':f'n{i}_{j}','ra':dec(ra),'dec':'0.000000',
                        'magnitude':dec(rng.uniform(12,18)),'calibrator':False})
                for j in range(rng.randint(1,3)):
                    ob.append({'id':f'q{i}_{j}','ra':dec(ra),'dec':'0.000000',
                        'sigma':'0.001000','instrumental_magnitude':dec(rng.uniform(12,18))})
            if not ob:ob.append({'id':'z','ra':'200.000000','dec':'0.000000',
                'sigma':'0.001000','instrumental_magnitude':'13.000000'})
            if not refs:refs.append({'id':'r0','ra':'300.000000','dec':'0.000000',
                'magnitude':'15.000000','calibrator':False})
            return {'reference':refs,'threshold':'1',
                'anchor_tolerance':rng.choice(['0.3','0.5','1']),
                'exposures':[{'id':'e0','observations':ob}]}
        sites=[(rng.uniform(0,359),rng.uniform(-80,80)) for _ in range(rng.randint(3,6))]
        refs=[{'id':f'r{j}','ra':dec(rng.choice(sites)[0]),'dec':dec(rng.choice(sites)[1]),
            'magnitude':dec(rng.uniform(12,18)),'calibrator':rng.random()<0.6} for j in range(rng.randint(2,6))]
        exps=[]
        for e in range(rng.randint(1,2)):
            ob=[{'id':f'o{i}','ra':dec(rng.choice(sites)[0]),'dec':dec(rng.choice(sites)[1]),
                'sigma':dec(rng.uniform(0.05,2)),'instrumental_magnitude':dec(rng.uniform(12,18))}
                for i in range(rng.randint(1,5))]
            exps.append({'id':f'e{e}','observations':ob})
        return {'reference':refs,'threshold':rng.choice(['1','2']),
            'anchor_tolerance':rng.choice(['0.3','0.5','1']),'exposures':exps}

    def test_against_brute_force(self):
        rng=random.Random(20261124)
        seen=set()
        for trial in range(300):
            data=self.case(rng)
            got=process_exposures(data)
            want=self.brute(data)
            for exp in got['exposures']:
                be=want[exp['exposure']]
                cal=exp['calibration']
                seen.add(cal['status'])
                self.assertEqual(cal['status'],be['status'])
                self.assertEqual(cal['zero_point'],be['zero'])
                self.assertEqual(
                    sorted((a['observation'],a['reference'],a['offset']) for a in cal['anchors']),
                    be['anchors'])
                for d in exp['detections']:
                    bd=be['det'][d['observation']]
                    self.assertEqual(tuple(d['possible']),bd['possible'])
                    self.assertEqual(d['may_be_unassigned'],bd['may'])
                    self.assertEqual(d['corrected_magnitude'],bd['corrected'])
                    sel=bd['possible'][0] if len(bd['possible'])==1 and not bd['may'] else None
                    self.assertEqual(d['selected'],sel)
        self.assertEqual(seen,{'calibrated','no_consensus','insufficient_anchors'})

if __name__=='__main__':
    unittest.main()
