from math import sin,cos,atan2,sqrt,radians,degrees
from .checks import fields,array,ident,unique,number

def process(data):
 fields(data,('reference','observations','threshold'));threshold=number(data['threshold']); reference={}
 if threshold<=0:raise ValueError('positive threshold required')
 for x in array(data['reference'],1,1000):
  fields(x,('id','ra','dec','magnitude'));k=ident(x['id']);
  if k in reference:raise ValueError('duplicate reference')
  reference[k]=(number(x['ra']),number(x['dec']),number(x['magnitude']))
 observations=[]
 for x in array(data['observations'],1,2000):
  fields(x,('id','ra','dec','sigma'));k=ident(x['id']);
  if any(o['id']==k for o in observations):raise ValueError('duplicate observation')
  observations.append({'id':k,'ra':number(x['ra']),'dec':number(x['dec']),'sigma':number(x['sigma'])})
  if observations[-1]['sigma']<=0:raise ValueError('sigma must be positive')
 for ra,dec,_ in reference.values():
  if not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid sky coordinates')
 for o in observations:
  if not 0<=o['ra']<360 or not -90<=o['dec']<=90:raise ValueError('invalid sky coordinates')
 matches=[]
 for o in observations:
  candidates=[]
  for rid,(ra,dec,mag) in reference.items():
   a,b=radians(float(o['dec'])),radians(float(dec))
   h=sin((a-b)/2)**2+cos(a)*cos(b)*sin(radians(float(o['ra']-ra))/2)**2
   h=min(1.0,max(0.0,h));dist=degrees(2*atan2(sqrt(h),sqrt(1-h)))
   if dist<=threshold*o['sigma']: candidates.append((dist,rid,mag))
  candidates.sort(key=lambda z:(z[0],z[1]))
  matches.append({'observation':o['id'],'candidates':[{'reference':rid,'distance':str(d),'magnitude':str(m)} for d,rid,m in candidates],'selected':candidates[0][1] if len(candidates)==1 else None,'status':'matched' if len(candidates)==1 else ('ambiguous' if candidates else 'unmatched')})
 return {'matches':matches,'summary':{'reference':len(reference),'observations':len(observations),'matched':sum(m['status']=='matched' for m in matches),'ambiguous':sum(m['status']=='ambiguous' for m in matches),'unmatched':sum(m['status']=='unmatched' for m in matches)}}
