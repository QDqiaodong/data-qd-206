from fractions import Fraction
from .checks import fields,array,ident,unique,number

def process(data):
 fields(data,('reference','observations','threshold'));threshold=number(data['threshold']); reference={}
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
 matches=[]
 for o in observations:
  candidates=[]
  for rid,(ra,dec,mag) in reference.items():
   dist=((o['ra']-ra)**2+(o['dec']-dec)**2)**Fraction(1,2)
   if dist<=threshold*o['sigma']: candidates.append((dist,rid,mag))
  candidates.sort(key=lambda z:(z[0],z[1]))
  matches.append({'observation':o['id'],'candidates':[{'reference':rid,'distance':str(d),'magnitude':str(m)} for d,rid,m in candidates],'selected':candidates[0][1] if len(candidates)==1 else None,'status':'matched' if len(candidates)==1 else ('ambiguous' if candidates else 'unmatched')})
 return {'matches':matches,'summary':{'reference':len(reference),'observations':len(observations),'matched':sum(m['status']=='matched' for m in matches),'ambiguous':sum(m['status']=='ambiguous' for m in matches),'unmatched':sum(m['status']=='unmatched' for m in matches)}}
