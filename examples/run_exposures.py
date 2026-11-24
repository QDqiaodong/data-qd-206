"""离线示例:多次曝光整批关联与测光零点校正,输出写入 run/exposures.json。"""
from star_crossmatch import process_exposures
from star_crossmatch.io import publish

DATA={
 'reference':[
  {'id':'C1','ra':'50','dec':'0','magnitude':'15','calibrator':True},
  {'id':'C2','ra':'60','dec':'0','magnitude':'15','calibrator':True},
  {'id':'C3','ra':'70','dec':'0','magnitude':'15','calibrator':True},
  {'id':'N1','ra':'80','dec':'0','magnitude':'15','calibrator':False},
  {'id':'N2','ra':'80','dec':'0','magnitude':'17','calibrator':False}],
 'threshold':'1',
 'anchor_tolerance':'0.5',
 'exposures':[
  {'id':'night1','observations':[
   {'id':'k1','ra':'50','dec':'0','sigma':'0.001','instrumental_magnitude':'13'},
   {'id':'k2','ra':'60','dec':'0','sigma':'0.001','instrumental_magnitude':'13'},
   {'id':'k3','ra':'70','dec':'0','sigma':'0.001','instrumental_magnitude':'13'},
   {'id':'q','ra':'80','dec':'0','sigma':'0.001','instrumental_magnitude':'13'},
   {'id':'z','ra':'85','dec':'0','sigma':'0.001','instrumental_magnitude':'13'}]},
  {'id':'night2','observations':[
   {'id':'k1','ra':'50','dec':'0','sigma':'0.001','instrumental_magnitude':'13.4'},
   {'id':'k2','ra':'60','dec':'0','sigma':'0.001','instrumental_magnitude':'13.4'},
   {'id':'k3','ra':'70','dec':'0','sigma':'0.001','instrumental_magnitude':'13.4'},
   {'id':'w','ra':'80','dec':'0','sigma':'0.001','instrumental_magnitude':'13.4'}]}]}

if __name__=='__main__':
    result=process_exposures(DATA)
    for exp in result['exposures']:
        cal=exp['calibration']
        print(exp['exposure'],cal['status'],'zero_point=',cal['zero_point'])
        for d in exp['detections']:
            print(' ',d['observation'],d['status'],'selected=',d['selected'],
                'possible=',d['possible'],'corrected=',d['corrected_magnitude'])
    publish('run/exposures.json',result)
    print('written run/exposures.json')
