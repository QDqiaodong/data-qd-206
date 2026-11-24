from star_crossmatch import process
from star_crossmatch.io import load,publish
if __name__=='__main__':
 r=process(load('fixtures/input.json'));publish('run/sample.json',r);print(r['summary'])
