import json,sys
for f in sys.argv[1:]:
    print('==',f)
    for line in open(f):
        d=json.loads(line); print(d['q0'],'PU',d['PU'],'SL',d['SL0'],'P0',d['P0'],'idx',d['index'],'ANGout',d.get('ANG_out'),'A1ANG',d.get('A1_ANG'),d.get('profiles',''))
