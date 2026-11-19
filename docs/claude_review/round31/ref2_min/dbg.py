exec(open('flats.py').read().split('level = {0')[0])
ks={}
for r in range(70):
    nb=add([],GEN[r]); k=key(nb); ks.setdefault(k,[]).append(r)
print([v for v in ks.values() if len(v)>1]); print(GEN[0], add([],GEN[0]))
