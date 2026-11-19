from cap2 import run
for L in list(range(13,37,2))+[35,36,37]:
    try:
        r,_,_=run(3,2,[(1,0),(0,1)],[1,1],L,verbose=False); print(L,r['D'],r['max_pair'],r['critical'],r['NQ=R^2'],r['norms_equal'])
    except AssertionError as e: print(L,'assert',e)
for fs in ([(1,0),(0,1),(1,1)],[(1,0),(0,1),(1,1),(1,2)]):
  for L in (61,200):
    r,_,_=run(3,2,fs,[1]*len(fs),L,verbose=False); print(len(fs),L,r['g'],r['D'],r['pair_degs'],r['critical'],r['norms_equal'],r['distinct'])
