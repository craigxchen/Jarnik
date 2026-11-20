import itertools
exec(open('dmin.py').read().split('for P in map')[0])
# brute force small P
for P in (10,20,30):
    ps=[p for p in primes(P,2*P) if p%4==1]; q=1
    for p in ps:q*=p
    from math import isqrt
    pts=[(x,y) for x in range(-isqrt(q),isqrt(q)+1) for y in [isqrt(q-x*x),-isqrt(q-x*x)] if x*x+y*y==q]
    pts=list(set(pts))
    bf=min((a[0]-b[0])**2+(a[1]-b[1])**2 for a,b in itertools.combinations(pts,2))
    print("bruteforce P",P,"npts",len(pts),"dmin2",bf); run(P,4,1)
for P in (200,250,300): run(P,8,5)
# lift count: for odd l, a>=1, q with l not dividing q: number of classes mod l^{a+1} of norm q lifting a fixed class mod l^a
for l in (3,5,7,11):
  for a in (1,2):
    for q in (13,17,29*37):
        if q%l==0: continue
        M=l**(a+1); m=l**a
        cls={}
        for x in range(M):
            for y in range(M):
                if (x*x+y*y-q)%M==0: cls.setdefault((x%m,y%m),0); cls[(x%m,y%m)]+=1
        assert set(cls.values())=={l}, (l,a,q,set(cls.values()))
    # a=0 count
    for q in (13,17):
        if q%l: print("l",l,"q",q,"#sol mod l",sum(1 for x in range(l) for y in range(l) if (x*x+y*y-q)%l==0))
print("lift counts all equal l: OK")
