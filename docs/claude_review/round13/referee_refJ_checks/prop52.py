# Prop 5.2(2) check: the coordinate-span point u has consecutive gaps 2/n; G membership needs |u_x-u_y|>=e^{-eps_xy/2} PER PAIR.
# Construct pair data with max eps >= 2 ln(n/2) but where this u is NOT in G (small n illustration and general k).
import mpmath as mp, itertools
mp.mp.dps=30
for k in (4,16):
    n=2**k
    # u values
    def u(x):
        b=[(x>>i)&1 for i in range(k)]
        return -b[0]+sum(mp.mpf(2)**(-(i)) *b[i] for i in range(1,k))
    vals=sorted(u(x) for x in range(n)) if k<=12 else None
    if vals: print("k",k,"min gap",min(vals[i+1]-vals[i] for i in range(n-1)),"2/n",mp.mpf(2)/n, "max|u|",max(abs(v) for v in vals))
    # pair data: all eps = E_big except eps for one pair of adjacent-valued points = 0
    Ebig=2*mp.log(n/2)+1
    print(f"k={k}: max eps={mp.nstr(Ebig,8)} >= 2ln(n/2)={mp.nstr(2*mp.log(n/2),8)}; one pair with eps=0 needs gap >= 1; u's adjacent pair (x=0, x with value 2/n) has gap 2/n<1 -> u notin G")
    # every L-point (span of coordinate columns, base 0) : values subset sums; claim that some u in L lies in G for that data
print("Average constraint: sum eps <= N2*abar => min eps <= abar. Corrected (2) needs min eps>=2ln(n/2):")
for k in range(16,21):
    n=2**k; d=n-1
    Wmin=3*n*mp.log(n)-14*n
    print(k, "abar_ne",mp.nstr(Wmin/(2*d),6),"2ln(n/2)",mp.nstr(2*mp.log(n/2),6), "W needed for abar>=2ln(n/2): W/Wmin >=", mp.nstr(2*mp.log(n/2)*2*d/Wmin,5))
