import mpmath as mp
mp.mp.dps=50
def par(n,r):
    d=n-1; w=mp.sqrt(8*d*mp.log(4*r)); j=int(mp.floor(w/2)); m=int(mp.ceil(4*w)); return d,w,j,m
rows=[(2900,2900),(4096,4095),(10**4,10**4),(10**4,10**8),(2**16,2**16-1),(10**5,10**5),(10**6,10**6),(10**8,10**8),(10**9,10**9)]
def trinom(bmax):
    # central trinomial coefficients T_b exact: (b+1)T_{b+1} = (2b+1)T_b + 3b T_{b-1}
    T=[1,1]
    for b in range(1,bmax): T.append(((2*b+1)*T[b]+3*b*T[b-1])//(b+1))
    return T
T=trinom(5000)
# check recurrence correctness on small values
assert T[:6]==[1,1,3,7,19,51]
print("Theorem 3.1 table: n r m b log10 mu{Q!=0}; and p_b>=0.4/sqrt(b) check")
for n,r in rows:
    d,w,j,m=par(n,r)
    if n>10**6: continue
    b=d//m; pb=mp.mpf(T[b])/mp.mpf(3)**b
    print(n,r,m,b, mp.nstr(m*mp.log10(1-pb),8))
bad=[b for b in range(2,3001) if mp.mpf(T[b])/mp.mpf(3)**b < 0.4/mp.sqrt(b)]
print("p_b<0.4/sqrt b for b in 2..3000:",bad[:5])
print("\nCW proxy (C=1, sigma=sqrt(m!) as author; also sigma=(2/3)^{m/2}sqrt(m!)), W=Wmin")
for n,r in rows[:7]:
    d,w,j,m=par(n,r); W=3*n*mp.log(n)-14*n
    L0=W/4-mp.log(mp.mpf(d)/(m+1))
    lnsig=mp.loggamma(m+1)/2
    for lab,ls in [('sqrt(m!)',lnsig),('(2/3)^(m/2)',lnsig+m/2*mp.log(mp.mpf(2)/3))]:
        val=mp.log(m)+(m*mp.log(d)-L0-ls)/m
        print(n,r,m,lab,mp.nstr(val/mp.log(10),4))
print("\nLedger: n m L0 Lc L0-Lc-m ln n  L0-Lc-L16 (W=Wmin,C=1)")
for n,r in rows:
    if (n,r)==(10**4,10**8) or n==10**8: continue
    d,w,j,m=par(n,r); W=3*n*mp.log(n)-14*n; abar=W/(2*d)
    L0=W/4-mp.log(mp.mpf(d)/(m+1)); Lc=m*mp.log(2*d*mp.e/m)
    print(n,m,mp.nstr(L0,6),mp.nstr(Lc,6),mp.nstr(L0-Lc-m*mp.log(n),6),mp.nstr(L0-Lc-m*(8*abar+mp.log(2)),6))
print("\nProp 4.2 partition q_max at theta=1/2 (equal c, N=floor(d/q))")
for n,r in rows:
    d,w,j,m=par(n,r); best=0
    for q in range(1,400):
        N=d//q; c=(m-1)//(2*q)
        tail=2*mp.e**(-3*mp.mpf(c+1)**2/(4*N))
        if tail<1 and (1-tail)**q>=mp.mpf(1)/2: best=q
    print(n,r,m,best)
