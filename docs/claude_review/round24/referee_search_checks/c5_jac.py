import mpmath as mp, itertools, sys
mp.mp.dps=80
sys.argv=[sys.argv[0]]+sys.argv[1:]
from c4_construct import curve, places, ray, funcs
def esym(vals):
    e=[mp.mpc(1)]
    for x in vals:
        e=[ (e[k] if k<len(e) else 0) - (x*e[k-1] if k>=1 else 0) for k in range(len(e)+1)]
    return [ (-1)**k*e[k] for k in range(len(e))]  # e_k
def system(rho):  # rho: dict T->complex, T subset of [4]
    out=[]
    for j in (1,2,3,4):
        Rj=[rho[T] for T in rho if j in T]
        e=esym(Rj)
        for k in (1,2,3,5,6,7,8): out.append(mp.im(e[k]))
    return out
for m in ['-1.5','-2.2']:
    C=curve(mp.mpf(m)); pl=places(C)
    rho={}
    for w,v in pl:
        A=ray(v)
        if 5 not in A: rho[A]=w
    assert len(rho)==15
    Ts=sorted(rho)
    F0=system(rho); print('m',m,'max |Im e_k| (k!=4):',mp.nstr(max(abs(x) for x in F0),3))
    kap4=[mp.im(esym([rho[T] for T in rho if j in T])[4]) for j in (1,2,3,4)]
    print('  Im e_4 (kappa_j):',[mp.nstr(x,5) for x in kap4])
    h=mp.mpf(10)**-35
    J=mp.matrix(28,30)
    for i,T in enumerate(Ts):
        for part,dz in ((0,h),(1,1j*h)):
            r2=dict(rho); r2[T]=rho[T]+dz
            r1=dict(rho); r1[T]=rho[T]-dz
            Fp=system(r2); Fm=system(r1)
            for a in range(28): J[a,2*i+part]=(Fp[a]-Fm[a])/(2*h)
    U,S,V=mp.svd_r(J)
    sv=sorted([S[i] for i in range(28)])
    print('  singular values (smallest 4):',[mp.nstr(x,4) for x in sv[:4]],' largest',mp.nstr(sv[-1],4))
    # family tangent: d rho/dm (match places by ray) and gauge tangent rho
    dm=mp.mpf(10)**-30
    C2=curve(mp.mpf(m)+dm); pl2=places(C2); rho2={}
    for w,v in pl2:
        A=ray(v)
        if 5 not in A: rho2[A]=w
    tang=mp.matrix(30,1); gauge=mp.matrix(30,1)
    for i,T in enumerate(Ts):
        d=(rho2[T]-rho[T])/dm; tang[2*i]=mp.re(d); tang[2*i+1]=mp.im(d)
        gauge[2*i]=mp.re(rho[T]); gauge[2*i+1]=mp.im(rho[T])
    print('  |J*tangent|/|tangent| =',mp.nstr(mp.norm(J*tang)/mp.norm(tang),3),' |J*gauge|/|gauge| =',mp.nstr(mp.norm(J*gauge)/mp.norm(gauge),3))
    # independence of tangent and gauge
    M2=mp.matrix(30,2)
    for i in range(30): M2[i,0]=tang[i]/mp.norm(tang); M2[i,1]=gauge[i]/mp.norm(gauge)
    print('  sing values [tangent,gauge]:',[mp.nstr(x,4) for x in mp.svd_r(M2,compute_uv=False)])
