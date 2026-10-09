"""Structured demand of design D_A on additive-quadruple characters c = e_0 - e_1 - e_x + e_(x+1).
At a small prime q (< 4M), level 1, kappa_y = sigma_q(y) + 2 (y mod p'(q)); the mod-p' parts cancel
unless p'(q) | x+1 (a carry), so c collides (unit 1) whenever its sigma-part vanishes -- about half the
small primes for random parity colourings.  Reports log m_struct(c)/M for random sigma (parities iid),
compared with the sharp.md bound log m_struct <= 11 M and theta(4M) ~ 4M."""
import math, random
rng=random.Random(8)
def primes_upto(n):
    s=bytearray([1])*(n+1); s[0:2]=b"\x00\x00"
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P=primes_upto(2**20)
pp={3:2,5:2,7:3}
for k in range(3,19):
    B=[q for q in P if 2**k<=q<2**(k+1) and q>7]; T=[p for p in P if 2**(k-2)<=p<2**(k-1)]
    for i,q in enumerate(B): pp[q]=T[(i*len(T))//len(B)]
for M in (200,1000,5000,20000):
    qs=[q for q in P if 2<q<4*M]
    theta4M=sum(math.log(q) for q in qs)
    best=0; tot=0; cnt=0
    for trial in range(20):
        x=rng.randrange(2,M-1)
        rows=[0,1,x,x+1]; c=[1,-1,-1,1]
        L=0.0
        for q in qs:
            sig=[rng.randint(0,1) for _ in rows]   # parity colouring sigma_q on the 4 rows (random model)
            p=pp[q]; m=q-1 if q%4==1 else q+1
            val=sum(ci*(si+2*(y%p)) for ci,si,y in zip(c,sig,rows))
            if val % m == 0:      # unit 1; other units impossible for |val| small unless m tiny
                L+=math.log(q)
        best=max(best,L); tot+=L; cnt+=1
    print(f"M={M:6d}: additive quadruples, level-1 structured log m: mean {tot/cnt/M:.2f} M, max {best/M:.2f} M  (theta(4M) = {theta4M/M:.2f} M; bound 11 M)")
