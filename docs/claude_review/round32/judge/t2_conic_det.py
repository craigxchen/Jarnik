# [J] Test 2 (slopes analyst, Section 3.2): for M=2d+1 lattice points on |z|^2=N,
#   det[ z_j^a (0<=a<=d), zbar_j^b (1<=b<=d) ] = +- N^{d(d+1)/2} V(z) / prod_j z_j^d   (V = prod_{i<j}(z_j - z_i))
# hence |V|^2 >= N^{d^2} = N^{(M-1)^2/4}; and check this never exceeds what the cut identity gives
# (parity part: prod|z_i-z_j|^8 >= N^{M(M-2)+1} * prod|c_ij|^8 with |c_ij|>=1).
import sympy as sp, math, random, itertools
I=sp.I
def pts(N):
    out=set(); a=0
    while a*a<=N:
        b2=N-a*a; b=math.isqrt(b2)
        if b*b==b2:
            for x,y in [(a,b),(-a,b),(a,-b),(-a,-b),(b,a),(-b,a),(b,-a),(-b,-a)]: out.add((x,y))
        a+=1
    return sorted(out,key=lambda p: math.atan2(p[1],p[0]))
random.seed(2)
ok=0; tot=0; ratios=[]
for N in [5*13*17*29, 5*13*17*29*37, 25*13*17*41]:
    P=pts(N)
    for d in (1,2,3):
        M=2*d+1
        for trial in range(6):
            if trial%2==0:
                st=random.randrange(len(P)); S=[P[(st+i)%len(P)] for i in range(M)]
            else: S=random.sample(P,M)
            z=[sp.Integer(x)+I*sp.Integer(y) for x,y in S]
            zb=[sp.Integer(x)-I*sp.Integer(y) for x,y in S]
            rows=[[zj**a for a in range(d+1)]+[zbj**b for b in range(1,d+1)] for zj,zbj in zip(z,zb)]
            D=sp.expand(sp.Matrix(rows).det(method='bareiss'))
            V=sp.expand(sp.prod([z[j]-z[i] for i in range(M) for j in range(i+1,M)]))
            rhs=sp.expand(sp.Integer(N)**(d*(d+1)//2)*V)
            lhs=sp.expand(D*sp.prod([zj**d for zj in z]))
            tot+=1
            if sp.simplify(lhs-rhs)==0 or sp.simplify(lhs+rhs)==0: ok+=1
            else: print("MISMATCH",N,d,S)
            absV2=sp.Abs(V)**2
            # log-margin of |V|^2 over N^{d^2}
            ratios.append(float(sp.log(absV2,10)-d*d*sp.log(N,10)))
print("identity checked",ok,"/",tot)
print("min log10(|V|^2/N^{d^2}) =",round(min(ratios),3),"(>=0 required)")
