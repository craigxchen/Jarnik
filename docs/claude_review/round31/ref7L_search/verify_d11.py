# independent check of the D=11, h=6 real five-point configuration from 40-digit roots
import mpmath as mp, re, itertools
mp.mp.dps=80
rows=[l.split(None,2) for l in open('d11_type12_roots.txt').read().splitlines()[1:]]
cols=[];rho=[]
for k,c,z in rows:
    cols.append([int(ch) for ch in c]); rho.append(mp.mpc(complex(0)) + mp.mpmathify(z.strip('()').replace(' ','').replace('j','j')))
def poly(roots):
    c=[mp.mpc(1)]
    for x in roots: c=[(c[j] if j<len(c) else 0)-x*(c[j-1] if j>=1 else 0) for j in range(len(c)+1)]
    return c
Z=[]
for i in range(5):
    roots=[rho[k] if cols[k][i] else mp.conj(rho[k]) for k in range(len(cols))]
    Z.append(poly(roots))
D=len(Z[0])-1
print("D",D)
for i,j in itertools.combinations(range(5),2):
    d=[abs(Z[i][k]-Z[j][k]) for k in range(D+1)]
    print(i,j,"|coef t^10..t^0| of z_i-z_j:",[mp.nstr(x,2) for x in d[1:]])
# norm check z*conj(z) equal
def mul(a,b): 
    r=[0]*(len(a)+len(b)-1)
    for x,u in enumerate(a):
        for y,w in enumerate(b): r[x+y]+=u*w
    return r
N=[mul(z,[mp.conj(c) for c in z]) for z in Z]
print("max norm discrepancy",mp.nstr(max(abs(N[i][k]-N[0][k]) for i in range(5) for k in range(2*D+1)),3))
print("min |Im rho|",mp.nstr(min(abs(mp.im(x)) for x in rho),5))
