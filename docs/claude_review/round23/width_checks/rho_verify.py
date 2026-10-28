# Exact DP for the recursive peeling bound rho_n(a,b) (n<=7) and check of the closed forms.
import sys
def peel(vals):
    v=sorted(vals,reverse=True); return max(x+i for i,x in enumerate(v))
def table(nmax,L):
    T={2:[[a+b for b in range(L+1)] for a in range(L+1)]}
    for n in range(3,nmax+1):
        P=T[n-1]
        T[n]=[[peel([P[a-c][b] for c in range(a+1)]+[P[a][b-c] for c in range(1,b+1)]) for b in range(L+1)] for a in range(L+1)]
    return T
L=int(sys.argv[1]); T=table(7,L)
F={3:lambda a,b:a+b+min(a,b),
   4:lambda a,b:a+b+min((a+b)//2,2*min(a,b)),
   5:lambda a,b:a+b+min(min(a,b)+max(a,b)//3,3*min(a,b))}
for n in (3,4,5):
    bad=[(a,b) for a in range(L+1) for b in range(L+1) if T[n][a][b]!=F[n](a,b)]
    print('n',n,'closed form matches rho on [0,%d]^2:'%L, not bad, bad[:5])
print('rho5(m,m)==floor(10m/3) for m<=%d:'%L, all(T[5][m][m]==(10*m)//3 for m in range(L+1)))
for n in (6,7):
    print('n',n,'rho(m,m):',[T[n][m][m] for m in range(1,min(L,30)+1)])
