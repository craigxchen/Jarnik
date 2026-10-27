# Exact Walsh-Fourier coefficients of Phi(eps)=|f.eps| for ternary sum-zero f with n nonzero entries (n=4,6,8),
# and the constants E_n=Phi^(0), A_n=sum_{I!=0}|Phi^(I)|(|I|-1), B_n=sum_{I!=0}|Phi^(I)| for the Weil bound.
import itertools
from fractions import Fraction as Fr
for n in [4,6,8,10]:
    f=[1]*(n//2)+[-1]*(n//2)
    eps=list(itertools.product([1,-1],repeat=n))
    Phi={e:abs(sum(a*b for a,b in zip(f,e))) for e in eps}
    E=Fr(sum(Phi.values()),2**n); A=Fr(0); B=Fr(0); byk={}
    for k in range(1,n+1):
        for I in itertools.combinations(range(n),k):
            co=Fr(sum(Phi[e]*eval('*'.join(str(e[i]) for i in I)) for e in eps),2**n)
            if co: 
                A+=abs(co)*(k-1); B+=abs(co); byk.setdefault(k,set()).add(abs(co))
    print(f"n={n} E_n={E}={float(E):.6f} A_n={A}={float(A):.4f} B_n={B}={float(B):.4f} |coef| by |I|={ {k:sorted(v) for k,v in byk.items()} }",flush=True)
