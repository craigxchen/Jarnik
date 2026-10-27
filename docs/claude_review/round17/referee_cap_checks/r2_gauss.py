# Referee pass 2: exact check of Lemma 2.1 / Cor 2.2 in Z[zeta_q] (own code). Patterns built literally as in
# round16/referee_barrier_checks/sdc_sym.py: P_a = {a} u {x: chi(a-x)=1}, P'_a = {a} u {x: chi(a-x)=-1}.
import random, sys
from flint import fmpz_poly, fmpq_mat
from sympy import isprime
out=[]
def P(*a):
    s=" ".join(map(str,a)); print(s); out.append(s); sys.stdout.flush()
def chi(a,q):
    a%=q
    if a==0: return 0
    return 1 if pow(a,(q-1)//2,q)==1 else -1
random.seed(20261027)
allok=True
for q in [p for p in range(7,200) if isprime(p) and p%4==3]:
    Phi = fmpz_poly([1]*q)
    def z(e): 
        c=[0]*q; c[e%q]=1; return fmpz_poly(c)
    g = sum((chi(x,q)*z(x) for x in range(1,q)), fmpz_poly([0]))
    ok = ((g*g) % Phi) == fmpz_poly([-q])
    Pa  = [[a]+[x for x in range(q) if chi(a-x,q)==1] for a in range(q)]
    Ppa = [[a]+[x for x in range(q) if chi(a-x,q)==-1] for a in range(q)]
    # set identities
    ok &= all(set(Pa[a])|set(Ppa[a])==set(range(q)) and set(Pa[a])&set(Ppa[a])=={a} and len(Pa[a])==(q+1)//2 for a in range(q))
    ok &= all(set(Ppa[a])==set([a]+[(a+y)%q for y in range(1,q) if chi(y,q)==1]) for a in range(q))
    for t in range(15):
        k=[random.randint(-5,5) for _ in range(q)]
        kh = sum((k[x]*z(x) for x in range(q)), fmpz_poly([0]))
        S=[sum(k[x] for x in Ppa[a]) for a in range(q)]; T=[sum(k[x] for x in Pa[a]) for a in range(q)]
        lhs1 = ((fmpz_poly([1])-g)*kh) % Phi; rhs1 = (2*sum((S[a]*z(a) for a in range(q)), fmpz_poly([0]))) % Phi
        lhs2 = ((fmpz_poly([1])+g)*kh) % Phi; rhs2 = (2*sum((T[a]*z(a) for a in range(q)), fmpz_poly([0]))) % Phi
        ok &= (lhs1==rhs1) and (lhs2==rhs2)
        ok &= all(S[a]+T[a]==sum(k)+k[a] for a in range(q))
    # W invertible over Q (rows (Wk)_a = k(P'_a))
    if q<=103:
        Wm = fmpq_mat(q,q,[1 if x in Ppa[a] else 0 for a in range(q) for x in range(q)])
        ok &= Wm.rank()==q
    allok &= ok
    P(f"q={q}: g^2=-q, (2.1) both, S_a+T_a=K+k_a, patterns, W invertible(q<=103): {ok}")
P("ALL OK:", allok)
open(__file__.replace('.py','_out.txt'),'w').write("\n".join(out)+"\n")
