# Tangent space check: span(1_I, columns of K_{S^c}^T restricted to I) == {x in Q^I : K_{S,I} x in Q*1}
# (exact rank computations over Q with sympy), plus the projection-after-normalisation claim.
import random
from sympy import Matrix
def paley(q):
    sq={x*x%q for x in range(1,q)}
    chi=[0]+[1 if x in sq else -1 for x in range(1,q)]
    return [[chi[(j-i)%q] for j in range(q)] for i in range(q)]
for q,rho,seed in ((43,0.9,1),(59,0.95,2),(103,0.99,3)):
    random.seed(seed)
    K=paley(q); A=[i for i in range(q) if K[i][0]==1]
    S=[i for i in A if random.random()<rho]; Sc=[i for i in range(q) if i not in S]
    I=list(range(1,q))
    KSI=Matrix([[K[i][j] for j in I] for i in S])
    # T = {x : (K_S x)_i - (K_S x)_{i0} = 0}
    Dm=Matrix([[KSI[a,c]-KSI[0,c] for c in range(len(I))] for a in range(1,len(S))])
    dimT=len(I)-Dm.rank()
    V=Matrix([[1]*len(I)]+[[K[i][j] for j in I] for i in Sc]).T
    rV=V.rank()
    inT = (Dm*V).is_zero_matrix
    print(f"q={q} |S|={len(S)} dim T={dimT} rank span(1,K_Sc^T)={rV} contained={inT} equal={inT and rV==dimT}")
