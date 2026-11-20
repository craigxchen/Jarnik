# Part (1): Paley core checks; exact integers
from itertools import combinations
def chi(a,m):
    a%=m
    if a==0: return 0
    return 1 if pow(a,(m-1)//2,m)==1 else -1
def S_of(m):
    return [[chi(v-i,m)-(1 if i==v else 0) for v in range(m)] for i in range(m)]
for m in [7,11,19,23,31,43,47]:
    S=S_of(m)
    ok=all(abs(S[i][v])==1 for i in range(m) for v in range(m))
    G=[[sum(S[i][t]*S[j][t] for t in range(m)) for j in range(m)] for i in range(m)]
    Gc=[[sum(S[t][i]*S[t][j] for t in range(m)) for j in range(m)] for i in range(m)]
    gram=all(G[i][j]==((m+1) if i==j else 0)-1 for i in range(m) for j in range(m))
    gramc=all(Gc[i][j]==((m+1) if i==j else 0)-1 for i in range(m) for j in range(m))
    rs=all(sum(S[i])==-1 for i in range(m)); cs=all(sum(S[i][v] for i in range(m))==-1 for v in range(m))
    sep=set(sum(S[i][v]!=S[j][v] for v in range(m)) for i,j in combinations(range(m),2))
    const=set()
    for i,j,k in combinations(range(m),3):
        const.add(sum(S[i][v]==S[j][v]==S[k][v] for v in range(m)))
    print(m, m%8, "pm1",ok,"SS^t",gram,"S^tS",gramc,"S1=-1",rs,"1S=-1",cs,"sep",sep,"(m+1)/2=",(m+1)//2,"const cols per triple",const,"(m-3)/4=",(m-3)/4)
