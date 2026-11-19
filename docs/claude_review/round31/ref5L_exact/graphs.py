# Which isotropic edge sets S in K5 satisfy the peeling condition at level s in {1,2,3}:
# exist distinct a,j,k with deg_S(a)>=s, deg_{S-a}(j)>=s-1, deg_{S-a-j}(k)>=s-2.
from itertools import combinations, permutations
V=range(5); E=list(combinations(V,2))
def deg(S,v,removed): return sum(1 for e in S if v in e and not (set(e)&removed))
def peel(S,s):
    for a,j,k in permutations(V,3):
        if deg(S,a,set())>=s and deg(S,j,{a})>=s-1 and deg(S,k,{a,j})>=s-2: return True
    return False
for s in (1,2,3):
    row=[]
    for m in range(0,11):
        subs=list(combinations(E,m)); good=sum(peel(S,s) for S in subs)
        row.append(f"{m}:{good}/{len(subs)}")
    print('s=',s,' '.join(row))
# list the failing 6- and 7-edge graphs at s=3 up to isomorphism (degree sequences)
for m in (6,7,8):
    bad=[S for S in combinations(E,m) if not peel(S,3)]
    print(m,'edges failing s=3:',len(bad), sorted({tuple(sorted(deg(S,v,set()) for v in V)) for S in bad}))
S4=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]; print('K4 passes s=3:',peel(S4,3))
star=[(0,1),(0,2),(0,3),(0,4)]; print('4-star passes s=1,2,3:',[peel(star,s) for s in (1,2,3)])
five=[(0,1),(0,2),(0,3),(0,4),(1,2)]; print('five-block graph passes s=1,2,3:',[peel(five,s) for s in (1,2,3)])
