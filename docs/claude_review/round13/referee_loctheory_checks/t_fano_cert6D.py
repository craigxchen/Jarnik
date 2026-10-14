# CERTIFY: J(Y_3(2))_7 = grid + lines + quads + width powers + sum over GL(3,2)-images of the 6-direction family with
# kernel z = e_a + e_b.  Generators of the image ideal computed EXACTLY over Q (fmpz nullspace); sandwich via mod-p rank.
import itertools, sys, time, pickle
from rlib import *
p=2147483647
D=int(sys.argv[1])
Y,h,dJ6,dJ7=pickle.load(open('fano_Y.pkl','rb'))
hD=h_func_exact(Y,D,[5]*7); hDm=h_func_exact(Y,D-1,[5]*7)
dJD=len(monos_eq(7,D,[5]*7))-(hD-hDm)
print('exact target dim (J/grid)_D',dJD)
gens=pickle.load(open('fano_gens.pkl','rb'))
cap=[5]*7; D=int(sys.argv[1])
cols=monos_eq(7,D,cap); cidx={c:i for i,c in enumerate(cols)}
avec=[a for a in itertools.product((0,1),repeat=3) if any(a)]
def mat_apply(Mx,a): return tuple(sum(Mx[r][c]*a[c] for c in range(3))%2 for r in range(3))
perms=[]
for Mx in itertools.product(itertools.product((0,1),repeat=3),repeat=3):
    imgs=[mat_apply(Mx,a) for a in avec]
    if len(set(imgs))==7 and all(any(x) for x in imgs): perms.append([avec.index(x) for x in imgs])
z=[1,1,0,0,0,0,0]
K,nul=flint.fmpz_mat([z]).nullspace()
F=[[int(K[r,c]) for r in range(7)] for c in range(nul)]
print("directions (rows):",F)
img=sorted(set(tuple(sum(a*b for a,b in zip(v,y)) for v in F) for y in Y))
print("|img|",len(img)); sys.stdout.flush()
rho_up=next(e for e in range(40) if h_func_modp(img,e)==len(img))
print("rho(img) <=",rho_up,"(mod-p upper bound)"); sys.stdout.flush()
T0=time.time(); newg={}
for e in range(1,D+1):
    sp=J_space(img,e)      # EXACT over Q
    if sp: newg[e]=[pullback(g,F,7) for g in sp]
    print("exact J(img)_%d dim"%e,len(sp),"%.0fs"%(time.time()-T0)); sys.stdout.flush()
base=ideal_degree_rows(newg,7,D,cap)
def rows_to_list(rows):
    out=[]
    for r in rows:
        v=[0]*len(cols)
        for k,x in r.items(): v[cidx[k]]=x%p
        out.append(v)
    return out
def basis_of(lst):
    M=flint.nmod_mat(lst,p); R,rk=M.rref()
    return [[int(R[i,c]) for c in range(len(cols))] for i in range(rk)]
B=basis_of(rows_to_list(ideal_degree_rows(gens,7,D,cap)))
print("start",len(B))
done=set()
for pm in perms:
    key=tuple(sorted([pm[0],pm[1]]))   # image of kernel support {0,1}
    rows=[]
    for r in base:
        nr={}
        for ex,v in r.items():
            ne=[0]*7
            for i,a in enumerate(ex): ne[pm[i]]=a
            nr[tuple(ne)]=v
        rows.append(nr)
    B=basis_of(B+rows_to_list(rows))
    if len(B)>=dJD: break
print("mod-p rank of exact generators:",len(B),"target (exact)",dJD, "CERTIFIED EQUAL" if len(B)==dJD else "not closed")
