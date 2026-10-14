# EVIDENCE (mod-p J of images): which subfamilies of {y1-y2, y3..y7} (with all GL(3,2) images) close the Fano gap?
import itertools, sys, time, pickle
from rlib import *
p=2147483647
Y,h,dJ6,dJ7=pickle.load(open('fano_Y.pkl','rb'))
gens=pickle.load(open('fano_gens.pkl','rb'))
cap=[5]*7; D=7
cols=monos_eq(7,D,cap); cidx={c:i for i,c in enumerate(cols)}
avec=[a for a in itertools.product((0,1),repeat=3) if any(a)]
print("coordinate labels",avec)
def mat_apply(Mx,a): return tuple(sum(Mx[r][c]*a[c] for c in range(3))%2 for r in range(3))
perms=[]
for Mx in itertools.product(itertools.product((0,1),repeat=3),repeat=3):
    imgs=[mat_apply(Mx,a) for a in avec]
    if len(set(imgs))==7 and all(any(x) for x in imgs): perms.append([avec.index(x) for x in imgs])
def J_space_p(X,e):
    ms=monos_le(len(X[0]),e)
    M=flint.nmod_mat([[evalmono(q,m)%p for m in ms] for q in X],p)
    N,nul=M.nullspace()
    top=[i for i,m in enumerate(ms) if sum(m)==e]
    if nul==0: return []
    R=flint.nmod_mat([[int(N[i,c]) for i in top] for c in range(nul)],p)
    Rr,rk=R.rref()
    return [{ms[top[c]]:int(Rr[r,c]) for c in range(len(top)) if int(Rr[r,c])} for r in range(rk)]
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
B0=basis_of(rows_to_list(ideal_degree_rows(gens,7,D,cap)))
e=lambda i:[1 if q==i else 0 for q in range(7)]
d0=[1,-1,0,0,0,0,0]
others=[e(i) for i in range(2,7)]
for size in range(1,6):
    for S in itertools.combinations(range(5),size):
        F=[d0]+[others[i] for i in S]
        img=sorted(set(tuple(sum(a*b for a,b in zip(v,y)) for v in F) for y in Y))
        newg={}
        for ee in range(1,D+1):
            sp=J_space_p(img,ee)
            if sp: newg[ee]=[pullback(g,F,7) for g in sp]
        base=ideal_degree_rows(newg,7,D,cap)
        if not base: 
            print("b=%d"%(size+1),"S",S,"no gens <=7"); continue
        B=list(B0)
        for pm in perms:
            rows=[]
            for r in base:
                nr={}
                for ex,v in r.items():
                    ne=[0]*7
                    for i,a in enumerate(ex): ne[pm[i]]=a
                    nr[tuple(ne)]=v
                rows.append(nr)
            B=basis_of(B+rows_to_list(rows))
            if len(B)>=dJ7: break
        print("b=%d"%(size+1),"extra coords",[avec[i+2] for i in S],"|img|",len(img),"rank",len(B),"/",dJ7); sys.stdout.flush()
