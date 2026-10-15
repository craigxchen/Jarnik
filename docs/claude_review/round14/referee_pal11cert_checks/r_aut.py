# Referee: signed-permutation automorphisms of A (pal11, L=2), by backtracking with pair statistics, then exact check.
import pickle, numpy as np
d=pickle.load(open('r_box.pkl','rb')); A=d['A'].astype(int); n=11
N=np.zeros((n,n,3,3),int)
for x in range(n):
    for y in range(n):
        for s in (-1,0,1):
            for t in (-1,0,1):
                N[x,y,s+1,t+1]=np.count_nonzero((A[:,x]==s)&(A[:,y]==t))
Aset=set(map(tuple,A)); sols=[]
def img(x,eps,s): return eps*s
def ok_pair(x,px,ex,y,py,ey):
    for s in (-1,0,1):
        for t in (-1,0,1):
            if N[x,y,s+1,t+1]!=N[px,py,ex*s+1,ey*t+1]: return False
    return True
def rec(x,perm,sgn,used):
    if x==n:
        P=np.array(perm); S=np.array(sgn)
        B=np.zeros_like(A); B[:,P]=A*S
        if set(map(tuple,B))==Aset: sols.append((tuple(perm),tuple(sgn)))
        return
    for px in range(n):
        if px in used: continue
        for ex in (1,-1):
            if not ok_pair(x,px,ex,x,px,ex): continue
            if all(ok_pair(x,px,ex,y,perm[y],sgn[y]) for y in range(x)):
                perm.append(px); sgn.append(ex); used.add(px)
                rec(x+1,perm,sgn,used)
                perm.pop(); sgn.pop(); used.discard(px)
rec(0,[],[],set())
print('automorphisms (signed perms of k\' coords preserving A):',len(sols))
for s in sols: print(s)
