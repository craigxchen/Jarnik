# reg(S) for finite S in Z^n: least d with Poly_{<=d}|_S = F^S.
# rank mod p <= rank over Q: full rank mod p at degree d certifies reg(S) <= d over Q.
# Rank deficiency mod p at d-1 does NOT certify reg >= d over Q; lower bounds must come from explicit f.
import itertools
from flint import nmod_mat
P=(1<<61)-1
def monos(n,d):
    out=[]
    def rec(i,rem,cur):
        if i==n-1:
            for a in range(rem+1): out.append(tuple(cur+[a]))
            return
        for a in range(rem+1): rec(i+1,rem-a,cur+[a])
    rec(0,d,[]); return out
def rank_at(pts,d,p=P):
    n=len(pts[0]); ms=monos(n,d)
    rows=[]
    for x in pts:
        pw=[[pow(xi%p,a,p) for a in range(d+1)] for xi in x]
        row=[]
        for a in ms:
            v=1
            for i in range(n): v=v*pw[i][a[i]]%p
            row.append(v)
        rows.append(row)
    return nmod_mat(rows,p).rank()
def reg(pts,dstart=0,p=P):
    pts=[tuple(x) for x in set(map(tuple,pts))]
    N=len(pts)
    if N==1: return 0
    d=dstart
    while rank_at(pts,d,p)<N: d+=1
    # ensure minimal: step down
    while d>0 and rank_at(pts,d-1,p)==N: d-=1
    return d
def width(pts,v):
    vals=[sum(a*b for a,b in zip(x,v)) for x in pts]; return max(vals)-min(vals)
