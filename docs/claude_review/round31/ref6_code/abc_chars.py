"""(e) abc support test  w0(v) >= W/4 - (delta/2) w1(v) - log(|c|_1 C) - O(1)  on light characters, equal weights.
FF(q,b): construct.md fully flipped Paley profile; DC(q): skew double core (record Thm 5.1: w0 >= M/2 columns)."""
import sys
sys.path.insert(0,'.')
from profiles import flipped, double_core, paley_H
def chars(S, cs):
    r=len(S[0]); out=[]
    for name,c in cs:
        assert sum(c)==0
        v=[sum(c[x]*S[x][j] for x in range(len(S))) for j in range(r)]
        assert all(t%2==0 for t in v); v=[t//2 for t in v]
        w0=sum(1 for t in v if t)/r; w1=sum(abs(t) for t in v)/r
        out.append((name,sum(abs(t) for t in c),w0,w1))
    return out
for q in (43,):
    for b in (3,4,6,10):
        S=flipped(q,b); Hn=paley_H(q); M=q+1
        cs=[("label col H(.,a)",[Hn[x][a] for x in range(M)]) for a in (2,5)]
        for name,l1,w0,w1 in chars(S,cs):
            print(f"FF(q={q},b={b}) {name}: |c|_1={l1}  w0/W={w0:.3f}  w1/W={w1:.3f}   (abc needs w0/W >= 1/4 - (delta/2) w1/W - o(1))")
    S=double_core(q,True); M=q+1
    import itertools
    K=[[None]]
    cs=[("pair e_x-e_y",[1 if x==1 else (-1 if x==2 else 0) for x in range(M)])]
    for name,l1,w0,w1 in chars(S,cs):
        print(f"DC(q={q}) {name}: |c|_1={l1} w0/W={w0:.3f} w1/W={w1:.3f}")
