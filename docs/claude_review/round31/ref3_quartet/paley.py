"""Own construction of the ordinary rows of the round-8 skew double core (lw.md 5.1), checks of the
cut counts used in III.2, the five-mark structure used in III.7, the lcm statement (III.3) and the 32-pattern count."""
import random
from itertools import combinations, product
from collections import Counter

def chi(x,q):
    x%=q
    if x==0: return 0
    return 1 if pow(x,(q-1)//2,q)==1 else -1

def columns(q):
    # rows = F_q (border row deleted); s1_a(x)=[x=a]+chi(a-x); s2_a = s1_a negated at row a
    cols=[]
    for a in range(q):
        s1=[(1 if x==a else chi(a-x,q)) for x in range(q)]
        s2=[(-1 if x==a else chi(a-x,q)) for x in range(q)]
        cols+=[s1,s2]
    return cols

def unoriented(c): return min(tuple(c),tuple(-v for v in c))

def qtype(col,quad):
    b=[col[i] for i in quad]
    if sum(b)!=0: return None
    return [j for j in (1,2,3) if b[j]==b[0]][0]   # partner of first mark

def check_counts(q):
    cols=columns(q); r=2*q
    assert len({unoriented(c) for c in cols})==r
    sep=[sum(c[i]!=c[j] for c in cols) for i,j in combinations(range(q),2)]
    assert set(sep)=={q+1}
    mins=[]
    for quad in combinations(range(q),4):
        cnt=Counter(qtype(c,quad) for c in cols)
        mins.append(min(cnt[1],cnt[2],cnt[3]))
    return min(mins), r/16, r/8

def five_mark(q=7):
    cols=columns(q); hist=Counter(); bad=0
    for quad in combinations(range(q),4):
        hist[Counter(qtype(c,quad) for c in cols)[1]]+=1
    for five in combinations(range(q),5):
        parts=[]
        for quad in combinations(five,4):
            cl={1:[],2:[],3:[]}
            for j,c in enumerate(cols):
                t=qtype(c,quad)
                if t: cl[t].append(j)
            if all(len(v)==3 for v in cl.values()):
                parts.append(frozenset(frozenset(v) for v in cl.values()))
        assert len(parts)==2
        bad+= parts[0]==parts[1]
    return dict(hist), bad

def lcm_III3(trials=20000):
    # v_p(B(ij|kx)) = h3-h2 if the low/high pairing of the 4 exponents is {i,j}|{k,x}; else 0
    def vB(e,pair,other):
        four=list(pair)+list(other); vals=sorted(e[t] for t in four)
        # pairing low-low/high-high; ties: if h2==h3 exponent 0 anyway
        if vals[1]==vals[2]: return 0
        lowset={t for t in four if e[t]<=vals[1]}
        return vals[2]-vals[1] if lowset in (set(pair),set(other)) else 0
    fails=0
    for _ in range(trials):
        m=random.randint(5,8); E=random.randint(1,6)
        e=[random.randint(0,E) for _ in range(m)]
        i,j,x,y=0,1,2,3
        lhs=max(min(vB(e,(i,j),(k,x)),vB(e,(i,j),(k,y))) for k in range(4,m))
        if lhs>vB(e,(i,j),(x,y)): fails+=1
    return fails

def bits32():
    def kind(v):
        if sum(v)!=2: return '0'
        return 'ABC'[[j for j in (1,2,3) if v[j]==v[0]][0]-1]
    return Counter((kind(v[:4]),kind(v[:3]+v[4:])) for v in product((0,1),repeat=5))

if __name__=='__main__':
    for q in (7,11,19,23,43,47,67,71,79,83,103,107,127,131):
        print('q',q,'min class size over quartets, r/16, r/8:',check_counts(q))
    print('five-mark (q=7): quartet histogram by per-class size, # five-sets with equal 3-class partitions:',five_mark())
    print('III.3 lcm counterexamples in random prime-power trials:',lcm_III3())
    print('32-pattern count:',sorted(bits32().items()))
