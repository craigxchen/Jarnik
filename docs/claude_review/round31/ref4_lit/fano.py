# Incidence of Astra's trace construction: marks a in F_q^*, primes x in F_q^*, x marked by a iff Tr(ax)=0.
# Check: it is the point-hyperplane design of PG(d-1,2) (Fano plane for d=3); and evaluate (IV.2) at delta=E=0.
from itertools import combinations
def gf(d):
    polys={3:0b1011,4:0b10011,5:0b100101,6:0b1000011,7:0b10001001}
    P=polys[d]; q=1<<d
    def mul(x,y):
        r=0
        while y:
            if y&1: r^=x
            y>>=1; x<<=1
            if x&q: x^=P
        return r
    def tr(x):
        s=0; y=x
        for _ in range(d): s^=y; y=mul(y,y)
        return s
    return q,mul,tr
for d in [3,4,5,6,7]:
    q,mul,tr=gf(d)
    assert all(tr(x) in (0,1) for x in range(q))
    blocks={a:frozenset(x for x in range(1,q) if tr(mul(a,x))==0) for a in range(1,q)}
    sizes={len(B) for B in blocks.values()}
    inter={len(blocks[a]&blocks[b]) for a,b in combinations(range(1,q),2)}
    # F2-linearity: block u {0} closed under xor
    lin=all(((x^y) in blocks[a]) for a in blocks for x in blocks[a] for y in blocks[a] if x!=y)
    m=q-1; L=q-1  # layers: one per linear prime, s_h = #marks at x = q/2-1
    s_h=q//2-1
    lhs=0.75*(m-24)*L+L*(s_h-m/2)**2
    print(f"d={d} q={q}: block sizes {sizes}, pair intersections {inter}, F2-subspace minus 0: {lin}, (IV.2) LHS at delta=E=0: {lhs}  (violated iff >0)")
