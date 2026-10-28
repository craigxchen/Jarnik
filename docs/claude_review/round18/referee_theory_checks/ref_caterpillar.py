# Independent check of theory.md Lemma 1.1 (referee code, written from scratch).
# Enumerate all compatible systems of nontrivial splits of {0..n-1} (= all trees),
# test (a) two-chain condition for anchor 0, (b) for all anchors, (c) one-chain choice of sides,
# and (d) internal-degree<=2 caterpillar criterion computed from an explicit tree build.
import itertools, sys
def run(n):
    full=(1<<n)-1
    # splits labelled by side not containing 0, size 2..n-2
    S=[m for m in range(1,full+1) if not (m&1) and 2<=bin(m).count('1')<=n-2]
    def comp(a,b):
        return (a&b)==0 or (a&b)==a or (a&b)==b or ((full^a)&(full^b))==0
    idx={m:i for i,m in enumerate(S)}
    systems=[]
    def bt(start,cur):
        systems.append(list(cur))
        for i in range(start,len(S)):
            m=S[i]
            if all(comp(m,c) for c in cur):
                cur.append(m); bt(i+1,cur); cur.pop()
    bt(0,[])
    def chain(lst):
        l=sorted(lst,key=lambda m:bin(m).count('1'))
        return all((l[i]&l[i+1])==l[i] for i in range(len(l)-1))
    def cond_a(sys_,a):
        bit=1<<a
        sides=[ (m if not (m&bit) else full^m) for m in sys_]
        r=len(sides)
        for mask in range(1<<r):
            g1=[sides[i] for i in range(r) if mask>>i&1]; g2=[sides[i] for i in range(r) if not mask>>i&1]
            u1=0
            for x in g1: u1|=x
            u2=0
            for x in g2: u2|=x
            if chain(g1) and chain(g2) and (u1&u2)==0: return True
        return False
    def cond_c(sys_):
        r=len(sys_)
        for mask in range(1<<r):
            sides=[ (sys_[i] if mask>>i&1 else full^sys_[i]) for i in range(r)]
            if chain(sides): return True
        return False
    def cond_d(sys_):
        # build tree: vertices = maximal consistent ... use: internal-edge graph degree.
        # Two internal edges e,f are adjacent (share a vertex) iff no third split separates them,
        # i.e. no g in sys with g strictly 'between' e and f. For compatible splits viewed as
        # sides not containing 0 (rooted at 0): use rooted-tree: parent of side X = minimal side
        # strictly containing X (or root). Internal vertex v_X = lower endpoint of edge X.
        # Build rooted tree on sides plus root R (vertex where leaf 0 attaches).
        sides=list(sys_)
        par={}
        for X in sides:
            sup=[Y for Y in sides if Y!=X and (X&Y)==X]
            par[X]=min(sup,key=lambda m:bin(m).count('1')) if sup else 'R'
        # vertices: 'R' and lower endpoint of each edge X (call X). internal edge X joins X to par[X].
        deg={}
        for X in sides:
            deg[X]=deg.get(X,0)+1
            deg[par[X]]=deg.get(par[X],0)+1
        return all(d<=2 for d in deg.values())
    res={'trees':len(systems),'cat':0,'mism':0}
    for s in systems:
        a0=cond_a(s,0); ball=all(cond_a(s,a) for a in range(n)); c=cond_c(s); d=cond_d(s)
        if not (a0==ball==c==d): res['mism']+=1
        res['cat']+=c
    return res
for n in range(5,9):
    print(n,run(n)); sys.stdout.flush()
