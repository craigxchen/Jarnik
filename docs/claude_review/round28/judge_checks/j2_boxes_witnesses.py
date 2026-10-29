# Judge check J2 (unrefereed): rebuild Sylvester restriction boxes from scratch, verify zeta witnesses
# exactly (all 2^t vertices in A, generators primitive, nonzero), and the face claims of verify.md §3.
import itertools, json, sys
D='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28'
def dot(u,x): return bin(u&x).count('1')&1
def syl_cuts(k): return [frozenset(x for x in range(2**k) if dot(u,x)) for u in range(1,2**k)]
def restrict(cuts,J):
    Js=set(J); out=set()
    for S in cuts:
        T=frozenset(S&Js)
        if not T or T==Js: continue
        C=frozenset(Js-T)
        out.add(min(T,C,key=lambda s:(len(s),sorted(s))))
    return sorted(out,key=lambda s:(len(s),sorted(s)))
def box(J,cuts,L):
    d=len(J)-1; base=J[-1]
    pos={x:i for i,x in enumerate(J)}
    A=set()
    for k in itertools.product((-1,0,1),repeat=d):
        kf=list(k)+[-sum(k)]
        if abs(kf[-1])>1: continue   # base is a cube point too? (checked below against |A|)
        if all(abs(sum(kf[pos[x]] for x in S))<=L for S in cuts): A.add(k)
    return A
def box_nobasecube(J,cuts,L):
    d=len(J)-1; pos={x:i for i,x in enumerate(J)}
    A=set()
    for k in itertools.product((-1,0,1),repeat=d):
        kf=list(k)+[-sum(k)]
        if all(abs(sum(kf[pos[x]] for x in S))<=L for S in cuts): A.add(k)
    return A
def verify_witness(A,k0,gens):
    vs=[]
    for g,m in gens: vs+= [tuple(g)]*m
    t=len(vs)
    for T in itertools.product((0,1),repeat=t):
        p=tuple(k0[i]+sum(vs[j][i] for j in range(t) if T[j]) for i in range(len(k0)))
        if p not in A: return t,False
    return t,all(any(v) for v in vs)
def in_pair(v):
    nz=[x for x in v if x]
    return (len(nz)==1 and abs(nz[0])==1) or (len(nz)==2 and sorted(nz)==[-1,1])
res={}
c8=syl_cuts(3); c16=syl_cuts(4)
M8p_c=restrict(c8,list(range(8)))
A_M8p=box_nobasecube(list(range(8)),M8p_c,3); A_M8p_cb=box(list(range(8)),M8p_c,3)
print('M8p |A| (base unconstrained / base in cube):',len(A_M8p),len(A_M8p_cb))
J9a=list(range(9)); J9ap=[1,2,3,4,5,6,7,8,0]
cA=restrict(c16,J9a); cAp=restrict(c16,J9ap)
for nm,J,c in (('J9a',J9a,cA),('J9a_prime',J9ap,cAp)):
    A1=box_nobasecube(J,c,3); A2=box(J,c,3)
    print(nm,'cuts',len(c),'|A| base-free / base-in-cube:',len(A1),len(A2))
    res[nm]=(A1,A2)

# ---- faces (exact set comparison) ----
AJ=res['J9a'][0]; AJp=res['J9a_prime'][0]
# J9a: coords points 0..7, base 8. Face {k_base = 0} = {sum k' = 0}; drop coord of point 7 (M8' base).
faceA={k[:7] for k in AJ if sum(k)==0}
print('J9a face {k_8base=0} -> M8p subset?',faceA<=A_M8p,'size',len(faceA),'== M8p cap |k7|<=1:',faceA==A_M8p_cb)
# J9a': coords points 1..8, base 0. Face {k(point 8)=0}: coords points 1..7, base 0 = Sylvester H_8 on {0..7} with base 0.
faceAp={k[:7] for k in AJp if k[7]==0}
M8p_base0=box_nobasecube([1,2,3,4,5,6,7,0],restrict(c8,[1,2,3,4,5,6,7,0]),3)
print("J9a' face {k_8=0} == M8' (base 0 labelling):",faceAp==M8p_base0,len(faceAp))
# M8' with base 0 vs base 7: translation x->x^7 maps point list (0..6, base 7) to (7,6,...,1, base 0)
tr=[x^7 for x in range(8)]  # point x -> x^7
cmp={tuple(k[tr[i]-1] for i in range(7)) for k in M8p_base0}  # coordinate of point i in base-7 frame = coord of point i^7 in base-0 frame (index (i^7)-1)
print("M8' base-0 box is the translate of the base-7 box:",cmp==A_M8p)
# ---- witnesses ----
def load(path,spec):
    out=None
    for line in open(path):
        r=json.loads(line)
        if r.get('spec')==spec: out=r
    return out
w=load(D+'/compute_checks/zeta.jsonl','rsyl:4:0,1,2,3,4,5,6,7,8:3')
for cls in ('pair','all'):
    t,ok=verify_witness(AJ,w[cls]['k0'],w[cls]['gens'])
    pairs=all(in_pair(g) for g,m in w[cls]['gens'])
    print('J9a',cls,'witness order',t,'all vertices in A & nonzero gens:',ok,'all gens in L_pair:',pairs)
wm=load(D+'/compute_checks/zeta.jsonl','syl:3:3')
for cls in ('pair','all'):
    t,ok=verify_witness(A_M8p,wm[cls]['k0'],wm[cls]['gens'])
    print("M8'",cls,'witness order',t,'ok',ok,'L_pair:',all(in_pair(g) for g,m in wm[cls]['gens']))
wp=load(D+'/verify_checks/zeta.jsonl','rsyl:4:1,2,3,4,5,6,7,8,0:3')
print("J9a' verify keys:",list(wp.keys())[:12])
k0=tuple(wp['wit_pair'][0]); gens=[(g,1) for g in wp['wit_pair'][1]]
t,ok=verify_witness(AJp,k0,gens)
print("J9a' pair witness order",t,'ok',ok,'L_pair:',all(in_pair(g) for g,m in gens))
# subset criterion consequence: the J9a' face is the M8' box (as a subset of Z^8 lying in {k_8=0})
