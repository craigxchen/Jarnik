# Judge check J3 (unrefereed): for all 105 boxes, rebuild A with own code from the compute lane's profile
# generator (pats only), compare |A| with boxes.csv, and verify every zeta witness exactly:
# all 2^t vertices in A, every generator nonzero, pair witnesses in L_pair. Witness of order t => zeta >= t.
import sys, json, csv, itertools
D='/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28'
sys.path.insert(0, D+'/compute_checks')
import run_box
rows={r['spec']:r for r in csv.DictReader(open(D+'/boxes.csv'))}
Z={}
for line in open(D+'/compute_checks/zeta.jsonl'):
    r=json.loads(line); Z[r['spec']]=r
def mybox(n,pats,L):
    d=n-1; A=set()
    for k in itertools.product((-1,0,1),repeat=d):
        kf=list(k)+[-sum(k)]
        if all(abs(sum(kf[x] for x in S))<=L for S in pats): A.add(k)
    return A
def in_pair(v):
    nz=[x for x in v if x]; return (len(nz)==1 and abs(nz[0])==1) or (len(nz)==2 and sorted(nz)==[-1,1])
def check(A,k0,gens):
    vs=[]
    for g,m in gens: vs+=[tuple(g)]*m
    pts={tuple(k0)}
    for v in vs: pts=pts|{tuple(a+b for a,b in zip(p,v)) for p in pts}
    # pts = set of all subset sums (as a set); membership of every subset sum <=> pts subset A
    return len(vs), pts<=A and all(any(v) for v in vs), all(in_pair(v) for v in vs)
nbad=0; nok=0; lt=0
for spec,row in rows.items():
    n,pats,L=run_box.parse(spec)
    A=mybox(n,[list(S) for S in pats],L)
    if len(A)!=int(row['size']): print('SIZE MISMATCH',spec,len(A),row['size']); nbad+=1; continue
    z=Z.get(spec)
    if z is None: print('no zeta record',spec); nbad+=1; continue
    for cls,col in (('pair','zeta_pair'),('all','zeta_all')):
        w=z[cls]
        t,ok,pr=check(A,w['k0'],w['gens'])
        good= ok and t==int(row[col]) and (pr or cls=='all')
        if not good: print('BAD',spec,cls,t,ok,pr,row[col]); nbad+=1
        else: nok+=1
    if int(row['reach_pair_P1'])<int(row['zeta_pair']): lt+=1
print('verified witnesses',nok,'bad',nbad,'; boxes with reach_pair < zeta_pair certified by witness:',lt)
