# msolve rational parametrisation mod p; remove degenerate solutions; count the rest.
import pickle,sys,subprocess,time,ast,itertools,re
import flint
from ms import system
BASE="/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/"
def parse_poly(s,vars_,p):
    # expanded string "c*x^a*y^b+..." -> list of (c, {var:exp})
    s=s.replace("-","+-")
    out=[]
    for t in s.split("+"):
        t=t.strip()
        if not t: continue
        fac=t.split("*"); c=1; mon={}
        for f in fac:
            if re.fullmatch(r"-?\d+",f): c*=int(f)
            elif f.startswith("-"):
                c*=-1; f=f[1:]
                nm,_,ex=f.partition("^"); mon[nm]=mon.get(nm,0)+(int(ex) if ex else 1)
            else:
                nm,_,ex=f.partition("^"); mon[nm]=mon.get(nm,0)+(int(ex) if ex else 1)
        out.append((c%p,mon))
    return out
def evalpoly(terms,X,w,p,cache):
    val=flint.nmod_poly([0],p)
    for c,mon in terms:
        t=flint.nmod_poly([c],p)
        for nm,ex in mon.items():
            key=(nm,ex)
            if key not in cache: cache[key]=pow(X[nm],ex,w) if False else (X[nm]**ex)%w
            t=(t*cache[key])%w
        val=(val+t)%w
    return val
def count(cols,h,p,tag,slices=0,seed=1,npts=5,timeout=3000,verify=True):
    v,e=system(cols,h,slices,seed,npts)
    fn=BASE+"q_"+tag
    open(fn+".in","w").write(",".join(v)+"\n"+str(p)+"\n"+",\n".join(e)+"\n")
    t0=time.time()
    try:
        r=subprocess.run(["prlimit","--as=3500000000","msolve","-f",fn+".in","-o",fn+".out"],capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        return dict(status="TIMEOUT",time=time.time()-t0)
    s=open(fn+".out").read().strip().rstrip(':')
    d=ast.literal_eval(s)
    if d==[-1] or d[0]==-1: return dict(status="EMPTY(before sat)",nsol=0,time=time.time()-t0)
    if d[0]!=0 or (len(d)>1 and isinstance(d[1],list) and len(d[1])>1 and d[1][2]==-1):
        return dict(status="POSDIM",raw=str(d)[:100],time=time.time()-t0)
    P=d[1]; vars_=P[3]; lin=P[4]
    if sorted(lin)==[0]*(len(lin)-1)+[1]: j=lin.index(1)
    elif vars_[-1]=='A': j=len(vars_)-1
    else: return dict(status="NONUNIT_LINFORM",time=time.time()-t0)
    Q=P[5][1]; w=flint.nmod_poly(Q[0][1],p)
    nums=[flint.nmod_poly(x[0][1],p) for x in Q[2]]
    others=[x for i,x in enumerate(vars_) if i!=j]
    X={nm:(-num)%w for nm,num in zip(others,nums)}
    X[vars_[j]]=flint.nmod_poly([0,1],p)
    if vars_[j]=='A': e=e  # A not in equations
    deg0=w.degree()
    cache={}
    if verify:
        for eq in e:
            val=evalpoly(parse_poly(eq,vars_,p),X,w,p,cache)
            if val.degree()>=0: return dict(status="VERIFY_FAIL",time=time.time()-t0)
    n=len(cols)
    for k,l in itertools.combinations(range(n),2):
        for sx,sy in [(1,1),(1,-1),(1,1j)]:
            pass
    # degeneracy forms in x,y: r=x+y, s=x-y
    forms=[]
    for k,l in itertools.combinations(range(n),2):
        dx=(X[f"x{k}"]-X[f"x{l}"])%w
        for a,b in [(1,1),(-1,-1),(1,-1),(-1,1)]:
            forms.append((dx+a*X[f"y{k}"]-b*X[f"y{l}"])%w)
    for k in range(n): forms.append(X[f"y{k}"])
    ww=w
    for f in forms:
        g=ww.gcd(f)
        if g.degree()>0: ww=ww//g
        if ww.degree()==0: break
    return dict(status="OK",nsol_total=deg0,nsol_nondeg=ww.degree(),time=time.time()-t0,w=ww,X=X,wfull=w)
if __name__=="__main__":
    from cuttypes import CUTS
    from gb import cuts_to_columns
    D,h,p,sl=map(int,sys.argv[1:5])
    js=sys.argv[5] if len(sys.argv)>5 else "all"
    T=pickle.load(open(BASE+f'types_{D}_{h}.pkl','rb'))
    rng=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
    for j in rng:
        cols=cuts_to_columns(T[j],CUTS)
        r=count(cols,h,p,f"{D}_{h}_{j}_{p}_{sl}",slices=sl)
        r.pop('w',None); r.pop('X',None); r.pop('wfull',None)
        print(D,h,j,["".join(map(str,a)) for a,e in cols],r,flush=True)
