# cut-group (coefficient) formulation: for each cut type tau with multiplicity n, monic G_tau, H_tau of degree n.
import subprocess,time,random,sys,pickle
from cuttypes import CUTS
BASE="/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/"
def build(groups,h,char,npts=5,slices=0,seed=1,satall=True):
    # groups: list of (set_of_points_with_G, n)
    V=[];Gs=[];Hs=[]
    for k,(S,n) in enumerate(groups):
        g=[f"g{k}_{j}" for j in range(1,n+1)]; hh=[f"h{k}_{j}" for j in range(1,n+1)]
        V+=g+hh
        Gs.append("(t^%d"%n+"".join(f"+{c}*t^{n-j}" for j,c in enumerate(g,1))+")")
        Hs.append("(t^%d"%n+"".join(f"+{c}*t^{n-j}" for j,c in enumerate(hh,1))+")")
    k0=min(range(len(groups)),key=lambda k:groups[k][1])
    lines=[]
    lines.append("ideal I;")
    for i in range(1,npts):
        inc=[k for k,(S,n) in enumerate(groups) if i in S]
        d=sum(groups[k][1] for k in inc)
        lines.append(f"poly P={'*'.join(Gs[k] for k in inc)}-{'*'.join(Hs[k] for k in inc)};")
        lines.append(f"for(int m=1;m<{h};m++){{ I=I+jet(subst(diff(P,t),t,0),0); }}" if False else "")
        lines.append(f"matrix C=coeffs(P,t); for(int m=1;m<{h};m++){{ if ({d}-m+1>=1) {{ I=I+C[{d}-m+1,1]; }} }} kill C; kill P; kill m;")
    lines.append("I=I+("+"+".join(f"g{k}_1+h{k}_1" for k in range(len(groups)))+");")
    n0=groups[k0][1]
    lines.append(f"I=I+(resultant({Gs[k0]},{Hs[k0]},t)-1);")
    rnd=random.Random(seed)
    for _ in range(slices):
        lines.append("I=I+("+"+".join(f"{rnd.randint(1,9999)}*{v}" for v in V)+f"+{rnd.randint(1,9999)});")
    sat=[f"resultant({Gs[k]},{Hs[k]},t)" for k in range(len(groups)) if k!=k0]
    if satall:
        for a in range(len(groups)):
            for b in range(a+1,len(groups)):
                sat+= [f"resultant({Gs[a]},{Gs[b]},t)",f"resultant({Hs[a]},{Hs[b]},t)",f"resultant({Gs[a]},{Hs[b]},t)",f"resultant({Hs[a]},{Gs[b]},t)"]
        for k,(S,n) in enumerate(groups):
            if n>1: sat+= [f"resultant({Gs[k]},diff({Gs[k]},t),t)",f"resultant({Hs[k]},diff({Hs[k]},t),t)"]
    scr=f'''LIB "elim.lib";
option(redSB);
ring rr={char},(t,{",".join(V)}),dp;
{chr(10).join(lines)}
ring r2={char},({",".join(V)}),dp;
ideal I=imap(rr,I);
I=std(I);
string s0="before: dim="+string(dim(I))+" vdim="+string(vdim(I));
setring rr;
list F={",".join(sat) if sat else "1"};
setring r2;
list F=imap(rr,F);
int j;
for(j=1;j<=size(F);j++){{ if (deg(I[1])==0) {{break;}} I=sat(I,F[j]); I=std(I); }}
print(s0+"  after: dim="+string(dim(I))+" vdim="+string(vdim(I))+" one="+string(deg(I[1])==0));
quit;
'''
    return scr
def run(groups,h,char,tag,timeout=3600,**kw):
    scr=build(groups,h,char,**kw)
    fn=BASE+"cf_"+tag+".sing"; open(fn,"w").write(scr)
    t=time.time()
    try: out=subprocess.run(["prlimit","--as=3500000000","Singular","-q",fn],capture_output=True,text=True,timeout=timeout).stdout
    except subprocess.TimeoutExpired: out="TIMEOUT"
    return out.strip().replace("\n"," | ")[:400],time.time()-t
def groups_of(ty):
    return [(set(CUTS[k]),c) for k,c in enumerate(ty) if c]
if __name__=="__main__":
    D,h,char=map(int,sys.argv[1:4]); js=sys.argv[4] if len(sys.argv)>4 else "all"
    to=int(sys.argv[5]) if len(sys.argv)>5 else 3600
    T=pickle.load(open(BASE+f"types_{D}_{h}.pkl","rb"))
    rng=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
    for j in rng:
        g=groups_of(T[j])
        out,tm=run(g,h,char,f"{D}_{h}_{j}_{char}",timeout=to)
        print(D,h,j,[("".join(str(int(i in S)) for i in range(5)),n) for S,n in g],f"{tm:.1f}s",out,flush=True)
