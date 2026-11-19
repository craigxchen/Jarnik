# Build and run Singular saturation test for a 5-point layer type.
# columns: list of (vector a in [0,e]^5, e); points z_i = prod (t-r_k)^{a_k^i} (t-s_k)^{e_k-a_k^i}
import subprocess, sys, itertools, os, tempfile, time
def build(columns,h,char=0,extra_sat=True,gauge=True,order='dp',slices=0,seed=1,npts=5):
    n=len(columns)
    R=[f"r{k}" for k in range(n)]; S=[f"s{k}" for k in range(n)]
    eqs=[]
    for i in range(1,npts):
        for m in range(1,h):
            terms=[]
            for k,(a,e) in enumerate(columns):
                c=a[i]-a[0]
                if c: terms.append(f"({c})*({R[k]}^{m}-{S[k]}^{m})")
            if terms: eqs.append("+".join(terms))
    if gauge:
        eqs.append("+".join(f"{e}*({R[k]}+{S[k]})" for k,(a,e) in enumerate(columns)))
        eqs.append(f"{R[0]}-{S[0]}-1")
    import random
    rnd=random.Random(seed)
    for _ in range(slices):
        eqs.append("+".join(f"{rnd.randint(1,10**4)}*{v}" for v in R+S)+f"+{rnd.randint(1,10**4)}")
    forms=[f"{R[k]}-{S[k]}" for k in range(n)]
    for k,l in itertools.combinations(range(n),2):
        forms+= [f"{R[k]}-{R[l]}",f"{S[k]}-{S[l]}",f"{R[k]}-{S[l]}",f"{S[k]}-{R[l]}"]
    eqs.append("u*"+"*".join(f"({R[k]}-{S[k]})" for k in range(n))+"-1")
    vars_=",".join(R+S+["u"])
    scr=f'''LIB "elim.lib";
option(redSB);
ring rr={char},({vars_}),{order};
ideal I={",".join(eqs)};
I=std(I);
print("before: dim="+string(dim(I)));
list F={",".join(forms)};
int j;
for(j=1;j<=size(F);j++){{ I=sat(I,F[j]); I=std(I); if (deg(I[1])==0) {{break;}} }}
print("after: dim="+string(dim(I))+" vdim="+string(vdim(I))+" one="+string(deg(I[1])==0));
if (deg(I[1])!=0 and dim(I)==0) {{ print(I); }}
quit;
'''
    return scr
def run(columns,h,char=0,timeout=3600,tag="x",**kw):
    scr=build(columns,h,char,**kw)
    fn=f"/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/sing_{tag}.sing"
    open(fn,'w').write(scr)
    t=time.time()
    try:
        out=subprocess.run(["prlimit","--as=3500000000","Singular","-q",fn],capture_output=True,text=True,timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        out="TIMEOUT"
    return out, time.time()-t
def cuts_to_columns(counts,CUTS):
    cols=[]
    for k,c in enumerate(counts):
        for _ in range(c):
            cols.append((tuple(1 if i in CUTS[k] else 0 for i in range(5)),1))
    return cols
