# msolve GB (with u*prod(y)=1 and optional slices) then Singular saturation by remaining degeneracy forms
import pickle,sys,subprocess,time,re,os
from ms import system, forms
BASE="/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/"
def analyse(cols,h,char,tag,slices=0,seed=1,npts=5,timeout=3000,extra_forms=None):
    v,e=system(cols,h,slices,seed,npts)
    fn=BASE+"m_"+tag
    open(fn+".in","w").write(",".join(v)+"\n"+str(char)+"\n"+",\n".join(e)+"\n")
    t=time.time()
    try:
        r=subprocess.run(["prlimit","--as=3500000000","msolve","-f",fn+".in","-o",fn+".out","-g","2","-t","1"],capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        return "MSOLVE_TIMEOUT",time.time()-t
    txt=open(fn+".out").read()
    body=txt[txt.index("["):txt.rindex("]")+1].strip()[1:-1]
    gens=[g.strip() for g in body.split(",\n")]
    if gens==["1"] or gens==["1*1"]:
        return "GB=1 before saturation",time.time()-t
    F=forms(len(cols)) + (extra_forms or [])
    scr=f'''LIB "elim.lib";
option(redSB);
ring rr={char},({",".join(v)}),dp;
ideal I={",".join(gens)};
I=std(I);
string s0="before: dim="+string(dim(I))+" vdim="+string(vdim(I));
list F={",".join(F)};
int j;
for(j=1;j<=size(F);j++){{ if (deg(I[1])==0) {{break;}} I=sat(I,F[j]); I=std(I); }}
print(s0+"  after: dim="+string(dim(I))+" vdim="+string(vdim(I))+" one="+string(deg(I[1])==0));
quit;
'''
    open(fn+".sing","w").write(scr)
    try:
        out=subprocess.run(["prlimit","--as=3500000000","Singular","-q",fn+".sing"],capture_output=True,text=True,timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        out="SING_TIMEOUT"
    return out.strip().replace("\n"," | "),time.time()-t
if __name__=="__main__":
    from cuttypes import CUTS
    from gb import cuts_to_columns
    D,h,char,sl=map(int,sys.argv[1:5])
    js=sys.argv[5] if len(sys.argv)>5 else "all"
    T=pickle.load(open(BASE+f'types_{D}_{h}.pkl','rb'))
    rng=range(len(T)) if js=="all" else [int(x) for x in js.split(",")]
    for j in rng:
        cols=cuts_to_columns(T[j],CUTS)
        out,tm=analyse(cols,h,char,f"{D}_{h}_{j}_{char}_{sl}",slices=sl)
        print(D,h,j,["".join(map(str,a)) for a,e in cols],f"{tm:.1f}s",out,flush=True)
