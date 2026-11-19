import subprocess, itertools, random, time, re
def system(columns,h,slices=0,seed=1,npts=5,xy=True):
    n=len(columns)
    X=[f"x{k}" for k in range(n)]; Y=[f"y{k}" for k in range(n)]
    def rms(k,m):
        # r^m - s^m with r=x+y, s=x-y : 2*sum_{j odd} C(m,j) x^{m-j} y^j
        from math import comb
        t=[]
        for j in range(1,m+1,2):
            c=2*comb(m,j)
            t.append(f"{c}*{X[k]}^{m-j}*{Y[k]}^{j}" if m-j>0 else f"{c}*{Y[k]}^{j}")
        return "("+"+".join(t)+")"
    eqs=[]
    for i in range(1,npts):
        for m in range(1,h):
            from math import comb
            terms=[]
            for k,(a,e) in enumerate(columns):
                c=a[i]-a[0]
                if not c: continue
                for j in range(1,m+1,2):
                    cc=2*comb(m,j)*c
                    terms.append(f"{cc}*{X[k]}^{m-j}*{Y[k]}^{j}" if m-j>0 else f"{cc}*{Y[k]}^{j}")
            if terms: eqs.append("+".join(terms).replace("+-","-"))
    eqs.append("+".join(f"{e}*{X[k]}" for k,(a,e) in enumerate(columns)))
    eqs.append(f"{Y[0]}-1")
    rnd=random.Random(seed)
    for _ in range(slices):
        eqs.append("+".join(f"{rnd.randint(1,10**4)}*{v}" for v in X+Y[1:])+f"+{rnd.randint(1,10**4)}")
    dup=[(k,l) for k,l in itertools.combinations(range(n),2) if columns[k]==columns[l]]
    fac=list(Y)
    for k,l in dup:
        lam=rnd.randint(2,10**4)
        fac.append(f"(x{k}-x{l}+{lam}*y{k}+{lam}*y{l})")
    if dup:
        vars_extra=["u"]
    import sympy
    ue=sympy.expand(sympy.sympify("u*"+"*".join(fac)+"-1"))
    eqs.append(str(ue).replace("**","^").replace(" ",""))
    rnd2=random.Random(seed+12345)
    eqs.append("T-"+"-".join(f"{rnd2.randint(1,997)}*{v}" for v in X+Y+["u"]))
    return X+Y+["u","T"],eqs
def forms(n):
    # degeneracies (complexified r=x+y,s=x-y): r_k=r_l: x_k-x_l+y_k-y_l ; s_k=s_l ; r_k=s_l: x_k-x_l+y_k+y_l ; s_k=r_l
    F=[]
    for k,l in itertools.combinations(range(n),2):
        F+=[f"x{k}-x{l}+y{k}-y{l}",f"x{k}-x{l}-y{k}+y{l}",f"x{k}-x{l}+y{k}+y{l}",f"x{k}-x{l}-y{k}-y{l}"]
    return F
def singular_script(vars_,eqs,F,char):
    return f'''LIB "elim.lib";
option(redSB);
ring rr={char},({",".join(vars_)}),dp;
ideal I={",".join(eqs)};
I=std(I);
print("before: dim="+string(dim(I))+" vdim="+string(vdim(I)));
list F={",".join(F)};
int j;
for(j=1;j<=size(F);j++){{ if (deg(I[1])==0) {{break;}} I=sat(I,F[j]); I=std(I); }}
print("after: dim="+string(dim(I))+" vdim="+string(vdim(I))+" one="+string(deg(I[1])==0));
quit;
'''
def run_sing(columns,h,char,tag,slices=0,seed=1,timeout=3600,npts=5):
    v,e=system(columns,h,slices,seed,npts)
    scr=singular_script(v,e,forms(len(columns)),char)
    fn=f"/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round31/ref7L_search/s_{tag}.sing"
    open(fn,'w').write(scr)
    t=time.time()
    try:
        out=subprocess.run(["prlimit","--as=3500000000","Singular","-q",fn],capture_output=True,text=True,timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        out="TIMEOUT"
    return out.strip().replace("\n"," | "),time.time()-t
