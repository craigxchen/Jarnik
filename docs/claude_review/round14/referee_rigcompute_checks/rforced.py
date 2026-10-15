import sys, json, subprocess, time
from reflib import *
SETS={'first4':[5,13,17,29],'rand4':[53,61,89,173],'bal4':[97,101,109,113],'last4':[353,373,389,397],
      'first5':[5,13,17,29,37]}
def verify_free(ps,K,n,c,q):
    nu=nu_par(ps,q)
    if n%2==0: assert all(x%2==p for x,p in zip(c,nu))
    for v in monos(len(ps),K,half=True):
        assert sum(a*b for a,b in zip(v,c))%n!=0
def run(name,K,Qmax):
    ps=SETS[name]; r=len(ps)
    qs=[q for q in primerange(3,Qmax+1) if q not in ps]
    lines=[]
    for q in qs:
        n=m_of(q)//4; nu=nu_par(ps,q)
        lines.append(f"{r} {K} {n} "+" ".join(map(str,nu)))
    out=subprocess.run(["./rfree"],input="\n".join(lines)+"\n",capture_output=True,text=True).stdout.split("\n")
    forced=[];ab=[]
    for q,ln in zip(qs,out):
        t=ln.split()
        if t[0]=="FREE": verify_free(ps,K,m_of(q)//4,list(map(int,t[1:])),q)
        elif t[0]=="FORCED": forced.append(q)
        else: ab.append(q)
    return forced,ab
if __name__=="__main__":
    name,K=sys.argv[1],int(sys.argv[2])
    a=json.load(open(f"../rigcompute_checks/thr_{name}_K{K}.json"))
    t=time.time()
    f,ab=run(name,K,a['Qscan'])
    print(name,K,"Qscan",a['Qscan'],"forced",len(f),"max",max(f),"abort",ab,"agree",f==sorted(a['forced']),f"[{time.time()-t:.0f}s]",flush=True)
    if f!=sorted(a['forced']): print(" mine-only",sorted(set(f)-set(a['forced'])),"theirs-only",sorted(set(a['forced'])-set(f)))
    json.dump(dict(forced=f,abort=ab),open(f"rforced_{name}_K{K}.json","w"))
