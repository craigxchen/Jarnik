# Lemma 4.0 identity, own code, exact integers, q=3251 and small q; also h'(e_x)=q+1
import random
def run(q,trials):
    chi=[0]+[1 if pow(x,(q-1)//2,q)==1 else -1 for x in range(1,q)]
    Pp=[[a]+[x for x in range(q) if chi[(a-x)%q]==-1] for a in range(q)]
    Pm=[[a]+[x for x in range(q) if chi[(a-x)%q]==1] for a in range(q)]
    bad=0; minratio=None
    for t in range(trials):
        sup=random.sample(range(q),random.randint(1,min(8,q)))
        v=[0]*q
        for x in sup: v[x]=random.choice([-3,-2,-1,1,2,3])
        K=sum(v); n2=sum(x*x for x in v); n1=sum(abs(x) for x in v)
        S=[sum(v[x] for x in P) for P in Pp]; Sm=[sum(v[x] for x in P) for P in Pm]
        if 4*sum(s*s for s in S)!=(q+1)*(K*K+n2) or 4*sum(s*s for s in Sm)!=(q+1)*(K*K+n2): bad+=1
        h=sum(abs(s) for s in S)+sum(abs(s) for s in Sm)
        if 2*h*n1 < (q+1)*n2: bad+=1
        r=h/(q+1); minratio=r if minratio is None else min(minratio,r)
    e=[0]*q; e[5]=1
    he=sum(abs(sum(e[x] for x in P)) for P in Pp+Pm)
    print('q',q,'trials',trials,'identity/ineq failures',bad,"min h'/(q+1) sampled",minratio,"h'(e_x)=",he)
random.seed(7)
for q in (7,11,19,23,43,3251): run(q,300 if q==3251 else 2000)
