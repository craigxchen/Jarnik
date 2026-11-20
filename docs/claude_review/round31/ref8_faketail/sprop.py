# S-property sanity at h=2 (all primitive balanced v=e_a-e_b on I) and h=4 (sampled), q=103,127.
# Also: first index s<=h/2 with sum v_j j^s != 0 mod q for random balanced v (Newton claim), incl. extremal cases.
import random, math
def paley(q):
    sq={x*x%q for x in range(1,q)}
    chi=[0]+[1 if x in sq else -1 for x in range(1,q)]
    return [[chi[(j-i)%q] for j in range(q)] for i in range(q)]
for q in (103,127):
    random.seed(q); K=paley(q); A=[i for i in range(q) if K[i][0]==1]
    rho=0.9; S=set(i for i in A if random.random()<rho); Sc=[i for i in range(q) if i not in S]
    worstg=0; worstr=9
    def test(v):
        global worstg,worstr
        w=[sum(K[i][j]*x for j,x in v.items()) for i in Sc]
        g=0
        for x in w: g=math.gcd(g,x)
        n2=sum(x*x for x in v.values()); r=sum(x*x for x in w)/((1-rho)*q*n2/2)
        worstg=max(worstg,g); worstr=min(worstr,r)
    for a in range(1,q):
        for b in range(1,q):
            if a!=b: test({a:1,b:-1})
    for _ in range(20000):
        lab=random.sample(range(1,q),4); 
        if random.random()<.5: test({lab[0]:1,lab[1]:1,lab[2]:-1,lab[3]:-1})
        else: test({lab[0]:2,lab[1]:-1,lab[2]:-1})
    print(f"q={q} rho={rho}: max gcd((Kv)_Sc)={worstg}, min ||(Kv)_Sc||^2/((1-rho)q||v||^2/2)={worstr:.3f}")
    # Newton
    mx=0
    for _ in range(3000):
        t=random.randint(1,(q//2-1)//2); 
        P=[random.randint(1,q-1) for _ in range(t)]
        rest=[x for x in range(1,q) if x not in P]; N=[random.choice(rest) for _ in range(t)]
        v={}
        for x in P: v[x]=v.get(x,0)+1
        for x in N: v[x]=v.get(x,0)-1
        s=next(s for s in range(1,q) if sum(c*pow(j,s,q) for j,c in v.items())%q)
        assert s<=t; mx=max(mx,s/t)
    # extremal: P = t-th roots-of-unity coset vs another coset have p_s equal for s<t
    t=next(t for t in range(2,q) if (q-1)%t==0 and t>2)
    g=next(g for g in range(2,q) if all(pow(g,(q-1)//r,q)!=1 for r in range(2,q) if (q-1)%r==0 and all(r%u for u in range(2,r))))
    H=[pow(g,(q-1)//t*i,q) for i in range(t)]
    P=H; N=[x*g%q for x in H]
    v={x:1 for x in P}; 
    for x in N: v[x]=v.get(x,0)-1
    s=next(s for s in range(1,q) if sum(c*pow(j,s,q) for j,c in v.items())%q)
    print(f"  Newton: random max s/t={mx:.2f}; coset example t={t}: first nonzero moment s={s} (=t: bound attained)")
