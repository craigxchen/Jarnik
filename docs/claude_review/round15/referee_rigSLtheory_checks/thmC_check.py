# Test the "second modulus" reading of Theorem C on the A1 labels (one-modulus, first ones found).
import sys, json, rlib, onemod
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); F=json.load(open(sys.argv[3]))['F']; q0=int(sys.argv[4])
pis,V,data=onemod.setup(P,K,F,q0); r=len(P)
m,L=rlib.labels(pis,q0); n=m//4
res=onemod.run(P,K,F,q0,cap=3000,detail=False)
# need PU-outside labels: rerun enum PU and keep labels
consPU=[]
for v in V:
    D,N,g=data[v]; mk=0
    for s in range(4):
        if q0*q0*g[s]*g[s] > (N if s%2==0 else 2*N): mk|=1<<s
    consPU.append((v,mk))
out=rlib.run_enum(r,m,1,3000,[x%2 for x in L],consPU,rlib.classcodes(L,n))
labs=[tuple(map(int,l.split()[1:])) for l in out.splitlines() if l.startswith('L ')]
cnt=0
for c in labs:
    for v in V:
        s=sum(a*b for a,b in zip(v,c))%m
        if s%n: continue
        sp=s//n; D,N,g=data[v]; o=(sp+2)%4
        UF1=[abs(D[t])//g[t] for t in range(4)]          # cofactor rel F' (includes q0-part)
        # cofactor relative to all of F: remove q0-part too
        def strip(x):
            while x%q0==0: x//=q0
            return x
        UF=[strip(u) for u in UF1]
        Es = q0 if D[sp]%q0 else 1     # extra at q0 for unit sp (level 1 only)
        delta_o = q0**rlib.ordp(D[o],q0)  # dropped genuine q0-part of opposite unit (label collides at sp, not o)
        if UF[o]==1 and Es>UF[sp] and q0>UF1[sp]:
            mm=list(g); mm[sp]*=q0
            ok=onemod.J0ok(mm,N) is None
            if ok and delta_o>1:
                cnt+=1
                if cnt<=5: print('label',c,'v',v,'s',sp,'U^F_s',UF[sp],'U^F_opp',UF[o],'q0-part of D_opp',delta_o,'J0 ok -> no second modulus needed')
print('instances (v-level) where Thm C hypotheses hold, J0 passes, and the only drop is at q0 itself:',cnt)
