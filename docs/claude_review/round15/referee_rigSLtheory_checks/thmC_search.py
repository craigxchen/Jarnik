# Per-monomial search: v in V_K, q0 in F, unit s with q0 | D_{s+2}(Z_v) (genuine collision of opposite unit at q0),
# q0 not | D_s, U^F_{s+2}=1 (opposite unit saturated rel. all of F), q0 > U^F_s (extra beyond cofactor slack).
# Demand of a system that moves the q0-collision from unit s+2 to unit s and is genuine elsewhere:
# m_s = q0*g_s, m_{s+2} = g_{s+2}/q0^a, other units genuine. Check full J0 (Lemma 2.1) for v.
import sys,json,rlib,onemod
P=json.loads(sys.argv[1]); K=int(sys.argv[2]); F=json.load(open(sys.argv[3]))['F']
pis=[rlib.gpi(p) for p in P]; found=0; tested=0
for v in rlib.VK(len(P),K):
    z=rlib.Zv(pis,v); D=rlib.Ds(z); N=z[0]**2+z[1]**2
    gF=[]
    for s in range(4):
        g=1
        for q in F: g*=q**rlib.ordp(D[s],q)
        gF.append(g)
    UF=[abs(D[s])//gF[s] for s in range(4)]
    for s in range(4):
        o=(s+2)%4
        if UF[o]!=1: continue
        for q0 in F:
            if D[o]%q0 or D[s]%q0==0: continue
            if not q0>UF[s]: continue
            tested+=1
            a=rlib.ordp(D[o],q0)
            mm=list(gF); mm[s]=gF[s]*q0; mm[o]=gF[o]//q0**a
            if onemod.J0ok(mm,N) is None:
                found+=1
                if found<=5: print('v',v,'s',s,'q0',q0,'a',a,'U^F_s',UF[s],'J0 holds: (2.1) satisfied with the only drop at q0')
print('tested',tested,'J0-feasible with no second modulus:',found)
