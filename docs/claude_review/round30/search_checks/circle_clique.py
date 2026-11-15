# unit-circle locus rho_k=e^{i theta_k}, 0<theta_1<...<theta_D<pi. Im rho^m = sin(theta) U_{m-1}(cos theta):
# a difference delta must give a signed measure sum delta_k sin(theta_k) delta_{cos theta_k} with L vanishing moments
# (j=0..L-1), hence >= L sign changes in the order of cos theta (T-system / Descartes). Max clique with
# balance + distance + sign changes.
import igraph as ig, sys, time
D=int(sys.argv[1]); c=int(sys.argv[2]); K=(D+2*c)//2; L=D-1-K; h=(L+1)//2
def sc(x,y):
    seq=[(1 if ((x>>k)&1) else -1) - (1 if ((y>>k)&1) else -1) for k in range(D)]
    seq=[s for s in seq if s]; return sum(1 for i in range(len(seq)-1) if seq[i]*seq[i+1]<0)
E=[]
pc=[bin(v).count('1') for v in range(2**D)]
for x in range(2**D):
    for y in range(x+1,2**D):
        p=pc[x&~y]; n=pc[y&~x]
        if p>=h and n>=h and p+n>=L+1 and sc(x,y)>=L: E.append((x,y))
g=ig.Graph(n=2**D,edges=E); t=time.time()
print("unit circle D",D,"c",c,"L",L,"max clique",g.clique_number(),"time %.1f"%(time.time()-t),flush=True)
