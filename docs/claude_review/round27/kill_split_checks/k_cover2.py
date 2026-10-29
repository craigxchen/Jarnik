# Longer, varied search at a given one-sided class (by size signature), reporting residual stagnation statistics.
import numpy as np, itertools, sys, time
sys.path.insert(0,'.')
_argv=sys.argv; sys.argv=['x','0']
exec(open('k_cover.py').read().split("classes=onesided_classes()")[0])
sys.argv=_argv
sig=eval(sys.argv[1]); T=float(sys.argv[2]); seed=int(sys.argv[3])
rng=np.random.default_rng(seed)
cls=[F for F in onesided_classes() if sorted(len(B) for B in F)==sig]
print("matching classes:",len(cls)); F,A,D=setup(cls[0])
t0=time.time(); stats=[]; found=0; tries=0
while time.time()-t0<T:
    tries+=1
    scale=10**rng.uniform(-0.7,0.9)
    s4=complex(*rng.normal(size=2))*10**rng.uniform(-0.5,0.8)
    if rng.random()<0.3:
        z=np.concatenate([rng.normal(size=15), rng.normal(size=4)]).astype(complex)*scale
        z+=1e-3j*rng.normal(size=19)
    else:
        z=np.concatenate([rng.normal(size=15)+1j*rng.normal(size=15), rng.normal(size=4)+1j*rng.normal(size=4)])*scale
    z,nf,it=newton(z,s4,A,D,it=300)
    # normalised residual
    stats.append(nf)
    if nf<1e-11:
        r=z[:15]; S=np.array([0,1,-1,s4])
        sep=min(min(abs(r[i]-r[j]) for i in range(15) for j in range(i)), min(abs(si-rj) for si in S for rj in r))
        cdeg=min(abs(z[15:19])); c1=min(abs(z[15:19]-1))
        print("  converged: sep=%.2e |c|min=%.2e |c-1|min=%.2e |r|max=%.2e"%(sep,cdeg,c1,max(abs(r)))); sys.stdout.flush()
        if sep>1e-4 and cdeg>1e-4 and c1>1e-4:
            found+=1
            sdom,sv=dominance_rank(z,s4,F,A,D,complex(0.37,0.21))
            print("  NONDEGENERATE solution; dominance SVs",sdom); np.save('k_cover2_sol_%d.npy'%seed,np.append(z,s4)); break
stats=np.array(stats)
print("tries",tries,"found",found,"residual quantiles",np.quantile(stats,[0,0.1,0.5,0.9]))
