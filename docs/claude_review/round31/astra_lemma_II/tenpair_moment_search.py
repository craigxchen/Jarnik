"""Nonexhaustive real search for the degree-one five-star system."""
import argparse
import itertools
import numpy as np
from scipy.optimize import least_squares

pairs = list(itertools.combinations(range(5), 2))
incidence = np.array([[float(i in edge) - float(4 in edge) for edge in pairs] for i in range(4)])
orders = np.arange(1, 6)[:, None]

def evaluate(v):
    z = v[:10] + 1j*v[10:]
    data = z[None, :] ** orders
    moments = (data.imag @ incidence.T) / orders
    return np.r_[moments.ravel(), v[:10].sum(), np.dot(v, v) - 10]

def jacobian(v):
    z = v[:10] + 1j*v[10:]
    partial = z[None, :] ** (orders - 1)
    dx = incidence[None, :, :] * partial.imag[:, None, :]
    dy = incidence[None, :, :] * partial.real[:, None, :]
    return np.vstack((np.c_[dx.reshape(-1,10),dy.reshape(-1,10)], np.r_[np.ones(10),np.zeros(10)], 2*v))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--seed', type=int, default=7)
    args=parser.parse_args()
    rng=np.random.default_rng(args.seed)
    best=0
    near=0
    for n in range(args.count):
        v=rng.normal(size=20)
        v[:10]-=v[:10].mean()
        v*=np.sqrt(10/np.dot(v,v))
        fit=least_squares(evaluate,v,jac=jacobian,method='lm',max_nfev=350,ftol=1e-12,xtol=1e-12,gtol=1e-12)
        err=np.linalg.norm(fit.fun,np.inf)
        if err>1e-9:
            continue
        near+=1
        z=fit.x[:10]+1j*fit.x[10:]
        roots=np.r_[z,z.conj()]
        sep=min(abs(a-b) for a,b in itertools.combinations(roots,2))
        if sep>best:
            best=sep
            if sep>.001:
                print(n,err,sep,fit.x.tolist(),flush=True)
    print("starts",args.count,"small residual",near,"best minimum separation",best)
if __name__=='__main__':
    main()
