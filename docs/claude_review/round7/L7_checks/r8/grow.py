# Build members of Mobius type E by adding events one at a time (continuation).
import numpy as np, itertools, sys
from fwd import *
def pair_collisions(g, a, b):
    A, B = g[a], g[b]
    c2 = A[0,0]*B[1,0] - B[0,0]*A[1,0]; c1 = A[0,0]*B[1,1] + A[0,1]*B[1,0] - B[0,0]*A[1,1] - B[0,1]*A[1,0]
    c0 = A[0,1]*B[1,1] - B[0,1]*A[1,1]
    return np.roots([c2, c1, c0])
def mob(gm, t): return (gm[0,0]*t + gm[0,1])/(gm[1,0]*t + gm[1,1])
def grow(E, rng, tries=30, verbose=False):
    n = 7
    order = list(E); rng.shuffle(order)
    # start: random member with first event
    cur = [order[0]]
    R = setup(cur, rng)
    for _ in range(50):
        V, ok = gn(rng.normal(size=2+28)+1j*rng.normal(size=2+28), cur, R)
        if ok and check(V, cur)[1] == 'ok': break
    else: return None
    for k in range(1, len(order)):
        e = order[k]; done = False
        m = len(cur); s = V[:m]; z = V[m:2*m]; C = V[2*m:]
        g = [to_mobius(C.reshape(n,4)[a]) for a in range(n)]
        cands = []
        for (a, b) in itertools.combinations(e, 2):
            for t in pair_collisions(g, a, b):
                if not any(abs(t - s[i]) < 1e-6 for i in range(m)): cands.append(t)
        rng.shuffle(cands)
        for t0 in cands[:tries]:
            z0 = mob(g[e[0]], t0)
            V0 = np.concatenate([np.append(s, t0), np.append(z, z0), C])
            new = cur + [e]
            Vn, ok = gn(V0, new, R)
            if ok and check(Vn, new)[1] == 'ok':
                V = Vn; cur = new; done = True; break
        if not done:
            if verbose: print('stuck at', k, e)
            return None
    return cur, V, R
if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]))
    types = [eval(l) for l in open('types12.txt')]
    idx = int(sys.argv[2])
    E = list(types[idx]); succ = 0
    for att in range(int(sys.argv[3]) if len(sys.argv)>3 else 20):
        r = grow(E, rng)
        if r is not None:
            cur, V, R = r; succ += 1
            np.save('m12_type%d_%d.npy' % (idx, succ), V)
            np.save('m12_type%d_%d_R.npy' % (idx, succ), np.array(R))
            open('m12_type%d_%d_E.txt' % (idx, succ), 'w').write(repr(cur))
            print('type', idx, 'success', succ, 'attempt', att, flush=True)
            if succ >= 2: break
    print('type', idx, 'successes', succ)
