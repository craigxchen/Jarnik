"""Direct check of (VI.4)-(VI.5), the descent to F_9 and the norm identity for q=27,k=2 at X=t^L:
F_j(t,t^L) = prod_{v in V} (t^L - v(t))^{1+chi(eps lambda_j(v))} computed in F_729[t] factor by factor."""
import numpy as np, itertools, sys, time
from bigfield import BF
from f9lib import F9P, pw, gcd9, divmod9
K = BF(6); q = 27; Q = 729
g = 1  # find generator element: exp[1]
gen = int(K.exp[1])
I = K.pow(gen, (Q-1)//4)                     # i^2 = -1
assert K.mul(I, I) == K.neg[1]
F27 = [0] + [K.pow(gen, 28*k) for k in range(26)]
assert len(set(F27)) == 27
chi27 = {0: 0}
for x in F27[1:]:
    chi27[x] = 1 if K.pow(x, 13) == 1 else -1
# F_9 inside F_729 in the basis 1, i
F9map = {}
one = 1
for a in range(3):
    for b in range(3):
        A = 0
        for _ in range(a): A = int(K.add[A, one])
        B = 0
        for _ in range(b): B = int(K.add[B, I])
        F9map[int(K.add[A, B])] = (a, b)
assert len(F9map) == 9
def lam_val(lam, x):
    s = 0
    for c, xb in zip(lam, x):
        t = xb if c % 3 == 1 else (K.neg[xb] if c % 3 == 2 else 0)
        s = int(K.add[s, t])
    return s

def product_point(lam, eps, L, k=2):
    poly = np.array([1], np.int64)   # coefficients in F_729, low first
    for x in itertools.product(F27, repeat=k):
        e = 1 + eps*chi27[lam_val(lam, x)]
        # v(t) = i*x0 + i*x1 t + ...
        vcoef = [K.mul(I, xb) for xb in x]
        for _ in range(e):
            # poly * (t^L - v)
            new = np.zeros(len(poly) + L, np.int64)
            new[L:] = poly
            for b, vb in enumerate(vcoef):
                if vb == 0: continue
                term = K.neg[K.mul(poly, vb)]
                seg = new[b:b+len(poly)]
                new[b:b+len(poly)] = K.add[seg, term]
            poly = new
    return poly

def to_f9(poly):
    re = []; im = []
    for c in poly:
        if int(c) not in F9map: return None
        a, b = F9map[int(c)]; re.append(a); im.append(b)
    return F9P(re, im)

if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    forms = [(1,0),(0,1),(1,1),(1,2)]
    pts = []
    t0 = time.time()
    for lam in forms:
        p = product_point(lam, 1, L)
        z = to_f9(p)
        assert z is not None, "coefficient outside F_9"
        pts.append(z)
    print('raw points computed in F_729[t]; all coefficients lie in F_9. deg', [z.deg() for z in pts], time.time()-t0)
    # compare with closed form (VI.5) computed in F_9[t]
    t = F9P([0, 1]); X = F9P([0]*L + [1])
    for (a, b), z in zip(forms, pts):
        w = F9P([-b, a]) * (0, 1)            # i*(-b + a t) spans ker of a x + b y
        wq1 = pw(w, q-1)
        A = pw(X, q) - wq1 * X
        # v0 with lambda(v0)=1: if a != 0: v0 = i*a^{-1}, else i*b^{-1} t
        v0 = F9P([pow(a, 1, 3) if a else 0, 0 if a else b % 3]) * (0, 1)  # a,b in {1,2}: inverse = itself mod 3
        c = pw(v0, q) - wq1 * v0
        F = A * pw(pw(A, (q-1)//2) - pw(c, (q-1)//2), 2)
        F2 = pw(X, 0)  # placeholder
        PV = pw(A, q) - pw(c, q-1) * A
        F3 = PV + pw(c, q-1) * A * (2, 0) - pw(c, (q-1)//2) * pw(A, (q+1)//2) * (2, 0)
        assert F == z and F3 == z
    print('closed form (VI.5) agrees with the direct product for all four marks')
    # norms
    PVs = set()
    norms = [z * z.conj() for z in pts]
    assert all(nm == norms[0] for nm in norms) and norms[0].is_f3()
    print('raw norms equal and in F_3[t], deg', norms[0].deg())
    G = pts[0]
    for z in pts[1:]: G = gcd9(G, z)
    print('full gcd degree', G.deg(), ' formal g*L =', 92*L)
    red = [divmod9(z, G) for z in pts]
    assert all(r.iszero() for _, r in red)
    red = [qq for qq, _ in red]
    D = red[0].deg()
    nm = [z * z.conj() for z in red]
    assert all(x == nm[0] for x in nm) and nm[0].is_f3()
    diffs = [(red[a] - red[b]).deg() for a, b in itertools.combinations(range(4), 2)]
    print('reduced D', D, 'parent deg', nm[0].deg(), 'diff degs', diffs, 'D/2', D/2)
    print('lead coeffs', [z.lead() for z in red], 'G lead', G.lead())
    # conjugate / unit relations
    for a in range(4):
        for b in range(4):
            if a != b: assert not (red[a] == red[b].conj())
    print('no point is the conjugate of another; gcd of reduced points:', end=' ')
    g2 = red[0]
    for z in red[1:]: g2 = gcd9(g2, z)
    print(g2.deg())
    np.save('four_mark_points.npy', np.array([np.concatenate([z.re, z.im]) for z in red]))
    import pickle; pickle.dump({'red': [(z.re, z.im) for z in red], 'G': (G.re, G.im), 'N': (nm[0].re, nm[0].im)}, open(f'four_mark_L{L}.pkl', 'wb'))
