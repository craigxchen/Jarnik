"""NUMERICAL sanity test: Klein ansatz for class-beta rational curves in Mbar_{0,6}.
Anchors 0,1,inf; movers a = A(t+1/t) (deg 2), b = B(t^2) (deg 2), c = C(t^2+t^-2) (Moebius).
Expected dimension 3 (fix the three z-values C^{-1}(0), C^{-1}(1), C^{-1}(inf))."""
import numpy as np, sys, itertools

def build(u, zc):
    v1, y1, v2, y2 = u[0:3], u[3:6], u[6:9], u[9:12]
    kA, kB, v0, y0 = u[12], u[13], u[14], u[15]
    return v1, y1, v2, y2, kA, kB, v0, y0

def Cmap(z, zc):
    return (z - zc[0])*(zc[1] - zc[2])/((z - zc[2])*(zc[1] - zc[0]))

def residual(u, zc):
    v1, y1, v2, y2, kA, kB, v0, y0 = build(u, zc)
    R = []
    R += list(v1**2 - 2 - zc)
    R += list(y1 + 1/y1 - zc)
    R += list(v2**2 - 2 - y2 - 1/y2)
    PA = [np.poly([v1[x], v2[x]]) for x in range(3)]
    PB = [np.poly([y1[x], y2[x]]) for x in range(3)]
    depA = kA*PA[0] - PA[2] - (kA-1)*PA[1]
    depB = kB*PB[0] - PB[2] - (kB-1)*PB[1]
    R += list(depA[1:]) + list(depB[1:])
    A = kA*np.polyval(PA[0], v0)/np.polyval(PA[2], v0)
    B = kB*np.polyval(PB[0], y0)/np.polyval(PB[2], y0)
    z0 = v0**2 - 2
    R += [z0 - y0 - 1/y0, A - Cmap(z0, zc), B - Cmap(z0, zc)]
    return np.array(R)

def newton(u, zc, iters=100):
    for it in range(iters):
        F = residual(u, zc)
        nF = np.linalg.norm(F)
        if not np.isfinite(nF): return None
        if nF < 1e-12: return u
        J = np.zeros((16, 16), complex); h = 1e-7
        for k in range(16):
            du = np.zeros(16, complex); du[k] = h
            J[:, k] = (residual(u+du, zc) - residual(u-du, zc))/(2*h)
        try: du = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError: return None
        st = 1.0
        while st > 1e-8:
            u1 = u + st*du; F1 = residual(u1, zc)
            if np.all(np.isfinite(F1)) and np.linalg.norm(F1) < nF*(1-0.05*st): break
            st /= 2
        if st <= 1e-8: return None
        u = u1
    return None

def movers_t(u, zc):
    v1, y1, v2, y2, kA, kB, v0, y0 = build(u, zc)
    PA = [np.poly([v1[x], v2[x]]) for x in range(3)]
    PB = [np.poly([y1[x], y2[x]]) for x in range(3)]
    def sub_v(P):
        deg = len(P)-1; out = np.zeros(1, complex)
        for k, co in enumerate(P):
            term = np.array([1.0+0j])
            for _ in range(deg-k): term = np.polymul(term, [1, 0, 1])
            term = np.polymul(term, np.concatenate([[1.0], np.zeros(k)]))
            out = np.polyadd(out, co*term)
        return out
    def sub_y(P):
        out = np.zeros(1, complex)
        for co in P: out = np.polyadd(np.polymul(out, [1, 0, 0]), [co])
        return out
    # c = C(z) with z = t^2+t^-2: numerator (z-zc0)(zc1-zc2), denominator (z-zc2)(zc1-zc0); times t^2
    zpoly = lambda r: np.array([1, 0, -r, 0, 1], complex)  # t^2*(z - r) = t^4 - r t^2 + 1
    a = (kA*sub_v(PA[0]), sub_v(PA[2]))
    b = (kB*sub_y(PB[0]), sub_y(PB[2]))
    c = ((zc[1]-zc[2])*zpoly(zc[0]), (zc[1]-zc[0])*zpoly(zc[2]))
    return {'a': a, 'b': b, 'c': c}

def collisions(M, names_extra=None, tol=1e-6):
    pts = {'0': (np.array([0j]), np.array([1+0j])), '1': (np.array([1+0j]), np.array([1+0j])),
           'i': (np.array([1+0j]), np.array([0j]))}
    pts.update(M)
    names = list(pts)
    cand = []
    for p, q in itertools.combinations(names, 2):
        N = np.polysub(np.polymul(pts[p][0], pts[q][1]), np.polymul(pts[q][0], pts[p][1]))
        while len(N) > 1 and abs(N[0]) < 1e-10*np.max(np.abs(N)): N = N[1:]
        cand += list(np.roots(N))
    reps = []
    for x in cand:
        if not any(abs(x-r) < tol*max(1, abs(r)) for r in reps): reps.append(x)
    events = []
    for t0 in reps:
        vals = {p: (np.polyval(pts[p][0], t0), np.polyval(pts[p][1], t0)) for p in names}
        groups = []
        for p in names:
            for g in groups:
                n1, d1 = vals[p]; n2, d2 = vals[g[0]]
                sc = max(abs(n1), abs(d1))*max(abs(n2), abs(d2))
                if abs(n1*d2-n2*d1) < 1e-5*sc: g.append(p); break
            else: groups.append([p])
        events.append((t0, [g for g in groups if len(g) >= 2]))
    return names, events

if __name__ == '__main__':
    from collections import Counter
    nst = int(sys.argv[1])
    rng0 = np.random.default_rng(12345)
    zc = rng0.normal(size=3) + 1j*rng0.normal(size=3)
    nsol = 0
    for seed in range(nst):
        rng = np.random.default_rng(seed)
        u = rng.normal(size=16) + 1j*rng.normal(size=16)
        u = newton(u, zc)
        if u is None: continue
        names, ev = collisions(movers_t(u, zc))
        cnt = Counter(); bad = 0
        P = set(names)
        for t0, big in ev:
            if len(big) != 1 or len(big[0]) > len(P)-2: bad += 1; continue
            S = tuple(sorted(big[0])); Sc = tuple(sorted(P - set(S)))
            cnt[min(S, Sc)] += 1
        print(seed, 'fibres', len(ev), 'bad', bad, 'splits', len(cnt), 'maxmult', max(cnt.values()) if cnt else 0, flush=True)
        if bad == 0 and len(cnt) == 25 and max(cnt.values()) == 1:
            np.save('klein6_sol_%d.npy' % seed, u); nsol += 1
    np.save('klein6_zc.npy', zc)
    print('nondegenerate', nsol)
