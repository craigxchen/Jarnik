# General graph ansatz: label a moves by x_a(t) = N_a(t)/M_a(t), deg <= d_a. Events e (tuples of labels)
# at (s_e, z_e). Incidence: N_a(s) - z M_a(s) = 0. Normalisation: random linear form on coeffs = 1.
import numpy as np, itertools
class GraphProblem:
    def __init__(self, degs, E, rng):
        self.d = list(degs); self.E = [tuple(e) for e in E]; self.n = len(degs)
        self.off = np.cumsum([0] + [2*(d+1) for d in self.d])
        self.R = [rng.normal(size=2*(d+1)) + 1j*rng.normal(size=2*(d+1)) for d in self.d]
        self.m = len(self.E); self.N = 2*self.m + self.off[-1]
    def coeffs(self, V, a):
        c = V[2*self.m + self.off[a]: 2*self.m + self.off[a+1]]; d = self.d[a]
        return c[:d+1], c[d+1:]
    def F(self, V):
        m = self.m; out = []
        for i, e in enumerate(self.E):
            s, z = V[i], V[m+i]
            for a in e:
                Nc, Mc = self.coeffs(V, a); pw = s**np.arange(len(Nc))
                out.append(Nc @ pw - z*(Mc @ pw))
        for a in range(self.n):
            c = V[2*m + self.off[a]: 2*m + self.off[a+1]]; out.append(self.R[a] @ c - 1)
        return np.array(out)
    def J(self, V):
        m = self.m; rows = []
        for i, e in enumerate(self.E):
            s, z = V[i], V[m+i]
            for a in e:
                Nc, Mc = self.coeffs(V, a); d = self.d[a]; k = np.arange(d+1); pw = s**k
                dpw = np.concatenate([[0], k[1:]*s**(k[1:]-1)])
                r = np.zeros(self.N, complex)
                r[i] = Nc @ dpw - z*(Mc @ dpw); r[m+i] = -(Mc @ pw)
                base = 2*m + self.off[a]; r[base:base+d+1] = pw; r[base+d+1:base+2*d+2] = -z*pw
                rows.append(r)
        for a in range(self.n):
            r = np.zeros(self.N, complex); r[2*m+self.off[a]:2*m+self.off[a+1]] = self.R[a]; rows.append(r)
        return np.array(rows)
    def gn(self, V, maxit=80, tol=1e-12):
        for it in range(maxit):
            f = self.F(V)
            if not np.all(np.isfinite(f)): return V, False
            if np.linalg.norm(f) < tol: return V, True
            V = V + np.linalg.lstsq(self.J(V), -f, rcond=None)[0]
            if np.linalg.norm(V) > 1e7: return V, False
        return V, np.linalg.norm(self.F(V)) < 1e-9
    def val(self, V, a, t):
        Nc, Mc = self.coeffs(V, a); pw = t**np.arange(len(Nc)); return (Nc @ pw)/(Mc @ pw)
    def pairpoly(self, V, a, b):
        Na, Ma = self.coeffs(V, a); Nb, Mb = self.coeffs(V, b)
        P = np.polynomial.polynomial
        return P.polysub(P.polymul(Na, Mb), P.polymul(Nb, Ma))
    def contacts(self, V, tol=1e-6):
        """List of (t, cluster) for all contacts; flag. Requires every pair poly to have full degree
        d_a+d_b (no contact at t=inf), every event a genuine |e|-cluster, all contact times distinct."""
        m = self.m; n = self.n
        for a in range(n):
            Nc, Mc = self.coeffs(V, a)
            # degree check: resultant nonzero ~ check that N, M have no common root
            rN = np.roots(Nc[::-1]) if self.d[a] > 0 else []
            if self.d[a] > 0 and len(rN) and min(abs(np.polyval(Mc[::-1], r)) for r in rN) < 1e-8*np.linalg.norm(Mc):
                return None, 'common_root'
        times = [(V[i], frozenset(e)) for i, e in enumerate(self.E)]
        for a, b in itertools.combinations(range(n), 2):
            pp = self.pairpoly(V, a, b)
            if abs(pp[-1]) < 1e-9*np.linalg.norm(pp): return None, 'deg_drop'
            for t in np.roots(pp[::-1]):
                if not any(abs(t - u) < tol*(1+abs(u)) and {a, b} <= S for (u, S) in times[:m]):
                    times.append((t, frozenset((a, b))))
        us = [t for t, S in times]
        for i in range(len(us)):
            for j in range(i+1, len(us)):
                if abs(us[i] - us[j]) < tol*(1+abs(us[i])): return None, 'coincide'
        for (t, S) in times:
            pos = [self.val(V, a, t) for a in range(n)]
            cl = [pos[a] for a in S]
            if max(abs(c - cl[0]) for c in cl) > 1e-5*(1+abs(cl[0])): return None, 'notcluster'
            others = [pos[a] for a in range(n) if a not in S] + [cl[0]]
            for i in range(len(others)):
                for j in range(i+1, len(others)):
                    if abs(others[i] - others[j]) < tol*(1+abs(others[i])): return None, 'extra'
        return times, 'ok'
def grow(degs, Efull, rng, tries=40, verbose=False, start=None):
    n = len(degs); order = list(Efull); rng.shuffle(order)
    cur = []; P = GraphProblem(degs, cur, rng)
    for _ in range(50):
        V, ok = P.gn(rng.normal(size=P.N) + 1j*rng.normal(size=P.N))
        if ok and P.contacts(V)[1] == 'ok': break
    else: return None
    for k, e in enumerate(order):
        done = False
        Pn = GraphProblem(degs, cur + [e], rng); Pn.R = P.R
        cands = []
        for (a, b) in itertools.combinations(e, 2):
            pp = P.pairpoly(V, a, b)
            for t in np.roots(pp[::-1]):
                if not any(abs(t - V[i]) < 1e-6 for i in range(P.m)): cands.append(t)
        rng.shuffle(cands)
        for t0 in cands[:tries]:
            z0 = P.val(V, e[0], t0)
            V0 = np.concatenate([np.append(V[:P.m], t0), np.append(V[P.m:2*P.m], z0), V[2*P.m:]])
            Vn, ok = Pn.gn(V0)
            if ok and Pn.contacts(Vn)[1] == 'ok':
                V = Vn; cur = cur + [e]; P = Pn; done = True; break
        if not done:
            if verbose: print('stuck at', k, e, flush=True)
            return (cur, V, P, False)
        if verbose and k % 5 == 0: print('added', k+1, flush=True)
    return (cur, V, P, True)

def try_add(P, V, cur, e, degs, rng, tries=40):
    Pn = GraphProblem(degs, cur + [e], rng); Pn.R = P.R
    cands = []
    for (a, b) in itertools.combinations(e, 2):
        pp = P.pairpoly(V, a, b)
        for t in np.roots(pp[::-1]):
            if not any(abs(t - V[i]) < 1e-6 for i in range(P.m)): cands.append(t)
    rng.shuffle(cands)
    for t0 in cands[:tries]:
        z0 = P.val(V, e[0], t0)
        V0 = np.concatenate([np.append(V[:P.m], t0), np.append(V[P.m:2*P.m], z0), V[2*P.m:]])
        Vn, ok = Pn.gn(V0)
        if ok and Pn.contacts(Vn)[1] == 'ok': return Pn, Vn
    return None
def grow2(degs, Efull, rng, verbose=False, init=None):
    """dynamic order: at each step try remaining events in random order until one can be added."""
    if init is None:
        cur = []; P = GraphProblem(degs, cur, rng)
        for _ in range(50):
            V, ok = P.gn(rng.normal(size=P.N) + 1j*rng.normal(size=P.N))
            if ok and P.contacts(V)[1] == 'ok': break
        else: return None
    else:
        cur, V, P = init
    rem = [e for e in Efull if e not in cur]
    while rem:
        order = list(rem); rng.shuffle(order); added = False
        for e in order:
            r = try_add(P, V, cur, e, degs, rng)
            if r is not None:
                P, V = r; cur = cur + [e]; rem.remove(e); added = True
                if verbose: print('added', len(cur), e, flush=True)
                break
        if not added:
            if verbose: print('stuck with', len(rem), 'remaining:', rem, flush=True)
            return (cur, V, P, False)
    return (cur, V, P, True)
