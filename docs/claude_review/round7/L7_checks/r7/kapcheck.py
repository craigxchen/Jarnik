"""NUMERICAL genuineness check of a cubic curve c(t) in the Kapranov psi_7 model of Mbar_{0,7}.
P^4 coords (c0..c4); p1..p5 = e0..e4, p6 = (1,1,1,1,1).  Labels 1..6 <-> points p_a.
Returns the dict of contacts: plane contacts, line/point incidences, hyperplane (pair) contacts."""
import numpy as np, itertools
Pn = np.polynomial.polynomial
PTS = {a: np.eye(5)[a-1] for a in range(1, 6)}
PTS[6] = np.ones(5)

def annihilator(span):
    """rows = linear forms vanishing on span of given vectors"""
    A = np.array(span)
    u, s, vh = np.linalg.svd(A)
    r = len(span)
    return vh[r:]

def common_roots(polys, tol=1e-7):
    """roots t (finite, plus inf) where all polynomials (low-first coeff arrays, deg<=3) vanish"""
    cand = []
    for p in polys:
        q = np.trim_zeros(p, 'b')
        if len(q) > 1:
            cand += list(Pn.polyroots(q))
    out = []
    scale = max(np.max(np.abs(p)) for p in polys) + 1e-300
    for r in cand:
        if all(abs(Pn.polyval(r, p)) < tol*scale*(1+abs(r))**3 for p in polys):
            if all(abs(r - o) > 1e-6 for o in out):
                out.append(r)
    # infinity: leading (t^3) coefficients all vanish
    if all(abs(p[3] if len(p) > 3 else 0) < tol*scale for p in polys):
        out.append(np.inf)
    return out

def remob(c, M=(1.3, -0.7, 0.45, 1.1)):
    """reparametrise cubic forms by t = (a s + b)/(g s + d) (avoids contacts at infinity)"""
    a, b, g, dd = M
    out = []
    for x in c:
        x = np.concatenate([np.array(x, complex), np.zeros(4-len(x))])
        y = np.zeros(4, complex)
        for k in range(4):
            y = y + x[k]*Pn.polymul(Pn.polypow([b, a], k), Pn.polypow([dd, g], 3-k))[:4] if True else y
        out.append(y)
    return out

def analyse(c):
    """c: list of 5 coefficient arrays (low first, length 4)."""
    c = remob(c)
    res = {}
    def forms(span):
        L = annihilator([PTS[a] for a in span])
        return [sum(L[i, j]*c[j] for j in range(5)) for i in range(L.shape[0])]
    for a in range(1, 7):
        r = common_roots(forms([a]))
        if r: res[('pt', a)] = r
    for ab in itertools.combinations(range(1, 7), 2):
        r = common_roots(forms(ab))
        if r: res[('line',) + ab] = r
    for abc in itertools.combinations(range(1, 7), 3):
        r = common_roots(forms(abc))
        if r: res[('plane',) + abc] = r
    for abcd in itertools.combinations(range(1, 7), 4):
        f = forms(abcd)[0]
        q = np.trim_zeros(f, 'b')
        rts = list(Pn.polyroots(q)) if len(q) > 1 else []
        rts += [np.inf]*(4 - len(q))
        res[('hyp',) + abcd] = rts
    return res

def summary(c, tol=1e-6):
    res = analyse(c)
    planes = {k[1:]: v for k, v in res.items() if k[0] == 'plane'}
    bad = [k for k in res if k[0] in ('pt', 'line')]
    contacts = []  # (t, divisor)
    for k, v in planes.items():
        for r in v: contacts.append((r, ('P',) + k))
    for k, v in res.items():
        if k[0] != 'hyp': continue
        for r in v:
            if any((r == q) if (r is np.inf or q is np.inf) else abs(r - q) < 1e-5 for q, _ in [(x, y) for x, y in contacts if y[0] == 'P' and set(y[1:]) <= set(k[1:])]):
                continue
            contacts.append((r, ('H',) + k[1:]))
    # distinctness
    ts = [x for x, _ in contacts]
    mind = np.inf
    for i in range(len(ts)):
        for j in range(i+1, len(ts)):
            a, b = ts[i], ts[j]
            if a is np.inf or b is np.inf:
                d = 0 if (a is np.inf and b is np.inf) else 1/(abs(b if a is np.inf else a)+1e-300)
            else:
                d = abs(a - b)
            mind = min(mind, d)
    return dict(planes=sorted(planes), nplanes=len(planes), bad=bad, ncontacts=len(contacts), mindist=mind,
                hyp=[y[1:] for _, y in contacts if y[0] == 'H'])
