"""NUMERICAL (referee): rational curves of class cub13 (L7 Prop 3; pair 6, triple 13, D_16 = 3) in the
frame anchors 1 -> 0, 2 -> inf, 7 -> 1, movers x3, x4, x5 (degree 2) and x6 (degree 3).
Gauge: the three contact points o123, o124, o125 (where 7 meets the cluster) are s = 0, 1, inf.
  x3 = (s-z34)(s-z35)/((s-p3)(s-p36))           x3(1) = 1
  x4 = (s-z4)(s-z34)/((s-p45)(s-p46))           x4(0) = 1
  x5 = k5 (s-z5)(s-z35)/((s-p45)(s-p56))        x5(0) = x5(1) = 1
  x6 = 1 + c s(s-1)/((s-p36)(s-p46)(s-p56))      (x6 = 1 at 0, 1, inf; zeros z6a,b,c = D_16 contacts)
  non-anchor triple collisions: x3=x4=x5 at q345, x3=x4=x6 at q346, x3=x5=x6 at q356, x4=x5=x6 at q456.
Unknowns (15): z34 z35 p3 p36 z4 p45 p46 z5 p56 k5 c q345 q346 q356 q456; equations 12; expected dim 3.
Then: (i) nondegeneracy = every split occurs exactly as prescribed (collision analysis of the 7-point
configuration at all roots of all pairwise differences), (ii) Jacobian rank 12 (local dim 3),
(iii) the evaluation map (family x P^1) -> M_{0,7} has rank 4 at a general point (covering).
"""
import numpy as np, itertools, sys
rng = np.random.default_rng(0)
NAMES = 'z34 z35 p3 p36 z4 p45 p46 z5 p56 k5 c q345 q346 q356 q456'.split()

def X(u, s):
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c = u[:11]
    x3 = (s-z34)*(s-z35)/((s-p3)*(s-p36))
    x4 = (s-z4)*(s-z34)/((s-p45)*(s-p46))
    x5 = k5*(s-z5)*(s-z35)/((s-p45)*(s-p56))
    x6 = 1 + c*s*(s-1)/((s-p36)*(s-p46)*(s-p56))
    return x3, x4, x5, x6

def F(u):
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c, q345, q346, q356, q456 = u
    R = [(1-z34)*(1-z35) - (1-p3)*(1-p36), z4*z34 - p45*p46,
         k5*z5*z35 - p45*p56, k5*(1-z5)*(1-z35) - (1-p45)*(1-p56)]
    a = X(u, q345); R += [a[0]-a[1], a[0]-a[2]]
    a = X(u, q346); R += [a[0]-a[1], a[0]-a[3]]
    a = X(u, q356); R += [a[0]-a[2], a[0]-a[3]]
    a = X(u, q456); R += [a[1]-a[2], a[1]-a[3]]
    return np.array(R)

def jac(f, u, h=1e-7):
    f0 = f(u); Jm = np.zeros((len(f0), len(u)), complex)
    for k in range(len(u)):
        d = np.zeros(len(u), complex); d[k] = h
        Jm[:, k] = (f(u+d) - f(u-d))/(2*h)
    return Jm

# polynomial (numerator/denominator) representation for collision analysis
def nd(u):
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c = u[:11]
    P = np.poly
    num = {3: P([z34, z35]), 4: P([z4, z34]), 5: k5*P([z5, z35])}
    den = {3: P([p3, p36]), 4: P([p45, p46]), 5: P([p45, p56])}
    d6 = P([p36, p46, p56]); num[6] = np.polyadd(d6, c*np.array([1, -1, 0])); den[6] = d6
    num[1] = np.array([0.0]); den[1] = np.array([1.0])          # label 1 at 0
    num[2] = np.array([1.0]); den[2] = np.array([0.0])          # label 2 at inf
    num[7] = np.array([1.0]); den[7] = np.array([1.0])          # label 7 at 1
    return num, den

def collisions(u, tol=1e-6):
    num, den = nd(u)
    L = [1, 2, 3, 4, 5, 6, 7]
    roots = []
    for i, j in itertools.combinations(L, 2):
        w = np.polysub(np.polymul(num[i], den[j]), np.polymul(num[j], den[i]))
        w = np.trim_zeros(w, 'f')
        if len(w) > 1: roots += list(np.roots(w))
    # include s = inf (leading behaviour): handled by gauge point o125 = inf
    pts = []
    for r in roots:
        if all(abs(r - p) > tol for p in pts): pts.append(r)
    events = []
    for r in pts + [np.inf]:
        if r is np.inf:
            vals = {}
            for i in L:
                n = np.trim_zeros(num[i], 'f'); d = np.trim_zeros(den[i], 'f')
                if len(n) == 0: vals[i] = 0
                elif len(d) == 0: vals[i] = np.inf
                elif len(n) > len(d): vals[i] = np.inf
                elif len(n) < len(d): vals[i] = 0
                else: vals[i] = n[0]/d[0]
        else:
            vals = {}
            for i in L:
                dv = np.polyval(den[i], r); nv = np.polyval(num[i], r)
                vals[i] = np.inf if abs(dv) < 1e-9*max(1, abs(nv)) else nv/dv
        groups = []
        for i in L:
            for g in groups:
                a, b = vals[i], vals[g[0]]
                if (a is np.inf or b is np.inf) and a is b: g.append(i); break
                if a is not np.inf and b is not np.inf and np.isfinite(a) and np.isfinite(b) and abs(a-b) < 1e-5*(1+abs(a)): g.append(i); break
            else: groups.append([i])
        big = [g for g in groups if len(g) >= 2]
        if big: events.append((r, big))
    return events

expected = {frozenset(s) for s in [(1,4),(1,5),(2,3),(1,2,3),(1,2,4),(1,2,5),(1,3,4),(1,3,5),(2,3,6),(2,4,5),(2,4,6),
            (2,5,6),(3,4,5),(3,4,6),(3,5,6),(4,5,6)]}
def canon(g):
    g = frozenset(g); c = frozenset(range(1, 8)) - g
    return g if 7 not in g else c


def diffs(u):
    z34, z35, p3, p36, z4, p45, p46, z5, p56, k5, c, q345, q346, q356, q456 = u
    special = [z34, z35, p3, p36, z4, p45, p46, z5, p56, 0, 1]
    d = [z34-p3, z34-p36, z35-p3, z35-p36, z4-p45, z4-p46, z34-p45, z34-p46, z5-p45, z5-p56, z35-p45, z35-p56,
         c, k5-1, z34-z35, p3-p36, p45-p46, p45-p56, p46-p56, p36-p46, p36-p56, z4-z34, z5-z35, z4-z5]
    qs = [q345, q346, q356, q456]
    d += [a-b for a, b in itertools.combinations(qs, 2)]
    d += [qq - sp for qq in qs for sp in special]
    return np.array(d)
ND = len(diffs(np.ones(15) + np.arange(15)*0.1j))
def G(v):
    u = v[:15]; w = v[15:]
    return np.concatenate([F(u), w*diffs(u) - 1])
