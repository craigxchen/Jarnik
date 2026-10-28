"""Exact actual clusters (normalised) for the Arakelov checks.

A cluster is returned as a dict with keys
  name, n0 (original norm), pts (normalised Gaussian integers, angle-sorted), N (normalised norm),
  M, C (exact-arc constant L/sqrt(R) of the ORIGINAL window, float), g (removed gcd).
"""
import math
import random
from sympy import factorint
from gauss import sub, mul, conj, norm, gcd_g, exact_div, lattice_points

def ang(z):
    return math.atan2(z[1], z[0])

def windows(n0, k):
    """All windows of k angle-consecutive lattice points of norm n0 (cyclic), with arc constant."""
    pts = sorted(lattice_points(n0), key=ang)
    P = len(pts)
    R = math.sqrt(n0)
    out = []
    if P < k:
        return out
    for s in range(P):
        w = [pts[(s + t) % P] for t in range(k)]
        a0 = ang(w[0]); a1 = ang(w[-1])
        d = (a1 - a0) % (2 * math.pi)
        L = R * d
        out.append((L / math.sqrt(R), w))
    return out

def best_window(n0, k):
    ws = windows(n0, k)
    return min(ws, key=lambda t: t[0])

def normalise(pts):
    g = pts[0]
    for z in pts[1:]:
        g = gcd_g(g, z)
    q = [exact_div(z, g) for z in pts]
    return q, g

def make(name, n0, pts, C):
    q, g = normalise(pts)
    N = norm(q[0])
    assert all(norm(z) == N for z in q)
    return dict(name=name, n0=n0, pts=q, N=N, M=len(q), C=C, g=g)

def standard_clusters(seed=1):
    out = []
    # the referee's primitive 10-point cluster and its 9-point subclusters
    C10, w10 = best_window(1176852625, 10)
    out.append(make('ten_1176852625', 1176852625, w10, C10))
    out.append(make('nine_drop_last', 1176852625, w10[:-1], None))
    out.append(make('nine_drop_first', 1176852625, w10[1:], None))
    # explicit record clusters (growth/families/check_families.py, check_sporadic)
    from decimal import Decimal
    sporadic = [
        ('eight_axis_5743586', 32988780417125, [(5743582, 6799), (5743583, 5894), (5743585, 3430), (5743586, 527),
                          (5743586, -527), (5743585, -3430), (5743583, -5894), (5743582, -6799)]),
        ('eight_axis_590749322', 348984761607530165, [(590749322-d, sg*y) for d, y in ((0, 12809), (1, 36682), (3, 60898), (4, 69929)) for sg in (1, -1)]),
        ('seven_sporadic', 235370525, [(14650, 4555), (14629, 4622), (14555, 4850), (14554, 4853), (14491, 5038),
                     (14453, 5146), (14443, 5174)]),
        ('six_axis', 469977431545, [(685547, 1656), (685548, 1171), (685549, 12), (685549, -12), (685548, -1171),
                        (685547, -1656)]),
        ('five_sporadic', 683617125425, [(706952, 428761), (706544, 429433), (706513, 429484), (706087, 430184),
                        (706009, 430312)]),
    ]
    for name, n0, pts in sporadic:
        assert all(x*x+y*y == n0 for x, y in pts)
        pts = sorted(pts, key=ang)
        R = math.sqrt(n0)
        d = (ang(pts[-1]) - ang(pts[0])) % (2*math.pi)
        out.append(make(name, n0, pts, R*d/math.sqrt(R)))
    # best 6- and 7-point windows on a few highly composite norms
    for n0 in (5**2*13*17*29*37*41, 5*13*17*29*37*41*53, 5**3*13**2*17*29):
        for k in (5, 6, 7):
            Ck, wk = best_window(n0, k)
            out.append(make('w%d_%d' % (k, n0), n0, wk, Ck))
    # random subsets (not arcs) to exercise the algebra on general configurations
    rng = random.Random(seed)
    for n0 in (5**4*13*17, 5*13**2*29*37, 65*85*29):
        pts = lattice_points(n0)
        for M in (4, 7, 10):
            if len(pts) >= M:
                sub_ = rng.sample(pts, M)
                out.append(make('rand%d_%d' % (M, n0), n0, sorted(sub_, key=ang), None))
    return out

if __name__ == '__main__':
    for c in standard_clusters():
        print(c['name'], c['M'], c['N'], factorint(c['N']), c['C'], c['g'])
