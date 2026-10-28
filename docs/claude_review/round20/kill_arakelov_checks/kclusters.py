"""Independent reconstruction of the proposal's configurations (same inputs, own code)."""
import math, random
from kgauss import *

def normalise(pts):
    g = pts[0]
    for z in pts[1:]:
        g = ggcd(g, z)
    return [gdiv_exact(z, g) for z in pts], g

def arcC(n0, pts):
    R = math.sqrt(n0)
    d = (ang(pts[-1]) - ang(pts[0])) % (2*math.pi)
    return R*d/math.sqrt(R)

def best_window(n0, k):
    pts = sorted(lattice_points(n0), key=ang)
    P = len(pts)
    best = None
    for s in range(P):
        w = [pts[(s+t) % P] for t in range(k)]
        c = arcC(n0, w)
        if best is None or c < best[0]:
            best = (c, w)
    return best

def configs(seed=1):
    out = []
    C10, w10 = best_window(1176852625, 10)
    out.append(('ten_1176852625', 1176852625, w10, C10, True))
    out.append(('nine_drop_last', 1176852625, w10[:-1], None, False))
    out.append(('nine_drop_first', 1176852625, w10[1:], None, False))
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
        out.append((name, n0, pts, arcC(n0, pts), True))
    for n0 in (5**2*13*17*29*37*41, 5*13*17*29*37*41*53, 5**3*13**2*17*29):
        for k in (5, 6, 7):
            c, w = best_window(n0, k)
            out.append(('w%d_%d' % (k, n0), n0, w, c, True))
    rng = random.Random(seed)
    for n0 in (5**4*13*17, 5*13**2*29*37, 65*85*29):
        pts = lattice_points(n0)
        for M in (4, 7, 10):
            if len(pts) >= M:
                sub_ = rng.sample(pts, M)
                out.append(('rand%d_%d' % (M, n0), n0, sorted(sub_, key=ang), None, False))
    return out
