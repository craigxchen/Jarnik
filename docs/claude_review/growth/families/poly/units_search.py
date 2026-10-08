# Broader fixed-shift quadratic-unit template search (E=4: four binary shifts 0,1,2,3, or shifts
# 0,1,2 with widths 1,2,1), over Gaussian recurrence sequences H_r=(A u_{r+1}+B u_r)+i(C u_{r+1}+D u_r),
# u_{r+1}=t u_r + s u_{r-1}, u_0=0,u_1=1 (unit lambda root of X^2-tX-s, s=+-1).
# Prefactors: rows of one parity class have leading phase offsets (2|a|-E)theta; we search Gaussian
# c_e (equal norm <= NB) with arg c_e + (2|a|-E)theta all equal mod pi/2.
import itertools, math
from gpoly import gmul
from golden import subspans
NB = 65
gauss = [(x, y) for x in range(-9, 10) for y in range(-9, 10) if 0 < x*x+y*y <= NB]
def seq(t, s, A, B, C, D, rmax):
    u = [0, 1]
    while len(u) < rmax + 3: u.append(t*u[-1] + s*u[-2])
    return [(A*u[r+1] + B*u[r], C*u[r+1] + D*u[r]) for r in range(rmax)]
def conj(z): return (z[0], -z[1])
def gpow(z, e):
    p = (1, 0)
    for _ in range(e): p = gmul(p, z)
    return p
def angle_mod(x, m=math.pi/2):
    x = x % m
    return x if x < m/2 else x - m
def find_prefactors(theta, offsets):
    # offsets: sorted list of distinct integers o (=2|a|-E); need c_o of equal norm with
    # arg c_o + o*theta equal mod pi/2.  Return dict o->c with minimal common norm, or None.
    best = None
    bynorm = {}
    for c in gauss: bynorm.setdefault(c[0]**2 + c[1]**2, []).append(c)
    for nrm in sorted(bynorm):
        cs = bynorm[nrm]
        o0 = offsets[0]
        for c0 in cs:
            target = math.atan2(c0[1], c0[0]) + o0*theta
            sel = {o0: c0}
            ok = True
            for o in offsets[1:]:
                found = None
                for c in cs:
                    if abs(angle_mod(math.atan2(c[1], c[0]) + o*theta - target)) < 1e-9:
                        found = c; break
                if found is None: ok = False; break
                sel[o] = found
            if ok: return sel
    return None
results = []
templates = [([0,1,2,3],[1,1,1,1]), ([0,1,2],[1,2,1])]
for (t, s) in [(1,1),(2,1),(3,1),(3,-1),(4,-1)]:
    lam = (t + math.sqrt(t*t + 4*s))/2
    for A, B, C, D in itertools.product(range(-2, 3), repeat=4):
        if A*D - B*C == 0: continue
        if (A, B) == (0, 0) or (C, D) == (0, 0): continue
        theta = math.atan2(C*lam + D, A*lam + B)
        H = seq(t, s, A, B, C, D, 80)
        if any(h == (0, 0) for h in H[20:]): continue
        for shifts, widths in templates:
            E = sum(widths)
            box = list(itertools.product(*[range(w+1) for w in widths]))
            for par in (0, 1):
                rows = [a for a in box if sum(a) % 2 == par]
                offs = sorted({2*sum(a) - E for a in rows})
                pre = find_prefactors(theta, offs)
                if pre is None: continue
                best = {}
                for n in range(40, 64):
                    pts = []
                    for a in rows:
                        z = pre[2*sum(a) - E]
                        for sft, w, x in zip(shifts, widths, a):
                            b = H[n + sft]
                            z = gmul(z, gmul(gpow(b, x), gpow(conj(b), w - x)))
                        pts.append(z)
                    try:
                        bb, N, g = subspans(pts)
                    except AssertionError:
                        bb = None
                    if bb is None: continue
                    for k, v in bb.items():
                        if k not in best or v < best[k]: best[k] = v
                if best:
                    results.append((best.get(len(rows), 9e9), t, s, (A, B, C, D), shifts, widths, par, best))
results.sort(key=lambda r: r[0])
seen = set()
for r in results[:400]:
    key = (r[1], r[2], r[4] == [0,1,2,3], round(r[0], 6))
    if key in seen: continue
    seen.add(key)
    print("C_full=%.6f t=%d s=%d ABCD=%s shifts=%s widths=%s par=%d  per-k: %s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], {k: round(v, 4) for k, v in sorted(r[7].items()) if k >= 5}))
    if len(seen) > 14: break
