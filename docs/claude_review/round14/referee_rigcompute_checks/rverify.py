# rverify.py -- referee's independent verifier of a joint witness (jwmix_*.pkl), integer-label interpretation:
# the stored top lift (integers) is the label at EVERY level a >= 1 (reduced mod m q^(a-1)); collision depth of v
# at q is 1 + v_q(<v,top>/n) when n | <v,top> (infinite if <v,top> = 0).  Unit from <v,top> mod m and iota.
# Checks: coverage of the referee's own forced set, parity, exact budgets for all v (both signs), every s;
# class masks w.r.t. the referee's own genuine logs; segment coherence; per-v four-unit angle compatibility.
import sys, pickle, json, math
from reflib import *

class R(pickle.Unpickler):
    def find_class(self, mod, name):
        raise pickle.UnpicklingError(f"blocked {mod}.{name}")

def per_v_angle_ok(a):
    """a[s] = nu_s * m_s / sqrt(N) (s=0..3).  Exists x in [0,pi) with |sin(x - s pi/4)| >= a[s] for all s?
    Feasible set for s: x - s pi/4 in [asin a, pi - asin a] mod pi.  Check on a fine grid + interval logic."""
    if max(a) > 1: return False
    ivs = []
    for s in range(4):
        lo = math.asin(a[s]) + s * math.pi / 4; hi = math.pi - math.asin(a[s]) + s * math.pi / 4
        ivs.append((lo, hi))
    # candidate points: endpoints mod pi
    cands = []
    for lo, hi in ivs:
        for t in (lo, hi):
            cands.append(t % math.pi)
    def inside(x, lo, hi):
        for k in (-2, -1, 0, 1, 2):
            y = x + k * math.pi
            if lo - 1e-12 <= y <= hi + 1e-12: return True
        return False
    for x in cands:
        if all(inside(x, lo, hi) for lo, hi in ivs): return True
    return False

def verify(fn, myforced, iota_mode="mine"):
    W = R(open(fn, "rb")).load()
    P, K = W["P"], W["K"]; r = len(P)
    labs = W["labels"]
    assert sorted(labs) == sorted(myforced), ("coverage", sorted(set(labs) ^ set(myforced)))
    pis = [gprime(p) for p in P]
    V = monos(r, K)
    dem = {}  # (v,s) -> product
    masks = {}; full = (1 << 2 ** (r - 1)) - 1
    maxdepth = {}
    for q in sorted(labs):
        top = [int(x) for x in labs[q]["lifts"][-1]]
        # lifts consistency
        lifts = labs[q]["lifts"]; m = m_of(q); n = m // 4
        for i in range(1, len(lifts)):
            Mi = m * q ** (i - 1)
            assert all((a - b) % Mi == 0 for a, b in zip(lifts[i], lifts[i - 1])), ("lift", q, i)
        nu = nu_par(P, q)
        assert all(x % 2 == p for x, p in zip(top, nu)), ("parity", q)
        L, iota, mm = genuine_L(P, q); assert mm == m
        k = (iota // n) % 4  # 1 or 3
        md = 0
        for v in V:
            d = sum(a * b for a, b in zip(v, top))
            if d % n: continue
            assert d != 0, ("infinite collision", q, v)
            s = ((d % m) // n) * (1 if k == 1 else 3) % 4
            depth = 1 + val(d // n, q)
            md = max(md, depth)
            dem[(v, s)] = dem.get((v, s), 1) * q ** depth
        maxdepth[q] = md
        masks[q] = class_mask([x % n for x in top], L, n, r)
    bad = [(v, s, mm) for (v, s), mm in dem.items() if mm * mm > capv(P, v, s)]
    # per-v angle compatibility
    angbad = []
    for v in V:
        N = norm_v(P, v)
        if not angle_ok_exact([dem.get((v, s), 1) for s in range(4)], N): angbad.append((v, [dem.get((v, s), 1) for s in range(4)]))
    # segment coherence
    F = sorted(labs); acc = full; segQ = None; first_break = None
    for q in F:
        acc &= masks[q]
        if acc == 0: first_break = q; break
        segQ = q
    return dict(P=P, K=K, Qclaim=W.get("Q"), budget_violations=len(bad), bad=bad[:5], angle_incompatible=len(angbad),
                angbad=angbad[:5], coherent_upto=segQ, first_noncoherent=first_break,
                outside_all=[q for q in F if masks[q] == 0], maxdepth=maxdepth)

if __name__ == "__main__":
    fn, forcedfn = sys.argv[1], sys.argv[2]
    myF = json.load(open(forcedfn))["forced"]
    res = verify(fn, myF)
    F = sorted(myF)
    Q = res["first_noncoherent"]
    above = [q for q in F if Q and q > Q]
    out_above = [q for q in above if q in res["outside_all"]]
    print(fn.split("/")[-1], "P", res["P"], "K", res["K"], "Qclaim", res["Qclaim"],
          "| budget violations", res["budget_violations"], res["bad"],
          "| coherent on F<=", res["coherent_upto"], "first non-coherent", Q,
          "| outside every eps-class:", len(res["outside_all"]), "of", len(F),
          f"; above Q: {len(out_above)}/{len(above)}",
          "| per-v angle-incompatible v:", res["angle_incompatible"], res["angbad"][:3])
