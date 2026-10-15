# rflex.py -- referee check of one-modulus FLEX witnesses (def_wit_<set>_K<K>.json):
# system = GENUINE labels at every odd prime q != q0 not dividing N (all levels; stronger than the author's
# collision-free labels at free moduli), witness lift (integers) at q0 at every level.
# Checks: (L, iota) of the file equal w*(my L, my iota) for a unit w (consistent coordinates), parity, exact budgets,
# and that the system is non-coherent with every column-sign class on F ∩ [3, q0] (mask AND = 0).
import sys, json, math
from reflib import *
from rverify import per_v_angle_ok
SETS = {'first4': [5, 13, 17, 29], 'first5': [5, 13, 17, 29, 37], 'first6': [5, 13, 17, 29, 37, 41],
        'first7': [5, 13, 17, 29, 37, 41, 53], 'first8': [5, 13, 17, 29, 37, 41, 53, 61],
        'rand4': [53, 61, 89, 173], 'bal4': [97, 101, 109, 113], 'last4': [353, 373, 389, 397]}

def odd_part_wo(x, q0):
    x = abs(x)
    while x % 2 == 0: x //= 2
    while x % q0 == 0: x //= q0
    return x

def check(name, K, which="first", Pfile=None, WI=0):
    D = json.load(open(f"../rigcompute_checks/def_wit_{name}_K{K}.json"))
    thr = json.load(open(f"../rigcompute_checks/thr_{name}_K{K}.json"))
    P = thr["P"]; r = len(P); F = sorted(thr["forced"]); pis = [gprime(p) for p in P]
    flex = [e for e in D if e["verdict"] == "FLEX"]
    if not flex: return None
    e = min(flex, key=lambda e: e["q0"]) if which == "first" else which(flex)
    q0 = e["q0"]; m = e["m"]; n = m // 4
    Lm, iota_m, mm = genuine_L(P, q0); assert mm == m
    ws = [w for w in range(m) if math.gcd(w, m) == 1 and (w * iota_m) % m == e["iota"] and
          all((w * a) % m == b for a, b in zip(Lm, e["L"]))]
    assert ws, "inconsistent coordinates"
    w = e["witnesses"][WI]; ell = w["ell"]; lift = w["lift"]
    assert all((a - b) % m == 0 for a, b in zip(lift, ell))
    nu = nu_par(P, q0); assert all(x % 2 == p for x, p in zip(lift, nu)), "parity"
    k = (e["iota"] // n) % 4
    worst = -1e9; nviol = 0; ncoll = 0; angbad = 0
    Fset = F + [q for q in thr.get("aborted", []) if q not in F]
    for v in monos(r, K):
        z = Zv(pis, v); Dv = Dvals(z)
        d = sum(a * b for a, b in zip(v, lift))
        dems = []
        for s in range(4):
            dem = 1
            for qq in Fset:
                if qq != q0:
                    t = val(Dv[s], qq)
                    if t: dem *= qq ** t
            if d % n == 0:
                assert d != 0
                sd = ((d % m) // n) * (1 if k == 1 else 3) % 4
                if sd == s:
                    dem *= q0 ** (1 + val(d // n, q0)); ncoll += 1
            c = capv(P, v, s)
            if dem * dem > c: nviol += 1
            dems.append(dem)
            worst = max(worst, math.log(dem * dem) - math.log(c))
        Nv = norm_v(P, v)
        if not angle_ok_exact(dems, Nv): angbad += 1
        pass
    # coherence on F ∩ [3, q0]: genuine masks below (+ the witness mask at q0)
    acc = (1 << 2 ** (r - 1)) - 1
    for q in F:
        if q > q0: break
        if q == q0:
            acc &= class_mask([x % n for x in ell], e["L"], n, r)
        else:
            Lq, _, mq_ = genuine_L(P, q); acc &= class_mask([x % (mq_ // 4) for x in Lq], Lq, mq_ // 4, r)
    return dict(name=name, K=K, q0=q0, violations=nviol, worst_log_slack=round(worst, 4), ncoll=ncoll,
                noncoherent=(acc == 0), angle_incompatible_v=angbad)

if __name__ == "__main__":
    for a in sys.argv[1:]:
        name, K = a.split(":"); K = int(K)
        try:
            D = json.load(open(f"../rigcompute_checks/def_wit_{name}_K{K}.json"))
            e = min([e for e in D if e["verdict"] == "FLEX"], key=lambda e: e["q0"])
            res = [check(name, K, WI=i) for i in range(len(e["witnesses"]))]
            print(name, K, e["q0"], "witnesses", len(res), "budget-ok", sum(x["violations"] == 0 for x in res),
                  "angle-compatible & noncoherent", sum(x["violations"] == 0 and x["angle_incompatible_v"] == 0 and x["noncoherent"] for x in res), flush=True)
        except FileNotFoundError as ex:
            print(a, "missing", ex)
