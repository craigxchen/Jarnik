"""Part V cluster over F_2 (norm A^2+AB+B^2 of F_4(t)/F_2(t)), built with own arithmetic.
Usage: python3 cluster.py r [L]   (m = 2^r, d = r+1; default L = (d-2) 2^(d-1))
Writes a JSON summary to cluster_m{m}.json
"""
import sys, json, itertools, random, time
from f24 import *

def parity(x):
    return bin(x).count('1') & 1

def setup(r, L=None):
    d = r + 1
    Q = 1 << d
    m = Q // 2
    if L is None:
        L = (d - 2) * (1 << (d - 1))
    forms = [a for a in range(Q) if a & 1]          # l(v) = parity(a & v), l(1) = 1
    H = {a: [v for v in range(Q) if parity(a & v) == 0] for a in forms}
    Yv = {v: ((1 << L) ^ v, 1) for v in range(Q)}    # Y + v,  Y = t^L + theta
    return d, Q, m, L, forms, H, Yv

def monic_polys4(k):
    for cs in itertools.product(range(4), repeat=k):
        yield poly4_from_coeffs(list(cs) + [1])

def small_primes(d):
    """monic irreducible f in F_4[t], not self-conjugate, with 2 deg f < d"""
    irr = []
    for k in range(1, (d - 1) // 2 + 1):
        for f in monic_polys4(k):
            ok = True
            for g in irr:
                if deg4(g) <= k // 2:
                    _, rr = divmod4(f, g)
                    if iszero4(rr):
                        ok = False
                        break
            # also self-conjugate irreducibles of lower degree could divide; test all monic of deg<=k/2
            if ok:
                for j in range(1, k // 2 + 1):
                    for g in monic_polys4(j):
                        _, rr = divmod4(f, g)
                        if iszero4(rr):
                            ok = False
                            break
                    if not ok:
                        break
            if ok:
                irr.append(f)
    return [f for f in irr if f[1] != 0 or True and conj4(f) != f]

def val4(P, f):
    k = 0
    while True:
        q, rr = divmod4(P, f)
        if not iszero4(rr):
            return k
        P = q
        k += 1

def main():
    r = int(sys.argv[1])
    Larg = int(sys.argv[2]) if len(sys.argv) > 2 else None
    t0 = time.time()
    d, Q, m, L, forms, H, Yv = setup(r, Larg)
    D = (m - 1) * L
    W = (Q - 2) * L
    out = dict(m=m, d=d, Q=Q, L=L, D=D, W=W)
    sp = small_primes(d)
    # orientation of small primes: pick one of each conjugate pair
    orient = []
    seen = set()
    for f in sp:
        if f in seen:
            continue
        fb = monic4(conj4(f))
        seen.add(f); seen.add(fb)
        orient.append(f)
    # valuations of Y+v at all small primes
    val = {}
    for f in sp:
        for v in range(Q):
            k = val4(Yv[v], f)
            if k:
                val[(f, v)] = k
    # self-check: no self-conjugate small irreducible divides any Y+v
    bigdeg = {}
    for v in range(Q):
        s = sum(deg4(f) * val.get((f, v), 0) for f in sp)
        bigdeg[v] = L - s
    out['n_small_primes'] = len(sp)
    out['small_valuations_nonzero'] = len(val)
    # points
    z = {}
    for a in forms:
        P = (1, 0)
        for v in H[a]:
            if v:
                P = mul4(P, Yv[v])
        z[a] = P
    Npoly = (1, 0)
    for v in range(2, Q):
        Npoly = mul4(Npoly, Yv[v])
    assert Npoly[1] == 0, "N not in F_2[t]"
    N = Npoly[0]
    assert deg2(N) == W
    normok = all(mul4(z[a], conj4(z[a])) == (N, 0) for a in forms)
    out['norms_ok'] = normok
    out['distinct'] = len(set(z.values())) == m
    out['deg_points'] = sorted(set(deg4(z[a]) for a in forms))
    pairs = list(itertools.combinations(forms, 2))
    ddiff = {}
    for a, b in pairs:
        ddiff[(a, b)] = deg4(add4(z[a], z[b]))
    out['max_diff_deg'] = max(ddiff.values())
    out['D_over_2'] = D / 2
    # valuations of points at small primes
    vz = {(f, a): sum(val.get((f, v), 0) for v in H[a] if v) for f in sp for a in forms}
    vN = {f: sum(val.get((f, v), 0) for v in range(2, Q)) for f in sp}
    # gcd degrees via valuations
    def gdeg(a, b):
        J = [v for v in H[a] if v and parity(b & v) == 0]
        return sum(bigdeg[v] for v in J) + sum(deg4(f) * min(vz[(f, a)], vz[(f, b)]) for f in sp)
    gd = {(a, b): gdeg(a, b) for a, b in pairs}
    # spot check with direct gcds
    random.seed(0)
    chk = random.sample(pairs, min(6, len(pairs)))
    direct_ok = True
    for a, b in chk:
        g = gcd4(z[a], z[b])
        if deg4(g) != gd[(a, b)]:
            direct_ok = False
    out['gcd_valuation_vs_direct_ok(6 pairs)'] = direct_ok
    # cut slack Q (orientation-free, via gcds) and via layers
    Qslack_g = m * m * W / 4 - sum(W - 2 * gd[p] for p in pairs)
    Qslack_l = 0
    layers_nontriv = []
    for f in orient:
        e = 2 * deg4(f)  # degree of the F_2 prime
        for l in range(1, vN[f] + 1):
            s = sum(1 for a in forms if vz[(f, a)] >= l)
            Qslack_l += e * (s - m / 2) ** 2
            if s != m // 2:
                layers_nontriv.append((e, l, s))
    out['Qslack_from_gcds'] = Qslack_g
    out['Qslack_from_layers'] = Qslack_l
    out['nontrivial_layers(e,depth,s)'] = layers_nontriv[:40]
    # (IV.2) with delta = 0 (D even) or -1/2 (D odd), E = 0
    K = D // 2
    delta = K - D / 2
    lhs = 0.75 * (m - 24) * W + Qslack_l
    rhs = m * (m + 3) * delta
    out['IV2_delta'] = delta
    out['IV2_LHS'] = lhs
    out['IV2_RHS'] = rhs
    out['IV2_violated'] = lhs > rhs
    out['BMY_slack_upper_bound_3eorb_minus_KB2'] = -(2 / 9) * (lhs - rhs)
    # transferred cut identity: sum deg c_ij + Q/2 = sum deg(z_i-z_j) - m(m-2)D/4
    sum_c = sum(ddiff[p] - gd[p] for p in pairs)
    out['sum_deg_c'] = sum_c
    out['max_deg_c'] = max(ddiff[p] - gd[p] for p in pairs)
    out['identity_lhs'] = sum_c + Qslack_l / 2
    out['identity_rhs'] = sum(ddiff.values()) - m * (m - 2) * D / 4
    out['chord_budget_mD/4'] = m * D / 4
    out['time'] = time.time() - t0
    # save data for later steps
    import pickle
    with open(f'cluster_m{m}_L{L}.pkl', 'wb') as fh:
        pickle.dump(dict(d=d, Q=Q, m=m, L=L, D=D, W=W, forms=forms, H=H, sp=sp, orient=orient,
                         val=val, bigdeg=bigdeg, vz=vz, vN=vN, z=z, N=N, ddiff=ddiff, gd=gd), fh)
    print(json.dumps(out, indent=1, default=str))
    with open(f'cluster_m{m}_L{L}.json', 'w') as fh:
        json.dump(out, fh, indent=1, default=str)

if __name__ == '__main__':
    main()
