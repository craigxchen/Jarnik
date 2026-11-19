# Test VI.2 on Part V clusters and VI.5 squarefree clusters, m = 4, 8, 16.
import sys, math, itertools
from f4lib import *
from clusters import *

def check(name, C, full_identity=True):
    m = C['m']; Y = C['Y']; V = C['V']
    marks = hyperplane_marks(V, C['E0basis'], Y)
    assert len(marks) == m
    zs = [z for _, _, z in marks]
    D = d4(zs[0])
    assert all(d4(z) == D for z in zs)
    assert all(lc4(z) == 1 for z in zs)
    assert len(set(zs)) == m
    N = mul4(zs[0], conj4(zs[0]))
    assert is_real(N)
    normsok = all(mul4(z, conj4(z)) == N for z in zs)
    maxdiff = max(d4(add4(zs[i], zs[j])) for i in range(m) for j in range(i + 1, m))
    N2 = N[0]
    dN = deg(N2)
    Nprime = 0
    for k in range(1, dN + 1, 2):  # derivative in char 2: odd-degree terms
        if (N2 >> k) & 1:
            Nprime |= 1 << (k - 1)
    sqf = gcd2(N2, Nprime) == 1
    G = zs[0]
    for z in zs[1:]:
        G = gcd4(G, z)
    print(f"== {name}: m={m} L={C['L']} D={D} deg N={dN} norms_equal={normsok} maxdiff={maxdiff} "
          f"(D/2={D/2}) critical={maxdiff <= D//2} squarefree={sqf} deg gcd(points)={d4(G)}"
          + (f" deg R={C['degR']}" if 'degR' in C else ""))
    # labels and blocks
    labels = sorted(marks[0][1].keys())
    assert len(labels) == m - 1
    # normalized: is there an all-plus row?
    allplus = [i for i, (_, s, _) in enumerate(marks) if all(s[v] == 1 for v in labels)]
    print("   all-plus rows:", allplus)
    f = {v: add4(Y, (v, 0)) for v in labels}
    fb = {v: conj4(f[v]) for v in labels}
    # block structure exactly reproduces the marks?
    for (_, s, z) in marks:
        assert prod4([f[v] if s[v] == 1 else fb[v] for v in labels]) == z
    # Hadamard orthogonality and balance
    chi = {v: [marks[i][1][v] for i in range(m)] for v in labels}
    for v in labels:
        assert sum(chi[v]) == 0
        for w in labels:
            if v < w:
                assert sum(a * b for a, b in zip(chi[v], chi[w])) == 0
    # disjoint supports? f_v coprime to fbar_v?
    cop = all(d4(gcd4(f[v], fb[v])) == 0 for v in labels)
    disj = all(d4(gcd4(mul4(f[v], fb[v]), mul4(f[w], fb[w]))) == 0 for v in labels for w in labels if v < w)
    print(f"   f_v coprime to conj: {cop}; disjoint supports of distinct labels: {disj}")
    # exact products
    if full_identity:
        ok1 = True
        for v in labels[: min(len(labels), 6)]:
            lhs = prod4([zs[i] for i in range(m) if chi[v][i] == 1] + [pow4(fb[v], m // 2)])
            rhs = prod4([zs[i] for i in range(m) if chi[v][i] == -1] + [pow4(f[v], m // 2)])
            ok1 &= (lhs == rhs)
        ok2 = True
        for v, w in list(itertools.combinations(labels, 2))[:6]:
            pos = [i for i in range(m) if chi[v][i] == 1 and chi[w][i] == -1]
            neg = [i for i in range(m) if chi[v][i] == -1 and chi[w][i] == 1]
            lhs = prod4([zs[i] for i in pos] + [pow4(mul4(fb[v], f[w]), m // 4)])
            rhs = prod4([zs[i] for i in neg] + [pow4(mul4(f[v], fb[w]), m // 4)])
            ok2 &= (lhs == rhs)
        print(f"   prod z^chi_v = r_v^(m/2): {ok1};  prod z^((chi_v-chi_w)/2) = (r_v/r_w)^(m/4): {ok2}")
    # depths and inequalities
    K = (D + 1) // 2
    dv = {v: d4(f[v]) for v in labels}
    assert sum(dv.values()) == D
    Bv = {v: add4(f[v], fb[v]) for v in labels}  # = (B,0)
    depth = {v: depth_ratio(f[v], fb[v]) for v in labels}
    assert all(depth[v] == dv[v] - d4(Bv[v]) for v in labels)
    sv = {v: dv[v] - D / m for v in labels}
    ineq1 = all(depth[v] >= D / m for v in labels)
    ineq1b = all(d4(Bv[v]) <= sv[v] + 1e-9 for v in labels)
    worst2 = None; ok2 = True
    for v, w in itertools.combinations(labels, 2):
        A_v, B_v = f[v]; A_w, B_w = f[w]
        num = clmul(A_v, B_w) ^ clmul(A_w, B_v)
        # cross-check numerator against f_v fbar_w + fbar_v f_w
        assert add4(mul4(f[v], fb[w]), mul4(fb[v], f[w])) == (num, 0)
        dep = dv[v] + dv[w] - deg(num)
        slack = sv[v] + sv[w] - deg(num)
        if worst2 is None or slack < worst2[0]:
            worst2 = (slack, v, w, deg(num), dep)
        ok2 &= (dep >= 2 * D / m - 1e-9)
    print(f"   depth(r_v) >= D/m all: {ineq1}; deg B_v <= s_v all: {ineq1b}; "
          f"min depth={min(depth.values())} vs D/m={D/m:.2f}; deg B_v set={sorted(set(d4(Bv[v]) for v in labels))}; "
          f"s_v range=({min(sv.values()):.2f},{max(sv.values()):.2f}) sum s_v={sum(sv.values()):.2f} = D/m={D/m:.2f}")
    print(f"   depth(r_v/r_w) >= 2D/m all: {ok2}; tightest pair slack s_v+s_w-deg(A_vB_w+A_wB_v) = {worst2[0]:.2f}")
    # actual depth of the products vs K
    if full_identity:
        v = labels[0]
        rv_pow_num = pow4(f[v], m // 2); rv_pow_den = pow4(fb[v], m // 2)
        dp = depth_ratio(rv_pow_num, rv_pow_den)
        print(f"   depth(r_v^(m/2)) = {dp} = (m/2)*depth(r_v) = {(m//2)*depth[v]}; K = {K}")
    s = 2 * D / (m * (m - 1))
    bound = m * (m - 1) * (math.log2(m - 1) - 3) / 6
    print(f"   s = 2D/(m(m-1)) = {s:.2f}; (m-1)/2 = {(m-1)/2} <= 2^(3s+2) = 2^{3*s+2:.1f}: {(m-1)/2 <= 2**(3*s+2)}; "
          f"(VI.2) bound {bound:.1f} <= D = {D}: {bound <= D}; ratio D/(m^2 log2 m) = {D/(m*m*math.log2(m)):.3f}")
    # count of labels with s_v <= s
    nsmall = sum(1 for v in labels if sv[v] <= s + 1e-9)
    print(f"   labels with s_v <= s: {nsmall} >= (m-1)/2 = {(m-1)/2}")

if __name__ == '__main__':
    for r in (2, 3, 4):
        check(f"Part V r={r}", partV(r), full_identity=True)
    for d in (3, 4, 5):
        C = sqfree_VI5(d)
        check(f"VI.5 d={d}", C, full_identity=True)
