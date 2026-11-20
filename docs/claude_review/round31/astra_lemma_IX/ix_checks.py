"""Integer checks for Part IX: PTE tests, local parity, and the Paley relaxation."""
from bisect import bisect_left
from collections import Counter, defaultdict
from itertools import combinations, product
from math import gcd, isqrt, prod


def prime_powers(k):
    k = abs(k)
    p = 2
    while p * p <= k:
        power = 1
        while k % p == 0:
            power *= p
            k //= p
        if power > 1:
            yield p, power
        p += 1
    if k > 1:
        yield k, k


def local_signatures(marks, factors):
    signatures = [0] * len(marks)
    for bit, modulus in enumerate(factors):
        residues = sorted({y % modulus for y in marks})
        assert len(residues) <= 2
        for i, y in enumerate(marks):
            signatures[i] |= (residues.index(y % modulus) << bit)
    return signatures


def pte(label, nodes, r, value):
    # node b/a, homogeneous moment vector (a,b), and coefficient a^2.
    moduli = [a * value - b for a, b in nodes]
    factors = {}
    for k in moduli:
        for p, power in prime_powers(k):
            factors[p] = max(factors.get(p, 1), power)
    total = prod(moduli)
    marks = []
    for flags in product((0, 1), repeat=len(nodes)):
        if sum(flag * a for flag, (a, b) in zip(flags, nodes)):
            continue
        if r == 2 and sum(flag * b for flag, (a, b) in zip(flags, nodes)):
            continue
        marks.append(sum(flag * a**2 * (total // k)
                         for flag, (a, b), k in zip(flags, nodes, moduli)))
    check_integer_family(label, marks, factors)


def check_integer_family(label, marks, factors):
    q = prod(factors.values())
    assert len(marks) == len(set(marks))
    sig = local_signatures(marks, factors.values())
    L = max(marks) - min(marks)
    edges = defaultdict(list)
    planes = set()
    for i, j in combinations(range(len(marks)), 2):
        h = sig[i] ^ sig[j]
        edges[h].append((i, j))
    for h, pairs in edges.items():
        H = prod(v for bit, v in enumerate(factors.values()) if not (h >> bit) & 1)
        K = q // H
        assert all((marks[i]-marks[j]) % H == 0 for i, j in pairs)
        assert len({(marks[i]+marks[j]) % K for i, j in pairs}) == 1
        assert len({(((marks[i]-marks[j])//H)**2) % K for i, j in pairs}) == 1
        assert len(pairs) <= (1 + 2*L // K) * (1 + L**2 // (H**2 * K))
        for (i,j), (k,l) in combinations(pairs, 2):
            plane = frozenset((sig[i],sig[j],sig[k],sig[l]))
            if len(plane) == 4:
                planes.add(plane)
    assert not planes or (2*L)**2 >= q
    code = set(sig)
    cubes = set()
    for plane in planes:
        x = next(iter(plane))
        for v in code - plane:
            cube = plane | frozenset(y ^ x ^ v for y in plane)
            if cube <= code:
                cubes.add(cube)
    assert not cubes or (4*L)**3 >= q*q
    H = {(i,j): prod(v for bit, v in enumerate(factors.values())
                     if not ((sig[i] ^ sig[j]) >> bit) & 1)
         for i,j in combinations(range(len(marks)),2)}
    for i,j,k in combinations(range(len(marks)),3):
        pairs = ((i,j),(j,k),(i,k))
        h = [H[pair] for pair in pairs]
        common = gcd(*h)
        blocks = [v // common for v in h]
        assert all(gcd(a,b) == 1 for a,b in combinations(blocks,2))
        t = [(marks[v]-marks[u])//v_h for (u,v),v_h in zip(pairs,h)]
        assert t[0]*blocks[0]+t[1]*blocks[1] == t[2]*blocks[2]
    print(label, "marks", len(marks), "planes", len(planes), "3-flats", len(cubes))


def maximal_modulus_family(label, marks):
    # For each binary partition take the gcd of its two sets of internal
    # differences. Their lcm gives the maximum modulus with two residues.
    factors = {}
    for flags in product((0,1), repeat=len(marks)-1):
        bits = (0,)+flags
        g = 0
        for i,j in combinations(range(len(marks)),2):
            if bits[i] == bits[j]:
                g = gcd(g,marks[i]-marks[j])
        for p,pf in prime_powers(g):
            factors[p] = max(factors.get(p,1),pf)
    check_integer_family(label,marks,factors)


def explicit_families():
    a,b = 1009,431
    maximal_modulus_family("parallelogram", [0,a,b,a+b])
    f = a*(a+b)
    g = b*(a+2*b)
    maximal_modulus_family("5 projective marks", [0,f,g,f+g,(a+b)*(a+2*b)])
    maximal_modulus_family("6 projective marks", [0,(a-b)*(a-2*b),a*(a-2*b),
                           -2*b*(a-b),-a*b,a*a-4*a*b+2*b*b])
    nodes = (1009,1013,1019)
    marks = [sum((prod(nodes)//v)*bit for v,bit in zip(nodes,bits))
             for bits in product((0,1),repeat=3)]
    maximal_modulus_family("cube",marks)


def primes_upto(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p)//p + 1)
    return [p for p in range(2,limit+1) if sieve[p]]


def fake(m, P):
    primes = [p for p in primes_upto(2*P) if p >= P and p % 4 == 1]
    T = (3*P//2)**2
    used = set()
    blocks = []
    for p in primes:
        if p in used:
            continue
        for j in range(bisect_left(primes, (T+p-1)//p), len(primes)):
            l = primes[j]
            if (p*l-T)*m*8 > T:
                break
            if l != p and l not in used:
                blocks.append(p*l)
                used.update((p,l))
                break
        if len(blocks) == m:
            break
    assert len(blocks) == m and len(used) == 2*m
    q = prod(blocks)
    residues = {i*i % m for i in range(1,m)}
    rows = [sum(1 << v for v in range(m) if (v-i) % m in residues) for i in range(m)]
    k = P//(4*m)
    labels = [k*i + int(i == m-1) for i in range(m)]
    diffs = {}
    for i,j in combinations(range(m),2):
        xor = rows[i] ^ rows[j]
        assert xor.bit_count() == (m+1)//2
        H = prod(blocks[v] for v in range(m) if not (xor >> v) & 1)
        d = H * (labels[j]-labels[i])
        assert d*d <= q and gcd(d,q) == H
        diffs[i,j] = d
    pair_xors = Counter(rows[i]^rows[j] for i,j in diffs)
    assert max(pair_xors.values()) == 1
    fails = sum(diffs[i,j]+diffs[j,l] != diffs[i,l] for i,j,l in combinations(range(m),3))
    assert fails
    print("Paley", m, "prime factors",len(used),"different pair sums",len(pair_xors),
          "failed signed triangles",fails)


def local_checks():
    # The two minima in the lattice obstruction are attained locally.
    for p, y, size, minimum in [(2,[0,4,2,6,1,3,5],4,5), (3,[0,3,6,1,4,2,5],3,2)]:
        cof = [(y[i]-y[j])//(p if j<size else 1) for i,j in combinations(range(7),2)]
        assert sum(c % p == 0 for c in cof) == minimum
    # The three nonconstant coefficients of each two-bit Boolean function.
    l1 = []
    for mask in range(16):
        coeff = [sum((mask >> x & 1)*(-1)**((u & x).bit_count() % 2)
                     for x in range(4)) for u in (1,2,3)]
        l1.append(sum(map(abs,coeff)))
    assert max(l1) == 3
    print("two-bit tables", sorted(Counter(l1).items()), "local cofactor minima", (5,2))


def gauss_mul(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def gauss_conj(z):
    return z[0], -z[1]


def gauss_div(z, w):
    x,y = gauss_mul(z,gauss_conj(w))
    norm = w[0]**2+w[1]**2
    assert x % norm == y % norm == 0
    return x//norm,y//norm


def gauss_product(factors):
    out = (1,0)
    for factor in factors:
        out = gauss_mul(out,factor)
    return out


def gaussian_triangles():
    primes = ((2,1),(3,2),(4,1))
    cap = 2
    rows = list(product(range(cap+1),repeat=len(primes)))
    def common_factor(indices, slots=range(len(primes))):
        factors = []
        for v in slots:
            exps = [rows[i][v] for i in indices]
            factors.extend([primes[v]]*min(exps))
            factors.extend([gauss_conj(primes[v])]*(cap-max(exps)))
        return gauss_product(factors)
    marks = [common_factor([i]) for i in range(len(rows))]
    assert len(set(marks)) == len(rows)
    num = forced = 0
    for i,j,k in combinations(range(len(rows)),3):
        triple = (i,j,k)
        gs = common_factor(triple)
        signed = []
        hidden = [v for v in range(len(primes)) if rows[i][v]==rows[j][v]!=rows[k][v]]
        known = [v for v in range(len(primes)) if v not in hidden]
        for u,v in ((i,j),(j,k),(i,k)):
            H = common_factor((u,v))
            delta = tuple(marks[v][a]-marks[u][a] for a in (0,1))
            real,imag = gauss_div(delta,H)
            assert real == 0 and imag and imag%2 == 0
            coef = imag//2
            r = gauss_div(H,gs)
            signed.append((coef*r[0],coef*r[1]))
            if (u,v)==(i,j):
                outside = gauss_div(common_factor((i,j),known), common_factor(triple,known))
                multiplier = (coef*outside[0],coef*outside[1])
            else:
                assert common_factor((u,v),hidden) == common_factor(triple,hidden)
        required = tuple(signed[2][a]-signed[1][a] for a in (0,1))
        assert required == signed[0]
        g = gauss_div(required,multiplier)
        recovered_norm = g[0]**2+g[1]**2
        assert recovered_norm == prod((primes[v][0]**2+primes[v][1]**2)
                                    **abs(rows[i][v]-rows[k][v]) for v in hidden)
        num += 1
        forced += bool(hidden)
    print("Gaussian cap 2 triangles",num,"with hidden slots for designated pair",forced)


if __name__ == "__main__":
    local_checks()
    explicit_families()
    for k in (4,5,8):
        a = [1]*((k+1)//2) + [-1]*(k//2)
        pte(f"{k} nodes", [(v,v*j) for j,v in enumerate(a)], 1, 1009)
    vectors = [(-3,-3),(-3,-2),(0,-1),(1,2),(1,3),(2,0),(2,1)]
    pte("7 nodes, 2 moments", [(a+7*b,b) for a,b in vectors], 2, 1009)
    fake(43, 40000)
    gaussian_triangles()
