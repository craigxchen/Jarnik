# Judge's check of the mod-four support bound (VIII.3, "An uncertainty estimate for congruences").
# (1) exact eigen-identity C y = g y for y(b) = (C-1)b/2 + b*omega, omega = (1+g)/2, g^2 = -q;
# (2) the reduction mod a prime above 2 of the lattice {y(b)}: dimension and exact minimum weight;
# (3) the transfer: supp(ybar) subset supp(v) u supp(Kv mod 4), and one of ybar_+-, ybar_- nonzero,
#     on random primitive balanced v at larger q; minimum support of Kv mod 4 against sqrt(q)+1-h.
import sys, math, random, itertools

def chi(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1

def conf(q):
    # indices 0..q-1 finite, q = infinity
    n = q+1
    C = [[0]*n for _ in range(n)]
    for i in range(q):
        for j in range(q):
            C[i][j] = chi(j-i, q)
        C[i][q] = 1; C[q][i] = -1
    return C

def matvec(C, b): return [sum(c*x for c, x in zip(row, b)) for row in C]

# elements of Z[omega] as pairs (x, y) = x + y*omega, omega^2 = omega - c, c = (1+q)/4
def check_identity(q, trials=20):
    C = conf(q); n = q+1; c = (1+q)//4
    for _ in range(trials):
        b = [random.randint(-3, 3) for _ in range(n)]
        if sum(b) % 2: b[0] += 1        # C b = b mod 2 iff sum b even (C = J - I mod 2)
        Cb = matvec(C, b)
        assert all((x-y) % 2 == 0 for x, y in zip(Cb, b))
        y = [((x - z)//2, z) for x, z in zip(Cb, b)]     # (C-1)b/2 + b*omega
        # C y
        Cy0 = matvec(C, [t[0] for t in y]); Cy1 = matvec(C, [t[1] for t in y])
        # g y with g = 2 omega - 1: (2w-1)(x + y w) = -x + (2x - y) w + 2y w^2 = -x - 2cy + (2x + y) w
        for k in range(n):
            x, z = y[k]
            gy = (-x - 2*c*z, 2*x + z)
            assert (Cy0[k], Cy1[k]) == gy, (q, k)
    return True

def code_min_weight_split(q):
    # q = 7 mod 8: 2 splits; mod L, omega -> root w0 of x^2 - x + c mod 2 (c even) -> w0 in {0,1}
    C = conf(q); n = q+1; c = (1+q)//4
    assert c % 2 == 0
    res = {}
    for w0 in (0, 1):
        gens = []
        bs = [[1 if t in (0, k) else 0 for t in range(n)] for k in range(1, n)] + [[2]+[0]*q]
        for b in bs:
            Cb = matvec(C, b)
            y = [((x - z)//2 + z*w0) % 2 for x, z in zip(Cb, b)]
            gens.append(sum(bit << k for k, bit in enumerate(y)))
        # row reduce over F2
        basis = []
        for g in gens:
            for bb in basis:
                g = min(g, g ^ bb)
            if g: basis.append(g)
        dim = len(basis)
        best = n
        # Gray code enumeration
        cur = 0
        for i in range(1, 1 << dim):
            j = (i & -i).bit_length() - 1
            cur ^= basis[j]
            w = bin(cur).count('1')
            if w < best: best = w
        res[w0] = (dim, best)
    return res

def code_min_weight_inert(q):
    # q = 3 mod 8: 2 inert, residue field F4 = F2[w]/(w^2+w+1); c odd
    C = conf(q); n = q+1; c = (1+q)//4
    assert c % 2 == 1
    # F4 elements as 2-bit (a + b w); vectors as pair of bitmasks (A, B)
    def mulw(A, B):  # (a + b w) * w = a w + b w^2 = a w + b (w + 1) = b + (a + b) w
        return (B, A ^ B)
    gens = []
    bs = [[1 if t in (0, k) else 0 for t in range(n)] for k in range(1, n)] + [[2]+[0]*q]
    for b in bs:
        Cb = matvec(C, b)
        A = sum((((x - z)//2) % 2) << k for k, (x, z) in enumerate(zip(Cb, b)))
        B = sum(((z) % 2) << k for k, z in enumerate(b))
        gens.append((A, B)); gens.append(mulw(A, B))
    # F2-basis of the F4-span (it's an F2-space of dim 2*dim_F4); encode as 2n-bit integers
    basis = []
    for A, B in gens:
        g = A | (B << n)
        for bb in basis: g = min(g, g ^ bb)
        if g: basis.append(g)
    dim2 = len(basis)
    mask = (1 << n) - 1
    best = n; cur = 0
    for i in range(1, 1 << dim2):
        j = (i & -i).bit_length() - 1
        cur ^= basis[j]
        w = bin((cur & mask) | (cur >> n)).count('1')
        if w < best: best = w
    return dim2 // 2, best

def transfer_test(q, hs=(2, 4, 6, 8, 10), trials=300):
    C = conf(q); n = q+1; c = (1+q)//4
    K = [row[:q] for row in C[:q]]
    out = []
    for h in hs:
        minsupp4 = q
        for _ in range(trials):
            # random primitive balanced v with ||v||_1 = h, supported on finite indices
            while True:
                v = [0]*q
                for t in range(h//2):
                    v[random.randrange(q)] += 1; v[random.randrange(q)] -= 1
                if sum(abs(x) for x in v) == h and math.gcd(*[abs(x) for x in v if x]) == 1: break
            Kv = [sum(K[i][j]*v[j] for j in range(q)) for i in range(q)]
            s4 = sum(1 for x in Kv if x % 4)
            minsupp4 = min(minsupp4, s4)
            supp = set(i for i in range(q) if v[i] or Kv[i] % 4)
            # reductions mod primes above 2: y = (Kv - v)/2 + v*omega, infinity coordinate 0
            ok = False
            if c % 2 == 0:
                for w0 in (0, 1):
                    yb = [((Kv[i]-v[i])//2 + v[i]*w0) % 2 for i in range(q)]
                    assert set(i for i in range(q) if yb[i]) <= supp
                    if any(yb): ok = True
                    if any(yb): assert sum(yb) >= math.sqrt(q) + 1, (q, "weight")
            else:
                ya = [((Kv[i]-v[i])//2) % 2 for i in range(q)]; yw = [v[i] % 2 for i in range(q)]
                for sgn in (0, 1):  # omega or 1-omega = omega + 1 mod 2
                    A = [(ya[i] + sgn*yw[i]) % 2 for i in range(q)]
                    sup = set(i for i in range(q) if A[i] or yw[i])
                    assert sup <= supp
                    if sup: ok = True; assert len(sup) >= math.sqrt(q) + 1
            assert ok
        out.append((h, minsupp4, round(math.sqrt(q)+1-h, 2)))
    return out

if __name__ == '__main__':
    random.seed(1)
    for q in (7, 11, 19, 23, 31, 43, 47):
        check_identity(q)
    print('identity C y(b) = g y(b): exact, q = 7..47')
    for q in (7, 23, 31, 47):
        print('split q', q, 'mod L above 2 (dim, min weight):', code_min_weight_split(q), 'sqrt(q)+1 =', round(math.sqrt(q)+1, 2))
    for q in (11, 19):
        print('inert q', q, 'mod 2 over F4 (dim, min weight):', code_min_weight_inert(q), 'sqrt(q)+1 =', round(math.sqrt(q)+1, 2))
    for q in (103, 127, 167, 199):
        print('transfer q', q, '(h, min supp Kv mod 4, bound):', transfer_test(q))
