def pool(maxsupp, maxentry, minsupp=2):
    out = []
    for s in range(minsupp, maxsupp + 1):
        for S in itertools.combinations(range(12), s):
            for vals in itertools.product([e for e in range(-maxentry, maxentry + 1) if e], repeat=s):
                if sum(vals) != 0 or vals[0] < 0: continue
                g = 0
                for v in vals: g = gcd(g, abs(v))
                if g != 1: continue
                z = [0]*12
                for x, v in zip(S, vals): z[x] = v
                out.append(z)
    return np.array(out, dtype=np.int64)
