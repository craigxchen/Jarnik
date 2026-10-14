s = open("cyc_rank.py").read()
old_start = s.index("        seen = set(); mons = []; tries = 0")
old_end = s.index("        w = np.array([pow(z")
s = s[:old_start] + "        mons = stratified(d, D, eps, rows + extra, perm, q, rng)\n" + s[old_end:]
s = s.replace("def certify(", '''from math import comb
import itertools as _it
def stratified(d, D, eps, need, perm, q, rng):
    """distinct sigma-orbit representatives of exponent vectors in {0,1,2}^d, degree <= D, parity eps,
    spread over degrees: equal quota per degree, all orbits of a degree taken when there are few."""
    degs = [t for t in range(D + 1) if t % 2 == eps]
    cnt = {t: sum(comb(d, j) * comb(d - j, t - 2*j) for j in range(0, t//2 + 1) if t - 2*j <= d - j) for t in degs}
    def canon(e):
        ee = e.copy(); best = tuple(ee)
        for _ in range(q - 1):
            ee = ee[perm]; best = min(best, tuple(ee))
        return best
    quota = {}; remaining = need; open_ = list(degs)
    while open_:
        share = -(-remaining // len(open_))
        small = [t for t in open_ if cnt[t] // q + 1 <= share]
        if not small:
            for t in open_: quota[t] = share
            break
        for t in small:
            quota[t] = cnt[t]; remaining -= cnt[t] // q + 1; open_.remove(t)
    mons = []
    for t in degs:
        seen = set(); want = quota[t]
        if cnt[t] <= 20000:
            for twos in range(0, t//2 + 1):
                ones = t - 2*twos
                if ones + twos > d: continue
                for P2 in _it.combinations(range(d), twos):
                    rest = [x for x in range(d) if x not in P2]
                    for P1 in _it.combinations(rest, ones):
                        e = np.zeros(d, dtype=np.int64); e[list(P2)] = 2; e[list(P1)] = 1
                        k = canon(e)
                        if k not in seen: seen.add(k); mons.append(np.array(k))
        else:
            w = np.array([comb(d, j) * comb(d - j, t - 2*j) if t - 2*j <= d - j else 0 for j in range(t//2 + 1)], float)
            tries = 0
            while len(seen) < want and tries < 50 * want:
                tries += 1
                j = rng.choice(len(w), p=w / w.sum())
                pos = rng.permutation(d); e = np.zeros(d, dtype=np.int64)
                e[pos[:j]] = 2; e[pos[j:j + t - 2*j]] = 1
                k = canon(e)
                if k not in seen: seen.add(k); mons.append(np.array(k))
    return mons
def certify(''', 1)
open("cyc_rank.py", "w").write(s)
