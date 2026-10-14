# EVIDENCE ONLY: images' J spaces computed mod p. Do 6-dim projections (kernel z in {-1,0,1}^7, all GL(3,2)-images)
# of Y_3(2) close the degree-7 gap over lines+quads+width powers?
import itertools, sys, time, pickle
from rlib import *
p = 2147483647
Y, h, dJ6, dJ7 = pickle.load(open('fano_Y.pkl', 'rb'))
gens = pickle.load(open('fano_gens.pkl', 'rb'))
cap = [5] * 7
D = 7
cols = monos_eq(7, D, cap)
cidx = {c: i for i, c in enumerate(cols)}
avec = [a for a in itertools.product((0, 1), repeat=3) if any(a)]


def mat_apply(Mx, a):
    return tuple(sum(Mx[r][c] * a[c] for c in range(3)) % 2 for r in range(3))


perms = []
for Mx in itertools.product(itertools.product((0, 1), repeat=3), repeat=3):
    imgs = [mat_apply(Mx, a) for a in avec]
    if len(set(imgs)) == 7 and all(any(x) for x in imgs):
        perms.append([avec.index(x) for x in imgs])


def J_space_p(X, e):
    ms = monos_le(len(X[0]), e)
    M = flint.nmod_mat([[evalmono(q, m) % p for m in ms] for q in X], p)
    N, nul = M.nullspace()
    top = [i for i, m in enumerate(ms) if sum(m) == e]
    if nul == 0:
        return []
    R = flint.nmod_mat([[int(N[i, c]) for i in top] for c in range(nul)], p)
    Rr, rk = R.rref()
    return [{ms[top[c]]: int(Rr[r, c]) for c in range(len(top)) if int(Rr[r, c])} for r in range(rk)]


def rows_to_list(rows):
    out = []
    for r in rows:
        v = [0] * len(cols)
        for k, x in r.items():
            v[cidx[k]] = x % p
        out.append(v)
    return out


def basis_of(lst):
    if not lst:
        return []
    M = flint.nmod_mat(lst, p)
    R, rk = M.rref()
    return [[int(R[i, c]) for c in range(len(cols))] for i in range(rk)]


B = basis_of(rows_to_list(ideal_degree_rows(gens, 7, D, cap)))
print("start rank", len(B), "target", dJ7)
sys.stdout.flush()
zs = [z for z in itertools.product((-1, 0, 1), repeat=7) if any(z)]
seen = set()
reps = []
for z in zs:
    if z in seen:
        continue
    orb = set()
    for pm in perms:
        w = [0] * 7
        for i in range(7):
            w[pm[i]] = z[i]
        orb.add(tuple(w))
        orb.add(tuple(-x for x in w))
    seen |= orb
    reps.append(z)
reps.sort(key=lambda z: sum(abs(x) for x in z))
print("orbit reps", len(reps))
sys.stdout.flush()
for z in reps:
    T0 = time.time()
    K, nul = flint.fmpz_mat([list(z)]).nullspace()
    F = [[int(K[r, c]) for r in range(7)] for c in range(nul)]
    img = sorted(set(tuple(sum(a * b for a, b in zip(v, y)) for v in F) for y in Y))
    newg = {}
    for e in range(1, D + 1):
        sp = J_space_p(img, e)
        if sp:
            newg[e] = [pullback(g, F, 7) for g in sp]
    base = ideal_degree_rows(newg, 7, D, cap)
    before = len(B)
    for pm in perms:
        rows = []
        for r in base:
            nr = {}
            for ex, v in r.items():
                ne = [0] * 7
                for i, a in enumerate(ex):
                    ne[pm[i]] = a
                nr[tuple(ne)] = v
            rows.append(nr)
        B = basis_of(B + rows_to_list(rows))
        if len(B) >= dJ7:
            break
    print("z", z, "|img|", len(img), "gens", sum(len(v) for v in newg.values()), "rank", before, "->", len(B),
          "/", dJ7, "%.0fs" % (time.time() - T0))
    sys.stdout.flush()
    if len(B) >= dJ7:
        print("CLOSED (mod p evidence)")
        break
