import numpy as np, itertools, random, math
from plib import chi_table, patterns, N_def, N_formula
random.seed(3)
def circ(q):
    ch = chi_table(q); return np.array([[ch[(a - t) % q] for a in range(q)] for t in range(q)], dtype=np.int64)
def sig(L, q, C):
    f = L[:, :q]; F = f.sum(1, keepdims=True); h = F + f; c = f @ C
    return np.maximum(np.abs(h), np.abs(c)).sum(1) - q
# (a) normal form completeness: brute-force all 4-term F=0 (no symmetry) at q=19,23
for q in [19, 23]:
    C = circ(q); rows = []
    for a, b in itertools.combinations(range(q), 2):
        for c, d in itertools.combinations([x for x in range(q) if x not in (a, b)], 2):
            v = [0]*(q+1); v[a] = v[b] = 1; v[c] = v[d] = -1; rows.append(v)
    L = np.array(rows, dtype=np.int64); s = sig(L, q, C)
    print(f"(a) q={q}: brute-force 4-term F=0 min sigma = {s.min()} over {len(L)} characters; values {sorted(set(s.tolist()))[:6]}...", flush=True)
# (b) Lemma 2.3 tests: zeros of f*chi vs (q-1)/2 - j0
def j0(f, q):
    T = [t for t in range(q) if f[t] % q]
    for j in range(q):
        if sum(f[t]*pow(t, j, q) for t in range(q) if f[t]) % q: return j
    return None
viol = 0; tests = 0; tight = 0
for q in [7, 11, 19, 23, 31, 43, 59, 83, 103, 131]:
    C = circ(q)
    for trial in range(400 if q > 11 else 2000):
        k = random.randint(1, (q+1)//2)
        T = random.sample(range(q), k)
        f = [0]*q
        for t in T: f[t] = random.choice([1, -1, 2, -2, 3, q, -q, 2*q+1])
        if all(f[t] % q == 0 for t in T): continue
        conv = np.array(f, dtype=np.int64) @ C
        z = int((conv == 0).sum()); J = j0(f, q); b = (q-1)//2 - J
        tests += 1; viol += (z > b); tight += (z == b)
print(f"(b) Lemma 2.3: {tests} random f (q=7..131, |T|<=(q+1)/2, values incl. multiples of q): violations {viol}, equality cases {tight}", flush=True)
# exhaustive ternary f at q=11 with |T| <= 6
q = 11; C = circ(q); viol = 0; cnt = 0
for v in itertools.product([-1, 0, 1], repeat=q):
    k = sum(1 for x in v if x)
    if k == 0 or k > (q+1)//2: continue
    z = int((np.array(v) @ C == 0).sum()); J = j0(list(v), q); cnt += 1; viol += (z > (q-1)//2 - J)
print(f"(b) Lemma 2.3 exhaustive ternary q=11, 1<=|T|<=6: {cnt} f, violations {viol}", flush=True)
# what happens beyond |T| = (q+1)/2 (hypothesis needed?): max zeros of f*chi with F = 0, ternary, |T| = (q+3)/2 at q=11
mx = 0
for v in itertools.product([-1, 0, 1], repeat=q):
    if sum(1 for x in v if x) == (q+3)//2 and sum(v) == 0:
        mx = max(mx, int((np.array(v) @ C == 0).sum()))
print(f"(b') q=11, F=0 ternary |T|=(q+3)/2=7: max #zeros = {mx} vs (q-3)/2 = {(q-3)//2}", flush=True)
# (c) dual pairs from the pattern definition
for q in [7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    ch = chi_table(q); f = [0]*q
    for x in range(2, q): f[x] = (ch[x] - ch[(x-1) % q])//2
    f[1] = 1; lam = f + [0]
    P = patterns(q); s = N_def(lam, P) - q
    assert all(abs(x) <= 1 for x in lam) and sum(lam) == 0
    print(f"(c) q={q}: dual pair sigma(def) = {s}, (q-1)/2 = {(q-1)//2}, ||lam||_1 = {sum(map(abs, lam))}, (q+1)/2 = {(q+1)//2}", flush=True)
for q in [3251, 3259, 5479]:
    ch = chi_table(q); f = [0]*q
    for x in range(2, q): f[x] = (ch[x] - ch[(x-1) % q])//2
    f[1] = 1; C = circ(q); L = np.array([f + [0]], dtype=np.int64)
    print(f"(c) q={q}: dual pair sigma(formula) = {sig(L, q, C)[0]}, (q-1)/2 = {(q-1)//2}, ||.||_1 = {int(np.abs(L).sum())}", flush=True)
# (d) base 4-terms e_x + e_y - e_z - e_inf: set of sigma values
for q in [7, 11, 19, 23, 31, 43, 59, 83, 103]:
    C = circ(q); rows = []
    for x, y in itertools.combinations(range(q), 2):
        for z in range(q):
            if z in (x, y): continue
            v = [0]*(q+1); v[x] = v[y] = 1; v[z] = -1; v[q] = -1; rows.append(v)
    s = sig(np.array(rows, dtype=np.int64), q, C)
    print(f"(d) q={q}: base 4-term sigma values {sorted(set(s.tolist()))}, (q-1)/2 = {(q-1)//2}", flush=True)
