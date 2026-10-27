# diag.py WITNESS.pkl -- mechanism diagnostics of a GB (phi = arg pi) witness (round 15).
# For each (v,s) on which the witness pays a NON-genuine collision (q with q^a || m_(v,s) but q not dividing D_s, or
# a higher power than in D_s), report D_s(Z_v), its F-part (what the genuine system pays), and the uncharged cofactor
# U = |D_s| / F-part (rigtheory Lemma 3.1) which is what pays for the deviation.
import sys, pickle, math, itertools
from collections import Counter
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import slverify as sv
from sympy import factorint

def run(fn):
    W = pickle.load(open(fn, "rb")); P, K = W["P"], W["K"]; r = len(P); F = sorted(W["F"]); pis = [sv.gp(p) for p in P]
    V = [v for v in sv.monos(r, K) if [t for t in v if t][-1] > 0]
    mv = {}
    for q in F:
        lab = W["labels"][q]; top = lab["top"]; A = lab["A"]; m = q - (1 if q % 4 == 1 else -1); n = m // 4
        tab = sv.dlog_table(tuple(lab["g1"]), q, m); iota = tab[(0, 1)]; Mtop = m * q ** (A - 1)
        for v in V:
            d = sum(a * b for a, b in zip(v, top)) % Mtop
            if d % n: continue
            s = [u for u in range(4) if (u * iota) % m == d % m][0]; b = 1
            while b < A and d % (n * q ** b) == 0: b += 1
            mv.setdefault((v, s), {})[q] = b
    rows = []; nongen_q = Counter(); tot_nongen = 0
    for (v, s), dq in mv.items():
        X, Y = sv.Zv(pis, v); D = (Y, X - Y, X, X + Y)[s]
        Fpart = 1
        for q in F:
            while D % (q ** (1 + int(math.log(abs(D), q)) if False else 1)) == 0 and False: pass
        fac = factorint(abs(D)); Fpart = math.prod(q ** e for q, e in fac.items() if q in F)
        extra = {q: b - fac.get(q, 0) for q, b in dq.items() if b > fac.get(q, 0)}
        if extra:
            tot_nongen += 1
            for q in extra: nongen_q[q] += 1
            m_tot = math.prod(q ** b for q, b in dq.items())
            rows.append((sum(map(abs, v)), v, s, D, Fpart, abs(D) // Fpart, extra, m_tot, fac))
    rows.sort(key=lambda t: (t[0], t[1], t[2]))
    print(f"{fn}: P={P} K={K}: {len(mv)} charged (v,s) pairs, {tot_nongen} with a non-genuine collision")
    print("  non-genuine collisions per modulus:", dict(sorted(nongen_q.items())))
    usage = Counter()
    for nrm, v, s, D, Fp, U, extra, mt, fac in rows:
        usage[nrm] += 1
    print("  by norm:", dict(sorted(usage.items())))
    big = Counter()
    for nrm, v, s, D, Fp, U, extra, mt, fac in rows:
        top = max(fac) if fac else 1
        big["2N-part only" if all(p == 2 or p in P for p in fac if p not in F) and U == math.prod(p ** e for p, e in fac.items() if p == 2 or p in P) else "large prime cofactor"] += 1
    print("  source of slack:", dict(big))
    for t in rows[:25]:
        nrm, v, s, D, Fp, U, extra, mt, fac = t
        print(f"   |v|={nrm} v={v} s={s} D_s={D} = {dict(fac)}  genuine F-part={Fp}  uncharged U={U}  fake m={mt}  extra={extra}")

if __name__ == '__main__':
    for fn in sys.argv[1:]: run(fn)
