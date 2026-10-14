"""J6 test: is the leading-form ideal J(A) of a T=1 box generated (in degrees <= m) by powers of
'width forms' <k',v>^{h_A(v)+1} (v integer direction, h_A(v) = width)?  Call that ideal J_w.
J_w is contained in J(A) always (prod_{t}(<k',v> - t) over the h+1 values vanishes on A).
Mod-p logic: cumulative dim_p(J_w)_{<=i} <= dim_Q(J_w)_{<=i} <= dim_Q J(A)_{<=i} = N_{<=i} - rank_Q(eval_{<=i})
<= N_{<=i} - rank_p(eval_{<=i}).  So equality of the outer two (mod p) proves J_w = J(A) in all degrees <= i over Q.
"""
import itertools, sys, math, time
import numpy as np
from flint import nmod_mat
PR = 2147483647

def monos(d, i):
    return [e for e in itertools.combinations_with_replacement(range(d), i)]

def mono_tuple(c, d):
    e = [0]*d
    for x in c: e[x] += 1
    return tuple(e)

def eval_rank_profile(A, d, m):
    """rank_p of evaluation of all monomials of degree <= i on A, i = 0..m (cumulative)"""
    A = np.array(A, dtype=np.int64)
    cols = []; ranks = []
    rows = []
    for i in range(m+1):
        for c in itertools.combinations_with_replacement(range(d), i):
            col = np.ones(len(A), dtype=np.int64)
            for x in c: col = col * A[:, x]
            rows.append([int(t) % PR for t in col])
        M = nmod_mat(len(rows), len(A), [v for r in rows for v in r], PR)
        ranks.append(M.rank())
    return ranks

def poly_pow_linear(v, e, d):
    """(sum v_x k_x)^e as dict exponent-tuple -> coeff mod PR"""
    res = {tuple([0]*d): 1}
    lin = {}
    for x in range(d):
        if v[x] % PR:
            t = [0]*d; t[x] = 1; lin[tuple(t)] = v[x] % PR
    for _ in range(e):
        new = {}
        for a, ca in res.items():
            for b, cb in lin.items():
                k = tuple(p+q for p, q in zip(a, b))
                new[k] = (new.get(k, 0) + ca*cb) % PR
        res = {k: c for k, c in new.items() if c}
    return res

def widths(A, d, R, hmax):
    A = np.array(A, dtype=np.int64)
    out = []
    rng = range(-R, R+1)
    vs = []
    for v in itertools.product(rng, repeat=d):
        nz = [t for t in v if t]
        if not nz or nz[0] < 0: continue
        if math.gcd(*[abs(t) for t in nz]) != 1: continue
        vs.append(v)
    V = np.array(vs, dtype=np.int64)
    for s in range(0, len(V), 20000):
        blk = V[s:s+20000]
        P = A @ blk.T
        h = P.max(axis=0) - P.min(axis=0)
        for t in np.nonzero(h <= hmax)[0]:
            out.append((tuple(int(z) for z in blk[t]), int(h[t])))
    return out

def jw_dims(gens, d, m):
    """dims of (J_w)_i for i=0..m, gens = list of (v,h): generator <k',v>^{h+1}"""
    dims = []; basis_prev = None
    for i in range(m+1):
        mi = [mono_tuple(c, d) for c in itertools.combinations_with_replacement(range(d), i)]
        idx = {e: t for t, e in enumerate(mi)}
        rows = []
        if basis_prev is not None:
            for r in basis_prev:
                for x in range(d):
                    row = {}
                    for e, c in r.items():
                        e2 = list(e); e2[x] += 1
                        row[tuple(e2)] = c
                    rows.append(row)
        for v, h in gens:
            if h + 1 == i:
                rows.append(poly_pow_linear(v, i, d))
        if not rows:
            dims.append(0); basis_prev = []; continue
        M = nmod_mat(len(rows), len(mi), [0]*(len(rows)*len(mi)), PR)
        for a, r in enumerate(rows):
            for e, c in r.items(): M[a, idx[e]] = c
        Rr, rk = M.rref()
        dims.append(rk)
        basis_prev = []
        for a in range(rk):
            row = {}
            for b in range(len(mi)):
                c = int(Rr[a, b])
                if c: row[mi[b]] = c
            basis_prev.append(row)
    return dims

def run(name, A, d, m, R=2):
    t0 = time.time()
    ranks = eval_rank_profile(A, d, m)
    H = [ranks[0]] + [ranks[i]-ranks[i-1] for i in range(1, m+1)]
    N = [math.comb(i+d-1, d-1) for i in range(m+1)]
    JA = [N[i]-H[i] for i in range(m+1)]
    gens = widths(A, d, R, m-1)
    Jw = jw_dims(gens, d, m)
    cumA = list(itertools.accumulate(JA)); cumW = list(itertools.accumulate(Jw))
    eq = [cumA[i] == cumW[i] for i in range(m+1)]
    s = (f"{name}: |A|={len(A)} d={d} m={m} R={R} #width-gens={len(gens)}\n  H_A    ={H}\n  dimJ(A)={JA}\n  dimJ_w ={Jw}\n"
         f"  cumulative equal (=> J(A)=J_w over Q up to that degree): {eq}   [{time.time()-t0:.1f}s]")
    print(s, flush=True)
    return s

def onecut(d, S, L):
    return [k for k in itertools.product((-1,0,1), repeat=d) if abs(sum(k[x] for x in S)) <= L]

def sylvester8(L):
    # points x in F_2^3, base = x=7 (last); k' = coords of other 7 points; k_base = -sum k'
    pts = list(range(8)); base = 7; others = [x for x in pts if x != base]
    def dot(a, x): return bin(a & x).count('1') % 2
    A = []
    for kp in itertools.product((-1,0,1), repeat=7):
        k = dict(zip(others, kp)); k[base] = -sum(kp)
        ok = all(abs(sum(k[x]*(1 if dot(a,x)==0 else -1) for x in pts)) <= 2*L for a in range(1, 8))
        if ok: A.append(kp)
    return A

if __name__ == '__main__':
    out = []
    which = sys.argv[1] if len(sys.argv) > 1 else 'small'
    if which == 'small':
        out.append(run('onecut d=5 S=all L=1', onecut(5, range(5), 1), 5, 6))
        out.append(run('onecut d=5 S=all L=2', onecut(5, range(5), 2), 5, 7))
        out.append(run('onecut d=6 S={0,1} L=1 (zetc n=7 example)', onecut(6, [0,1], 1), 6, 8))
        A2 = [k for k in itertools.product((-1,0,1), repeat=5) if abs(k[0]+k[1]+k[2]) <= 1 and abs(k[2]+k[3]+k[4]) <= 1]
        out.append(run('twocut d=5 S1={0,1,2} S2={2,3,4} L=1,1', A2, 5, 6))
        A3 = [k for k in itertools.product((-1,0,1), repeat=6) if abs(k[0]+k[1]+k[2]+k[3]) <= 1 and abs(k[2]+k[3]+k[4]+k[5]) <= 1]
        out.append(run('twocut d=6 S1={0..3} S2={2..5} L=1,1', A3, 6, 7))
    elif which == 'M8':
        A = sylvester8(2)
        out.append(run('M8 = Syl(8), L=2', A, 7, 8, R=2))
    elif which == 'M8L1':
        A = sylvester8(1)
        out.append(run('Syl(8), L=1', A, 7, 4, R=2))
    open(f"jw_test_{which}_out.txt", 'w').write("\n".join(out) + "\n")
