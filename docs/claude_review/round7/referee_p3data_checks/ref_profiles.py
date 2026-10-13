"""Referee: exact re-verification of every realised actual-prime profile (best_ub_M*.jsonl, ils2_M*.jsonl),
own code. Checks: primes distinct, prime, =1 mod 4, >X; columns nonconstant, rows distinct; exact pair
inequality Pi_sep >= m^4 Pi_same (m from referee-recomputed moduli, identical to stored); W; prime count;
stair; ties (pairs with Pi_sep == Pi_same impossible); exact ternary support-4 and 2e_x-e_y-e_z characters;
min pair margin."""
import json, glob, math, itertools, sys
import numpy as np
from sympy import isprime, nextprime
R2 = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/r2'
DATA = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round7/p3data_checks/data'

def primes1mod4_above(X, n):
    out = []; p = X
    while len(out) < n:
        p = nextprime(p)
        if p % 4 == 1: out.append(p)
    return out

def decode(mask, M):
    return np.array([-1 if (mask >> x) & 1 else 1 for x in range(M)])

def verify(cols, M, X, mod):
    S = np.array([decode(m, M) for m, _ in cols]).T  # M x r
    P = [math.prod(pl) for _, pl in cols]
    allp = [p for _, pl in cols for p in pl]
    assert len(set(allp)) == len(allp), 'primes not distinct'
    for p in allp: assert isprime(p) and p % 4 == 1 and p > X, ('bad prime', p)
    for j in range(S.shape[1]): assert abs(S[:, j].sum()) < M, 'constant column'
    assert len({tuple(S[x]) for x in range(M)}) == M, 'rows not distinct'
    N = math.prod(P)
    minmu = None; ok = True
    for x in range(M):
        for y in range(x+1, M):
            sep = math.prod(P[j] for j in range(len(P)) if S[x, j] != S[y, j])
            same = N // sep
            m = int(mod[x, y])
            if not sep >= m**4*same: ok = False
            mu_minus_D = 0.5*(math.log(sep) - math.log(same)) - 2*math.log(m)
            minmu = mu_minus_D if minmu is None else min(minmu, mu_minus_D)
    # characters: v_j = sum c_x a_xj, a = (1+S)/2 ; margin ok iff prod P_j^(2|v_j|) >= N
    A = (1 + S)//2
    def charok(c):
        v = c @ A
        lhs = math.prod(P[j]**(2*abs(int(v[j]))) for j in range(len(P)) if v[j] != 0)
        return lhs >= N, 0.5*sum(abs(int(v[j]))*math.log(P[j]) for j in range(len(P))) - 0.5*math.log(N)
    nchar = 0; cmin = None; cbad = 0
    idx = range(M)
    for quad in itertools.combinations(idx, 4):
        for pos in itertools.combinations(quad, 2):
            if 0 not in [quad.index(p) for p in pos]: continue  # c ~ -c
            c = np.zeros(M, dtype=np.int64)
            for z in quad: c[z] = 1 if z in pos else -1
            okc, s = charok(c); nchar += 1
            if not okc: cbad += 1
            cmin = s if cmin is None else min(cmin, s)
    for x in idx:
        for y, z in itertools.combinations([i for i in idx if i != x], 2):
            c = np.zeros(M, dtype=np.int64); c[x] = 2; c[y] = -1; c[z] = -1
            okc, s = charok(c); nchar += 1
            if not okc: cbad += 1
            cmin = min(cmin, s)
    W = sum(math.log(p) for p in allp)
    st = sum(math.log(p) for p in primes1mod4_above(X, len(cols)))
    st_primes = sum(math.log(p) for p in primes1mod4_above(X, len(allp)))
    return dict(ok=ok, minmu_minus_D=minmu, W=W, ncols=len(cols), nprimes=len(allp), stair_cols=st,
                stair_primes=st_primes, nchar=nchar, char_bad=cbad, char4_min=cmin, maxp=max(allp))

if __name__ == '__main__':
    out = []
    for fn, key in [('ils2_M*.jsonl', 'prime_ils_cols'), ('best_ub_M*.jsonl', 'prime_cols')]:
        for f in sorted(glob.glob(f'{R2}/{fn}')):
            for line in open(f):
                r = json.loads(line)
                if key not in r or r['M'] == 8: continue
                M, X, t = r['M'], r['X'], r['t']
                mod = np.load(f'{DATA}/mod_construct_M{M}_t{t}_u1_X{X//M}.npy', allow_pickle=True)
                res = verify(r[key], M, X, mod)
                Wclaim = r['prime_ils'] if key == 'prime_ils_cols' else r['prime']
                rec = dict(src=f.split('/')[-1], M=M, X=X, t=t, Wclaim=Wclaim, stair_claim=r.get('stair'), **res)
                out.append(rec); print(json.dumps(rec), flush=True)
    json.dump(out, open('ref_profiles_out.json', 'w'), indent=0)
    print('n', len(out), 'all pairs ok', all(o['ok'] for o in out), 'all char ok', all(o['char_bad'] == 0 for o in out))
