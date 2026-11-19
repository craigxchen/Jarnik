"""Smallest clusters inside the VI.4 family: q in {3,27}, k, all marks (form, sign), X = t^L, FULL gcd of the
chosen subfamily removed, exact cluster test deg(z_i - z_j) <= floor(D/2). Own F_9[t] arithmetic."""
import sys, itertools, time
from marks import raw_mark
from subsets import forms_list
from f9lib import F9P, gcd9, divmod9
q, k, Lmax, mmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
Lmin = int(sys.argv[5]) if len(sys.argv) > 5 else 1
forms = forms_list(k)
marks_idx = [(f, e) for f in forms for e in (1, -1)]
best = {}
for L in range(Lmin, Lmax+1):
    t0 = time.time()
    T = F9P([0]*L + [1])
    M = [raw_mark(f, e, q, T)[0] for f, e in marks_idx]
    Draw = M[0].deg()
    nM = len(M)
    dif = {}; pg = {}
    for a, b in itertools.combinations(range(nM), 2):
        d = (M[a]-M[b]).deg()
        if d < 0: dif[(a, b)] = None; continue
        dif[(a, b)] = d
    # pairwise gcds only where needed: pair can be in a cluster only if 2d <= Draw + gcd(a,b)
    ok = {}
    for (a, b), d in dif.items():
        if d is None: ok[(a, b)] = False; continue
        if 2*d <= Draw + Draw:   # cheap prefilter impossible to fail; compute gcd
            gab = gcd9(M[a], M[b]).deg()
            pg[(a, b)] = gab
            ok[(a, b)] = 2*d <= Draw + gab
    # cliques in the pair graph, then exact subset test with subset gcd
    adj = {a: set() for a in range(nM)}
    for (a, b), v in ok.items():
        if v: adj[a].add(b); adj[b].add(a)
    found = {}
    def extend(S, cand):
        m = len(S)
        if m >= 2:
            G = S_g[tuple(S)]
            Dr = Draw - G
            if all(2*(dif[(a, b)] - G) <= Dr for a, b in itertools.combinations(S, 2)):
                if m not in found or Dr < found[m][0]:
                    found[m] = (Dr, [marks_idx[s] for s in S])
            else:
                return   # removing more points only shrinks the gcd? no: subset gcd grows; keep exploring anyway below
        if m == mmax: return
        for c in sorted(cand):
            if c <= (S[-1] if S else -1): continue
            Snew = S + [c]
            if len(Snew) >= 2:
                prev = S_g.get(tuple(S))
                G = gcd9(S_poly[tuple(S)], M[c]) if prev is not None else None
                if len(S) == 1: G = gcd9(M[S[0]], M[c])
                S_poly[tuple(Snew)] = G
                S_g[tuple(Snew)] = G.deg()
            else:
                S_poly[tuple(Snew)] = M[c]
            extend(Snew, cand & adj[c])
    S_g = {}; S_poly = {}
    extend([], set(range(nM)))
    line = 'q=%d k=%d L=%d Draw=%d pairs ok=%d ' % (q, k, L, Draw, sum(ok.values()))
    for m in sorted(found):
        line += ' m=%d:D=%d' % (m, found[m][0])
        if m not in best or found[m][0] < best[m][0]:
            best[m] = (found[m][0], L, found[m][1])
    print(line, '[%.0fs]' % (time.time()-t0), flush=True)
print('BEST', {m: best[m] for m in sorted(best)})
