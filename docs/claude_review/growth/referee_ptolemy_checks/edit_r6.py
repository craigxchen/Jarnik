src=open('rchart6.py').read()
old_start=src.index("def analyse(S, z):")
old_end=src.index('if __name__ == "__main__":')
new='''def analyse(S, z):
    """Return ((sep, lead_sep, rmin, special, sep_big, rem), None) or (None, reason).  Pair points
    are obtained by polynomial division of the numerator of R_i - R_j by the known factors."""
    special = {}
    for i in rows:
        special[frozenset(r for r in rows if r != i)] = z[S['idx_a'][i]]
    for T in triples + quads:
        special[T] = point(S, z, T)
    big = list(special.values())
    sep_big = min(abs(x - y) for x, y in itertools.combinations(big, 2))
    worst_rem = 0.0
    for i, j in itertools.combinations(rows, 2):
        ai, aj = z[S['idx_a'][i]], z[S['idx_a'][j]]
        ri = 1.0 if i == 1 else z[S['idx_r'][i]]
        rj = 1.0 if j == 1 else z[S['idx_r'][j]]
        qi = np.zeros(S['nq']) if i == 1 else z[S['idx_q'][i]:S['idx_q'][i] + S['nq']]
        qj = np.zeros(S['nq']) if j == 1 else z[S['idx_q'][j]:S['idx_q'][j] + S['nq']]
        num = np.polyadd(np.polymul(np.polysub(qi, qj), np.polymul([1, -ai], [1, -aj])),
                         np.polysub(ri * np.array([1, -aj]), rj * np.array([1, -ai])))
        known = [special[T] for T in triples + quads if i in T and j in T]
        den = np.poly(known)
        quo, rem = np.polydiv(num, den)
        worst_rem = max(worst_rem, np.max(np.abs(rem)) / (1 + np.max(np.abs(num))))
        if len(quo) != 2 or abs(quo[0]) < 1e-9:
            return None, f'pair {i}{j}: degree drop'
        special[frozenset((i, j))] = -quo[1] / quo[0]
    vals = list(special.values())
    sep = min(abs(x - y) for x, y in itertools.combinations(vals, 2))
    leads = [0.0] + [z[S['idx_q'][i]] for i in rows[1:]]
    lsep = min(abs(x - y) for x, y in itertools.combinations(leads, 2))
    rmin = min([1.0] + [abs(z[S['idx_r'][i]]) for i in rows[1:]])
    return (sep, lsep, rmin, special, sep_big, worst_rem), None

'''
src=src[:old_start]+new+src[old_end:]
old_main=src[src.index('if __name__ == "__main__":'):]
new_main='''if __name__ == "__main__":
    trials = int(sys.argv[1]); seed = int(sys.argv[2])
    rng = np.random.default_rng(seed)
    good = 0
    for tr in range(trials):
        S = build(rng)
        z, nf = newton(S, S['z0'])
        if nf > 1e-9:
            print(f"trial {tr}: no convergence ({nf:.1e})", flush=True); continue
        out, why = analyse(S, z)
        np.save(f"r6sol_seed{seed}_tr{tr}.npy", z)
        if out is None:
            print(f"trial {tr}: converged; {why}", flush=True); continue
        sep, lsep, rmin, special, sep_big, rem = out
        ok = sep > 1e-5 and lsep > 1e-5 and rmin > 1e-5 and rem < 1e-8
        good += ok
        print(f"trial {tr}: converged; min sep {sep:.2e} (big {sep_big:.2e}), lead sep {lsep:.2e}, min|r| {rmin:.2e}, rem {rem:.1e} -> {'NONDEGENERATE' if ok else 'degenerate'}", flush=True)
    print("nondegenerate:", good, "/", trials)
'''
src=src.replace(old_main,new_main)
open('rchart6.py','w').write(src)
