import sys
p = sys.argv[1]
s = open(p).read()
old_start = s.index('def menu_realize(')
old_end = s.index('def realize(')
new_menu = '''def menu_realize(B, masks, lam, plist, extra=(), nopt=8, slack=0.25, time_limit=120, required=True, maxp=1, nsmall=12):
    """Each support cut chooses up to maxp primes among the nopt smallest primes with log p >= lam_s - slack
    plus the nsmall smallest primes overall (at least one prime if required); optional extra cuts may take
    up to maxp primes among the nsmall smallest primes; all primes distinct; mu >= B (+1e-7)."""
    M = B.shape[0]
    iu = pair_index(M)
    logs = np.log(np.array(plist, dtype=float))
    T = [(m, l) for m, l in zip(masks, lam) if l > 0]
    cuts, opts, req = [], [], []
    for m, l in T:
        idx = sorted(set([k for k in range(len(plist)) if logs[k] >= l - slack][:nopt]) | set(range(nsmall)))
        cuts.append(m); opts.append(idx); req.append(required)
    for m in extra:
        if m in cuts or nsmall == 0:
            continue
        cuts.append(m); opts.append(list(range(nsmall))); req.append(False)
    P, S = pair_coeffs(cuts, M)
    var = [(c, k) for c in range(len(cuts)) for k in opts[c]]
    nv = len(var)
    cost = np.array([logs[k] for c, k in var])
    Apair = np.zeros((len(iu[0]), nv))
    for v, (c, k) in enumerate(var):
        Apair[:, v] = P[:, c] * logs[k]
    rows = [Apair]; lo = [B[iu] + 1e-7]; hi = [np.full(len(iu[0]), np.inf)]
    Acut = np.zeros((len(cuts), nv))
    for v, (c, k) in enumerate(var):
        Acut[c, v] = 1
    rows.append(Acut); lo.append(np.array([1.0 if r else 0.0 for r in req])); hi.append(np.full(len(cuts), float(maxp)))
    used = sorted({k for c, k in var})
    Apr = np.zeros((len(used), nv))
    pos = {k: i for i, k in enumerate(used)}
    for v, (c, k) in enumerate(var):
        Apr[pos[k], v] = 1
    rows.append(Apr); lo.append(np.full(len(used), -np.inf)); hi.append(np.ones(len(used)))
    res = milp(cost, constraints=[LinearConstraint(csr_matrix(np.vstack(rows)), np.concatenate(lo), np.concatenate(hi))],
               integrality=np.ones(nv), bounds=Bounds(np.zeros(nv), np.ones(nv)),
               options=dict(time_limit=time_limit, disp=False))
    if res.x is None:
        return None
    cols = {}
    for v, (c, k) in enumerate(var):
        if res.x[v] > 0.5:
            cols.setdefault(cuts[c], []).append(plist[k])
    return list(cols.items())


'''
s = s[:old_start] + new_menu + s[old_end:]
old_r = s[s.index('def realize('):s.index('def realize_fb(')]
new_r = '''def realize(B, sc, plist, lp_support, time_limit=120):
    attempts = ((dict(nopt=8, slack=0.25, required=True, maxp=1, nsmall=0), ()),
                (dict(nopt=10, slack=0.6, required=False, maxp=2, nsmall=12), ()),
                (dict(nopt=12, slack=1.0, required=False, maxp=3, nsmall=24), tuple(lp_support)))
    for kw, extra in attempts:
        cols = menu_realize(B, sc['masks'], sc['lam'], plist, extra=extra, time_limit=time_limit, **kw)
        if cols is not None:
            return cols
    return None


'''
s = s.replace(old_r, new_r)
s = s.replace("def run_instance(args):\n    M, t, Xf, iters, tlim = args\n    t0 = time.time()",
              "def run_instance(args):\n    M, t, Xf, iters, tlim = args\n    t0 = time.time()\n    def log(msg):\n        with open('q_progress_M%d.txt' % M, 'a') as fp:\n            fp.write('t=%d X=%d %.0fs %s\\n' % (t, Xf * M, time.time() - t0, msg))")
for key, msg in (("    exact = (M <= 8)\n", "    log('lp done')\n    exact = (M <= 8)\n"),
                 ("    rec['sc_ncuts'] = int((sc['lam'] > 0).sum())\n", "    rec['sc_ncuts'] = int((sc['lam'] > 0).sum())\n    log('sc done')\n"),
                 ("    # placement for the no-sep prime profile", "    log('realised: %s' % (cols is not None))\n    # placement for the no-sep prime profile"),
                 ("    rec['sep_positions'] = name; rec['Wcut_sep'] = Wcs\n", "    rec['sep_positions'] = name; rec['Wcut_sep'] = Wcs\n    log('sep lp done')\n"),
                 ("    rec['unif_sep_mu'] = unif_sep(D)", "    log('sep realised: %s' % (cols_s is not None))\n    rec['unif_sep_mu'] = unif_sep(D)")):
    assert key in s, key
    s = s.replace(key, msg, 1)
open(p, 'w').write(s)
print('patched')
