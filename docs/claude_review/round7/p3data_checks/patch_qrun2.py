import sys
p = sys.argv[1]
s = open(p).read()
add = '''

def round_nearest(masks, lam, plist):
    """Give every support cut a set of distinct primes (one prime, or a product of two) whose log-sum is
    nearest to lam_s; cuts processed by decreasing weight.  Returns cols and total |deviation|."""
    import bisect
    logs = [math.log(p) for p in plist]
    used = set()
    cols, dev = [], 0.0
    order = sorted([(l, m) for m, l in zip(masks, lam) if l > 0], reverse=True)
    for l, m in order:
        best = None
        k = bisect.bisect_left(logs, l)
        for kk in range(max(0, k - 3), min(len(plist), k + 3)):
            if kk not in used:
                d = abs(logs[kk] - l)
                if best is None or d < best[0]:
                    best = (d, (kk,))
        if l >= 2 * logs[0]:
            for i in range(len(plist)):
                if logs[i] > l - logs[0]:
                    break
                if i in used:
                    continue
                tgt = l - logs[i]
                j = bisect.bisect_left(logs, tgt)
                for jj in (j - 1, j, j + 1):
                    if 0 <= jj < len(plist) and jj != i and jj not in used:
                        d = abs(logs[i] + logs[jj] - l)
                        if best is None or d < best[0]:
                            best = (d, (i, jj))
        for kk in best[1]:
            used.add(kk)
        cols.append((m, [plist[kk] for kk in best[1]]))
        dev += best[0]
    return cols, dev


def realize_safety(B, M, wmin, masks_all, lp_sup, nu, iters, tlim, exact, plist, log=None):
    """sc with safety delta, nearest-prime rounding, exact-in-float feasibility check; delta raised on failure."""
    iu = pair_index(M)
    for delta in (0.5, 1.0, 2.0, 3.0, 5.0):
        Bd = B + delta * (1 - np.eye(M))
        scd = solve_sc(Bd, M, wmin, masks_all, lp_sup, nu, max(2, iters // 2), tlim, exact)
        if scd['W'] is None:
            continue
        cols, dev = round_nearest(scd['masks'], scd['lam'], plist)
        mu = profile_margins(cols, M)
        ok = bool(((mu - B)[iu] >= 1e-9).all())
        if log:
            log('safety %.1f: sc %.2f rounded dev %.3f ok %s' % (delta, scd['W'], dev, ok))
        if ok:
            return cols, delta
    return None, None
'''
s = s.replace("\n\ndef profile_margins(cols, M):", add + "\n\ndef profile_margins(cols, M):", 1)
old = "    cols, rec['prime_floor'] = realize_fb(D, sc, plist, lp_sup, M, masks_all, nu, iters, tlim, exact)\n"
new = """    if M <= 8:
        cols, rec['prime_floor'] = realize_fb(D, sc, plist, lp_sup, M, masks_all, nu, iters, tlim, exact)
    else:
        cols, rec['prime_safety'] = realize_safety(D, M, wmin, masks_all, lp_sup, nu, iters, tlim, exact, plist, log)
"""
assert old in s; s = s.replace(old, new)
old = "    cols_s, rec['prime_sep_floor'] = realize_fb(Bs, scs, plist, lp_sup_s, M, ms if M <= 12 else None, nus, iters, tlim, exact)\n"
new = """    if M <= 8:
        cols_s, rec['prime_sep_floor'] = realize_fb(Bs, scs, plist, lp_sup_s, M, ms if M <= 12 else None, nus, iters, tlim, exact)
    else:
        cols_s, rec['prime_sep_safety'] = realize_safety(Bs, M, wmin, ms if M <= 12 else None, lp_sup_s, nus, iters, tlim, exact, plist, log)
"""
assert old in s; s = s.replace(old, new)
open(p, 'w').write(s)
print('patched')
