"""Sweep driver: for M, t in [tmin, tmax], all s = t mod 2, ..., t (s >= 0 suffices: k -> -k maps
S^M(p,n) onto S^M(n,p)), all partitions lam of M: isotypic regularity reg_lam(S_{s,t}) mod prime.

usage: python3 driver.py M tmin tmax prime_index(1|2) outfile [kind=shell|slice] [smax] [smin]
Appends one JSON line per (M, s, t, kind, prime); already-present entries are skipped."""
import sys, json, time, os
from combi import partitions, shell_orbits, slice_orbits, phi_floor, hook_dim
from iso import lam_profile, lam_profile_parity, LamData
from modla import P1, P2


def main():
    M, tmin, tmax, pi = map(int, sys.argv[1:5])
    out = sys.argv[5]
    kind = sys.argv[6] if len(sys.argv) > 6 else 'shell'
    smax = int(sys.argv[7]) if len(sys.argv) > 7 else 10 ** 9
    smin = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    p = P1 if pi == 1 else P2
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                r = json.loads(line)
                done.add((r['M'], r['s'], r['t'], r['kind'], r['prime']))
            except Exception:
                pass
    lds = {lam: LamData(lam) for lam in partitions(M)}
    for t in range(tmin, tmax + 1):
        for s in range(t % 2, min(t, smax) + 1, 2):
            if s < smin:
                continue
            key = (M, s, t, kind, p)
            if key in done:
                continue
            pp, nn = (t + s) // 2, (t - s) // 2
            orbs = shell_orbits(M, pp, nn) if kind == 'shell' else slice_orbits(M, pp, nn)
            dg = phi_floor(M, pp, nn)
            T = time.time()
            comps = []
            size = 0
            for lam in partitions(M):
                r = None
                if s == 0:
                    try:
                        r = lam_profile_parity(lam, orbs, p, ld=lds[lam], dguess=dg)
                    except Exception as ex:
                        print('parity variant failed, falling back:', lam, ex, flush=True)
                        r = None
                if r is None:
                    r = lam_profile(lam, orbs, p, ld=lds[lam], dguess=dg)
                size += hook_dim(lam) * r['N']
                if r['N'] == 0:
                    continue
                comps.append({'lam': list(lam), 'N': r['N'], 'reg': r['reg'],
                              'profile': r['profile'], 'ngen': r.get('ngen'),
                              'beyond_phi': r.get('beyond_guess', False),
                              'time': round(r['time'], 2), 'parity': r.get('parity', False),
                              'reg_even': r.get('reg_even'), 'reg_odd': r.get('reg_odd')})
            rec = {'M': M, 's': s, 't': t, 'p': pp, 'n': nn, 'kind': kind, 'prime': p,
                   'size': size, 'norb': len(orbs), 'phi_floor': dg,
                   'reg': max(c['reg'] for c in comps),
                   'argmax': [c['lam'] for c in comps if c['reg'] == max(x['reg'] for x in comps)],
                   'comps': comps, 'time': round(time.time() - T, 2)}
            with open(out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
            print(f"M={M} t={t} s={s} {kind} p={p}: |A|={size} reg={rec['reg']} ratio={rec['reg']/t:.4f} "
                  f"phi={dg} argmax={rec['argmax']} time={rec['time']}", flush=True)


if __name__ == '__main__':
    main()
