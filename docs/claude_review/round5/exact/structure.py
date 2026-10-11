"""Structure of a weighted support: pattern sizes, pair separation (weighted), decomposition of the
support into Hadamard cores (sets of M-1 balanced cuts separating every pair exactly M/2 times)."""
import sys, json, itertools
from fractions import Fraction
from core import pattern_of, rows_of, Group


def analyse(M, support, weights):
    T = [pattern_of(M, r) for r in support]
    w = [Fraction(x) for x in weights]
    sizes = sorted(min(len(r), M - len(r)) for r in support)
    sep = {}
    for x, y in itertools.combinations(range(M), 2):
        sep[(x, y)] = sum(wi for s, wi in zip(T, w) if ((s >> x) & 1) != ((s >> y) & 1))
    vals = sorted(set(sep.values()))
    out = {'M': M, '|T|': len(T), 'Omega': sum(w), 'cut sizes': sizes, 'pair separation values': [str(v) for v in vals],
           'group': Group(M, T).invs}
    # Hadamard cores inside the support (unweighted), only if M even and all cuts balanced
    cores = []
    if M % 2 == 0 and M <= 12:
        bal = [s for s in T if bin(s).count('1') in (M // 2, M // 2 - 1) and min(bin(s).count('1'), M - bin(s).count('1')) == M // 2]
        for comb in itertools.combinations(bal, M - 1):
            ok = True
            for x, y in itertools.combinations(range(M), 2):
                c = sum(1 for s in comb if ((s >> x) & 1) != ((s >> y) & 1))
                if c != M // 2:
                    ok = False; break
            if ok:
                cores.append([rows_of(s, M) for s in comb])
    out['hadamard_cores_inside'] = len(cores)
    return out, cores


if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    out, cores = analyse(d['M'], d['support'], d['weights'])
    print(json.dumps(out, default=str))
    for c in cores[:4]:
        print('  core:', c)
