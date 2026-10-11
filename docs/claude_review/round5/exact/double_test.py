"""Double-Hadamard supports: T = core(H) cup core(P H) for a row permutation P.
Prints group G_T and, if small, exact verification of sigma >= rho1 ||lambda||_1/2 at unit weights."""
import sys, random, time
from fractions import Fraction
from fixed import sylvester, paley, core_patterns
from core import Group, pattern_of, rows_of
from verify import verify
from inner import analyse, find_hom

name = sys.argv[1]; trials = int(sys.argv[2]); rho1 = Fraction(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
H = {'S8': lambda: sylvester(3), 'P8': lambda: paley(7), 'P12': lambda: paley(11), 'S16': lambda: sylvester(4)}[name]()
M = len(H)
c1 = core_patterns(H)
rng = random.Random(seed)
for t in range(trials):
    perm = list(range(M)); rng.shuffle(perm)
    H2 = [H[perm[x]] for x in range(M)]
    # renormalise H2 so first column is all ones (row sign changes)
    H2 = [[v * r[0] for v in r] for r in H2]
    c2 = core_patterns(H2)
    tries = 0
    while set(c1) & set(c2) and tries < 1000:
        tries += 1
        rng.shuffle(perm)
        H2 = [H[perm[x]] for x in range(M)]
        H2 = [[v * r[0] for v in r] for r in H2]
        c2 = core_patterns(H2)
    if set(c1) & set(c2):
        print('no disjoint second core found'); break
    T = c1 + c2
    G = Group(M, T)
    t0 = time.time()
    if not G.full_rank:
        print('trial', t, 'rank deficient'); continue
    if G.order > 5000:
        print('trial', t, 'group', G.invs, 'too large'); continue
    w = [Fraction(1)] * len(T)
    A = analyse(M, T, w, 0, rho1)
    if A.get('inconclusive'):
        print('trial', t, 'inconclusive coercivity'); continue
    eta = find_hom(G, list(A['bad'].keys())) if not A['int_viol'] else None
    print('trial', t, '|T|', len(T), 'group', G.invs, 'int_viol', len(A['int_viol']),
          'bad rational classes', len(A['bad']), 'avoiding hom', eta, 'min slack', A['min_slack'][0] if A['min_slack'] else None,
          f'{time.time()-t0:.1f}s', flush=True)
