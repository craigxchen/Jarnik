"""Exact evaluation of Hadamard-type supports in the normalised problem.

Supports (H normalised Hadamard of order M, labels a = 1..M-1, S_a = {x : H(x,a) = -1}):
  'core'      : {S_a}
  'twin-id'   : core + flipped patterns S_{a(x)} xor {x}, assignment a(x) = x for x >= 1 (Paley-I
                indexing gives the skew pairing), a(0) = 1
  'twin-rnd'  : core + flipped patterns for a random capacity-two assignment
For each support, V(T; rho0, rho1) = min over windings and weights (om_S >= 1) is computed with
fixed.evaluate_support (exact feasibility of the reported weights, exact LP-dual lower bound for the
final constraint pool).

usage: python3 hadamard_eval.py name rho0 rho1      (name in S4, S8, P8, P12)
"""
import sys, json, random
from fractions import Fraction
from fixed import evaluate_support, sylvester, paley, core_patterns, flip_pattern
from core import rows_of


def supports(H, seed=1):
    M = len(H)
    core = core_patterns(H)
    out = {'core': core}
    # identity-type assignment: row x -> label x (x >= 1), row 0 -> label 1
    assign = {x: x for x in range(1, M)}
    assign[0] = 1
    tw = list(core)
    for x in range(M):
        tw.append(flip_pattern(M, core[assign[x] - 1], x))
    out['twin-id'] = tw
    rng = random.Random(seed)
    labels = list(range(1, M)) + [rng.randrange(1, M)]
    rng.shuffle(labels)
    tw2 = list(core)
    for x in range(M):
        tw2.append(flip_pattern(M, core[labels[x] - 1], x))
    out['twin-rnd'] = tw2
    return out


if __name__ == '__main__':
    name = sys.argv[1]; rho0 = Fraction(sys.argv[2]); rho1 = Fraction(sys.argv[3])
    H = {'S4': lambda: sylvester(2), 'S8': lambda: sylvester(3), 'P8': lambda: paley(7),
         'P12': lambda: paley(11)}[name]()
    M = len(H)
    for kind, T in supports(H).items():
        assert len(set(T)) == len(T), kind
        r = evaluate_support(M, T, rho0, rho1, max_group=4000)
        print(name, kind, '|T|', len(T), 'group', r.get('group'), 'value', r.get('Omega_exact', r.get('value')),
              'lower', r.get('lower_exact'), 'weights', [str(w) for w in r.get('weights', [])][:20],
              r.get('reason', ''), flush=True)
