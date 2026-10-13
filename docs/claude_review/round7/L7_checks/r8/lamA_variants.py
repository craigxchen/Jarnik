# lambda_V(beta) (Lemma 5.3 count, mod p, general-position test curves) for other choices of the six
# event points among the eight cube points of the star arrangement (evidence; cf. lamA_frame.py).
import sys, random, itertools
import flint
import elamp, lamA_frame as L
from star import point_of, quadric_through, LAB
from arr import PLANES
cube = [(1 + i, 3 + j, 5 + k) for i in range(2) for j in range(2) for k in range(2)]
def setup(chosen):
    pts = [point_of([PLANES[i], PLANES[j], PLANES[k]]) for (i, j, k) in chosen]
    Qs = quadric_through(pts)
    assert len(Qs) == 4
    L.qb = [[[L.C(elamp.md(Q[i, j])) for j in range(4)] for i in range(4)] for Q in Qs]
def run(chosen, seed):
    setup(chosen); random.seed(seed)
    lm = L.lam_of(elamp.good_mob(()))
    tot = sum(L.lam_of(elamp.good_mob(T)) for T in itertools.combinations(range(7), 3))
    return lm, tot - 32 * lm
if __name__ == '__main__':
    for omit in [((2, 3, 6), (2, 4, 6)), ((1, 3, 5), (2, 4, 5)), ((1, 3, 5), (2, 4, 6)), ((1, 4, 6), (2, 3, 5))]:
        chosen = [c for c in cube if c not in omit]
        name = 'omit ' + ' & '.join(''.join(LAB[a] for a in o) for o in omit)
        try:
            print(name, [run(chosen, s) for s in (1, 2)], flush=True)
        except Exception as e:
            print(name, 'failed:', repr(e)[:120], flush=True)
