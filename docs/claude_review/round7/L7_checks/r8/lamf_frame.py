# Validation of lamA_frame.py's method on the forgetful foliation F_i (leaves: x_i moves, others fixed).
# Lemma 4.1 gives lambda_{f_i} = -psi_i + sum_j D_ij exactly; lamf.lam_f evaluates it on class vectors.
# Here the same quantity is computed by the local-frame method on Mobius test curves (mod p).
import sys, random, itertools
import flint
from elamp import P, good_mob
from lamA_frame import lam_from_br, C
from lamf import lam_f, split_key, LAB
n = 7
def brackets(mob, i):
    X = flint.nmod_poly([0, 1], P)
    N = [C(m[1]) + C(m[0]) * X for m in mob]      # x_a(t) = (al t + be : ga t + de)
    D = [C(m[3]) + C(m[2]) * X for m in mob]
    dN = [C(0)] * n; dD = [C(0)] * n
    dN[i] = C(1); dD[i] = C(0)                     # leaf direction: x_i -> (N_i + u : D_i)
    br0 = {}; br1 = {}
    for a in range(n):
        for b in range(n):
            if a != b:
                br0[(a, b)] = N[a] * D[b] - N[b] * D[a]
                br1[(a, b)] = dN[a] * D[b] + N[a] * dD[b] - dN[b] * D[a] - N[b] * dD[a]
    return br0, br1
def cls_of(T):
    c = {split_key(pr): 2 for pr in itertools.combinations(range(n), 2)}
    if T:
        for pr in itertools.combinations(T, 2): c[split_key(pr)] = 1
        c[split_key(T)] = 1
    return c
if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    bad = 0
    for i in range(n):
        tot = 0; lm = None
        for T in [()] + list(itertools.combinations(range(n), 3)):
            mob = good_mob(T)
            got = lam_from_br(*brackets(mob, i))
            exp = lam_f(i, cls_of(T))
            if got != exp: bad += 1; print('MISMATCH i', i, 'T', T, got, exp)
            if T: tot += got
            else: lm = got
        print('i =', LAB[i], ' lambda_{f_i}(mu) =', lm, ' frame lambda_{f_i}(beta) =', tot - 32 * lm, flush=True)
    print('mismatches:', bad)
