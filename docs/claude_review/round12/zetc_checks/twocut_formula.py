"""Two-cut blocks: compare nu_A(u) (mod p) with the 'atom' upper bound
F(u) = min over choices of atoms (block subspaces R^X with cost nu(R^X), column directions 1_{S_j} with cost 2L_j,
coordinates with cost 2) whose sum contains u, of the total cost.
If nu = F for all tested u, the sublevel sets of nu agree (on the tests) with unions of sums of atoms, which are
rational subspaces (evidence for ZE-TC^rat on two-cut blocks).
"""
import itertools, random, sys
from fractions import Fraction
from tclib import box_points, nu_mod, cube_monomials


def block(a, b, c, L1, L2):
    d = a + b + c
    P = list(range(0, a)); Q = list(range(a, a + b)); R = list(range(a + b, d))
    S1 = P + Q; S2 = Q + R
    pts = box_points(d, [(S1, L1), (S2, L2)])
    return d, P, Q, R, S1, S2, pts


def subsets_nonempty(regs):
    out = []
    for r in range(1, len(regs) + 1):
        for comb in itertools.combinations(regs, r):
            X = sorted(set().union(*[set(x) for x in comb]))
            out.append(X)
    return out


def generic_point(X, d, rng):
    u = [0] * d
    for x in X:
        u[x] = rng.randrange(10**6, 10**7)
    return u


def atom_costs(d, P, Q, R, S1, S2, pts, mons, rng):
    regs = [r for r in (P, Q, R) if r]
    blocks = {}
    for X in subsets_nonempty(regs):
        nu, _ = nu_mod(pts, generic_point(X, d, rng), mons=mons)
        blocks[tuple(X)] = nu
    return blocks


def F_bound(u, d, S1, S2, L1, L2, blocks):
    best = 2 * sum(1 for x in u if x != 0)
    blist = [()] + list(blocks.keys())
    # choose up to 2 blocks (their union of coordinates is free) -- blocks are unions of regions
    choices = [((), 0)]
    bl = list(blocks.items())
    for i in range(len(bl)):
        choices.append((bl[i][0], bl[i][1]))
        for j in range(i + 1, len(bl)):
            choices.append((tuple(sorted(set(bl[i][0]) | set(bl[j][0]))), bl[i][1] + bl[j][1]))
    for free, cost in choices:
        freeset = set(free)
        rest = [x for x in range(d) if x not in freeset]
        for use1 in (0, 1):
            for use2 in (0, 1):
                ccost = cost + use1 * 2 * L1 + use2 * 2 * L2
                if ccost >= best:
                    continue
                # minimise support of u - a 1_S1 - b 1_S2 on rest
                cand_a = {Fraction(0)}
                cand_b = {Fraction(0)}
                if use1:
                    cand_a |= {Fraction(u[x]) for x in rest if x in S1}
                if use2:
                    cand_b |= {Fraction(u[x]) for x in rest if x in S2}
                if use1 and use2:
                    # also a = u_x - b for x in S1&S2 etc: enumerate b, then a from residuals
                    pass
                m = None
                for bb in cand_b if use2 else [Fraction(0)]:
                    resid_a = {Fraction(0)}
                    if use1:
                        resid_a |= {Fraction(u[x]) - (bb if x in S2 else 0) for x in rest if x in S1}
                    for aa in resid_a:
                        s = sum(1 for x in rest if Fraction(u[x]) - (aa if x in S1 else 0) - (bb if x in S2 else 0) != 0)
                        m = s if m is None else min(m, s)
                if use1 and use2:
                    # also enumerate a first, then b from residuals
                    for aa in {Fraction(0)} | {Fraction(u[x]) for x in rest if x in S1}:
                        for bb in {Fraction(0)} | {Fraction(u[x]) - (aa if x in S1 else 0) for x in rest if x in S2}:
                            s = sum(1 for x in rest if Fraction(u[x]) - (aa if x in S1 else 0) - (bb if x in S2 else 0) != 0)
                            m = min(m, s)
                best = min(best, ccost + 2 * m)
    return best


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    ntests = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    tot = bad = 0
    for d in range(3, dmax + 1):
        mons = cube_monomials(d)
        for a in range(0, d + 1):
            for b in range(1, d + 1 - a):
                c = d - a - b
                if a < c:
                    continue
                for L1 in range(1, a + b):
                    for L2 in range(1, b + c):
                        dd, P, Q, R, S1, S2, pts = block(a, b, c, L1, L2)
                        blocks = atom_costs(d, P, Q, R, S1, S2, pts, mons, rng)
                        nb = 0
                        for t in range(ntests):
                            vals = rng.choice([[-1, 0, 1], [0, 1, 2], [-2, -1, 0, 1, 2], [0, 1, 2, 3, 5, 8]])
                            u = [rng.choice(vals) for _ in range(d)]
                            if all(x == 0 for x in u):
                                continue
                            nu, _ = nu_mod(pts, u, mons=mons)
                            F = F_bound(u, d, S1, S2, L1, L2, blocks)
                            tot += 1
                            if nu != F:
                                bad += 1; nb += 1
                                if nb <= 3:
                                    print('DIFF a=%d b=%d c=%d L=(%d,%d) u=%s nu=%d F=%d' % (a, b, c, L1, L2, u, nu, F), flush=True)
                        print('block a=%d b=%d c=%d L=(%d,%d) |A|=%d blocks=%s diffs=%d' % (a, b, c, L1, L2, len(pts), blocks, nb), flush=True)
    print('total', tot, 'diffs', bad)
