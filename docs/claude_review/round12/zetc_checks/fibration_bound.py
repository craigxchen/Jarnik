"""Rigorous upper bound for rho of Syl(16, L) via a three-column superset and the fibration (Newton) bound.

Columns a, b, a+b of the Sylvester profile: the 15 nonzero points split by (a.x, b.x) into
  R_ab  = {(1,1)}: in S_a and S_b        (4 points)
  R_a   = {(1,0)}: in S_a and S_{a+b}    (4 points)
  R_b   = {(0,1)}: in S_b and S_{a+b}    (4 points)
  R_0   = {(0,0)}: in no column           (3 points, free)
so Syl(16,L) is contained in  A' = Block x {-1,0,1}^3,
  Block = {k in {-1,0,1}^12 : |s_ab + s_a| <= L, |s_ab + s_b| <= L, |s_a + s_b| <= L}   (s_* = sum over the region).
rho(Syl) <= rho(A') = rho(Block) + 6.

Fibration bound (proved in the write-up, §5.3): Block is the disjoint union over beta = (x,y,z) in
B = {|x|,|y|,|z| <= 4, |x+y|,|x+z|,|y+z| <= L} of products of slices Lambda_x x Lambda_y x Lambda_z, whose functions
have degree <= delta(beta) = 12 - |x| - |y| - |z| (Theorem 2.1 proof).  For any ordering of B with Newton degrees
e_i (minimal degree of a polynomial on Z^3 vanishing at the earlier points and not at beta_i),
rho(Block) <= max_i (e_i + delta(beta_i)).  Newton degrees are computed EXACTLY (fmpz rank profiles).
"""
import itertools, sys
import flint

L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
B = [(x, y, z) for x in range(-4, 5) for y in range(-4, 5) for z in range(-4, 5)
     if abs(x + y) <= L and abs(x + z) <= L and abs(y + z) <= L]
delta = {b: 12 - abs(b[0]) - abs(b[1]) - abs(b[2]) for b in B}


def newton_degrees(order, emax=26):
    """e_i for the given ordering: min e such that row i is a pivot row of the evaluation matrix of
    monomials of degree <= e (rows in the given order)."""
    e_of = [None] * len(order)
    for e in range(0, emax + 1):
        mons = [(i, j, k) for i in range(e + 1) for j in range(e + 1 - i) for k in range(e + 1 - i - j)]
        # rows = points (in order); row rank profile via rref of the transpose
        M = flint.fmpz_mat(len(mons), len(order), [b[0] ** a[0] * b[1] ** a[1] * b[2] ** a[2] for a in mons for b in order])
        res_ = M.rref(); R = res_[0]; rank = res_[-1]
        piv = set()
        r = 0
        for col in range(len(order)):
            if r < rank and R[r, col] != 0:
                piv.add(col)
                r += 1
        for i in piv:
            if e_of[i] is None:
                e_of[i] = e
        if all(v is not None for v in e_of):
            break
    return e_of


# ordering: by delta descending (|beta|_1 ascending), ties broken lexicographically
order = sorted(B, key=lambda b: (-delta[b], b))
e = newton_degrees(order)
bound = max(ei + delta[b] for ei, b in zip(e, order))
print('L=%d |B|=%d  max(e_i + delta_i) = %d  => rho(Block) <= %d, rho(Syl(16,%d)) <= %d' % (L, len(B), bound, bound, L, bound + 6))
worst = sorted(((ei + delta[b], b, ei) for ei, b in zip(e, order)), reverse=True)[:5]
print('worst terms (value, beta, e_i):', worst)
