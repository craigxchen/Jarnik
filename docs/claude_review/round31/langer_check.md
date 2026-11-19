# Langer's Theorem 0.1 checked in the original (coordinator, round 31)

The Part IV referees and judge could confirm the cited theorem only through secondary quotations, because arxiv.org
was not reachable. After the network allowlist was extended the coordinator read the arXiv version of the paper:
A. Langer, *Logarithmic orbifold Euler numbers of surfaces with applications*, arXiv:math/0012180v1 (19 Dec 2000),
published in Proc. London Math. Soc. 86 (2003) 358-396. The PDF is not stored in the repository.

**Theorem 0.1 (page 2), as printed.** "Let X be a normal projective surface with a Q-divisor D = sum a_i D_i,
0 <= a_i <= 1. Assume that the pair (X, D) is log canonical and a multiple of K_X + D is effective. Then
3 e_orb(X, D) >= (K_X + D)^2. Moreover, if equality holds then K_X + D is nef."

This is the form used in (IV.2): coefficients in [0,1], log canonical, K_X + D Q-effective (not nef).

**Global orbifold Euler number (page 2; Definition 3.4).**
e_orb(X, D) = e_top(X) - sum_i a_i e_top(D_i - Sing(X, D)) + sum_{x in Sing(X, D)} (e_orb(x; X, D) - 1).

**Local value at a normal crossing (Theorem 8.3, three lines L_1, L_2, L_3 through 0 in C^2 with weights
a_1 <= a_2 <= a_3).** e_orb(0; C^2, D) = (1 - a_1 - a_2)(1 - a_3) when a_3 >= a_1 + a_2. With a_1 = 0 this gives
(1 - a_2)(1 - a_3) at a node of two components with arbitrary weights in [0,1], and with a_1 = a_2 = 0 it gives
1 - a_3 at a smooth point of one component.

**Consequence.** At a node of D_i and D_j the term (1 - a_i)(1 - a_j) - 1, together with removing the node from both
curves, contributes + a_i a_j relative to e_top(X) - sum a_i e_top(D_i). So for a simple-normal-crossing weighted
boundary e_orb = e_top(X) - sum_i a_i e_top(D_i) + sum_{nodes} a_i a_j, which is the formula used in Part IV and by
the referees. The judge's item "Not checked" in `README_IV.md` section 1.1 and the caveat "conditional on the
citation" are thereby resolved for the arXiv version; the published version was not compared with it.
