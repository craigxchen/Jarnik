# notes (parked derivations), round 19 search

## Excess identity (derived from the referee formula, ref_gain.py)
Notation: a_T = #edges of A inside T (all visible T, pairs included), c_T likewise, mu_T = min (0 for pairs,
since A,C disjoint), B = C - lam A, b_e = mult of pair e in B, interior Q-irreducible Y with mult m_Y,
Y.beta = 31 deg Y - sum_{|T|>=3} nu_T(Y). S = {T, |T|>=3 : a_T = c_T, nu_T(B) > mu_T}.
N0 = #{T visible : a_T != c_T}.
Code: Lb = 31D - sum mu_T; Fb = #{T in Ts: max(a,b,c) > mu} + #{e in A u C u pairB} + sum_Y Y.beta.
Then (degree count 31 D = 31 sum b_e + 31 sum m_Y deg Y, sum_T nu_T(B) = 30 sum b_e + sum_Y m_Y sum_T nu_T(Y),
nu_T(B) = mu_T off S):
   gain = -N0 + sum_e (b_e - 1) + sum_{T in S} (nu_T(B) - mu_T - 1) + sum_Y (m_Y - 1) Y.beta.     (*)
= -N0 + sum over components of F_1 of (mult - 1)(beta-degree).

## lam != +-1
(1) no pair factor of B_lam (restriction to x_i=x_j: products of differences proportional => ratio +-1)
(2) S empty (initial forms along L_T are products of differences in y and restricted x => ratio +-1)
=> gain = -N0 + sum_Y (m_Y-1) Y.beta; = -N0 if B_lam squarefree.
(3) B_lam non-squarefree => Disc_t(c - lam a) = 0 on every line (binary-form discriminant, formal degree D,
vanishes also if the restriction is identically 0). So Lambda(A,C) \ {+-1} subset rational roots of Delta_l
for ANY one line l with Delta_l not identically 0. Intersection over two lines still a superset.

## lam = +-1, fast path
T in S needs in_T(C) = s in_T(A); A,C disjoint => inside-T edge multisets both empty => mu_T = 0 on S.
Squarefree B (tested on a line: B|l squarefree of degree D => B squarefree) => gain = -N0 + sum_{T in S}(nu_T(B)-1).
