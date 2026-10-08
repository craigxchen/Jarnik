# Independent check of the triangle-cover Runge budget

The conditional theorem in
[six_point_triangle_cover_runge_budget.md](six_point_triangle_cover_runge_budget.md)
checks out. This audit verifies its stated hypotheses and comparisons; it
does not assert that every full fair configuration meets those hypotheses.

- At any point in `B_i`, the parity of the valuation of `F^w` is exactly
  `w_i`. Thus the Kummer classes are independent even over Qbar, and an
  m-dimensional subspace gives geometric degree `q=2^m`. The local
  geometric squareclass group at a boundary point has rank one: its
  nontrivial class is the uniformizer. Thus each active `B_i` has inertia
  two, with no extra geometric inertia from unit factors. At `ell=0`
  all valuations are even, so there is no geometric ramification there.
- Riemann--Hurwitz gives
  `2g(X)-2=8q+(4a)(q/2)=q(8+2a)`, hence `g(X)=1+q(4+a)`.
  The pole support has `(40-4a)q+4a(q/2)=q(40-2a)` points. These
  count distinct poles; active poles have higher orders after pullback.
- A finite Galois field k splitting every upstairs pole exists and also
  splits its downstairs image. For the specialized radicals,
  `L=k L0` is elementary abelian over k with `h=[L:k]<=q`.
  Splitting all geometric poles makes the L-pole orbit count equal to
  the geometric count even when h is smaller than q.
- At a good prime dividing only `M_i(s,t)`, primitive parameters and
  the excluded resultants ensure that ell, the other minors, and Delta
  are units. The point P consequently has good reduction on C into B_i.
  For a radical which is a unit there, its square root evaluated at a
  k-rational upstairs pole is a nonzero element of k. Reduction of P
  to that boundary point and the exclusions make the specialized unit
  congruent to this square. Odd-residue-characteristic Hensel lifting
  makes it a square in `k_v`.
- If i is inactive, all radicals are local squares, so v has h places
  above it in L. If i is active, choose one basis vector evaluating to
  one in coordinate i; the other basis vectors are local squares. The
  local degree is therefore at most two, giving at least h/2 places.
  Odd depth and unramifiedness of k at p force ramified quadratic local
  degree, residue degree one, and exactly h/2 places. Even depth need
  not force ramification, but does not invalidate the lower bound.
- When B_i is one degree-four Q-orbit, conjugation in `Gal(k/Q)` carries
  a reduction of rational P at one place to each of its four boundary
  points. These must occur at four distinct places: otherwise P would
  reduce to two distinct good boundary sections at the same place.
  Thus an active partition contributes at least 2h finite places and
  an inactive one at least 4h. Distinct rational primes make these
  contributions disjoint, yielding `h(40-2a)`.
- Odd valuation rows spanning `W*` detect every nonzero product of
  specialized radicals, including over k because p is unramified there.
  They imply h=q. This is sufficient, not necessary. Only when h=q
  does the stated finite-place lower bound equal the full geometric
  pole count; archimedean places then preclude Runge's strict surplus.

A fixed nonzero rescaling `c psi` changes none of these good-prime
valuations after the primes supporting c are included in the exception
set. That set also must retain the model, root-separation, boundary-unit,
and field-ramification exclusions. None is a cost independent of the
moving frame unless proved separately. The argument does not count
`q/2` places directly over Q_p; residual unit extensions there can enlarge
the local degree to four.

The principal limits are correctly retained: degree-two exceptional pole
orbits, specialization-degree drops, missing good private primes, even
valuation rows of insufficient rank, and point-dependent twists are not
excluded. Logarithmic fair prime mass supplies none of the required
odd-depth or outside-exception-set assertions on its own.

The separate
[constant-two-cover archimedean obstruction](six_point_constant_two_cover_archimedean_obstruction.md)
already states the matching full-fair caution. It proves `G>=A_12` (indeed
S_12) only for one specified first-ruling witness, and conditions its
cover obstruction on that group bound and the ten biquadratic pole
orbits. It makes no claim that all fair weights satisfy either condition.
The Q-defined multiplication-by-two cover with nonconstant group scheme
remains outside both arguments' claimed scope.
