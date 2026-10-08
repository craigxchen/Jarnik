# Pair-norm squareclasses: what the parity code does and does not force

This note isolates the information carried by the squareclasses of the
ordinary pair norms.  It uses the exact relation

`d_ab d_bc = q_abc^2 d_ac`,

where `q_abc` is the positive rational cancellation factor from the
Gaussian allocation.  The angular hypotheses are deliberately kept
separate: the construction below is an arithmetic allocation model and is
not asserted to be a short-arc circle cluster.

## What the full angular hypotheses add

The countermodel below concerns allocation information alone. The full
angular hypotheses do give a Sidon conclusion: the three-partition proof in
[the independent audit](angular_pair_squareclass_independent_audit.md)
excludes repeated pair squareclasses for fixed-unit C<=2 clusters. An older
four-wise baseline already implies the same qualitative packing bound after
subdividing arbitrary fixed-C arcs. Thus this construction identifies which
hypotheses are necessary; it does not obstruct that angular proof.

## The parity graph

Let `rho_a` be the vector of prime-allocation exponents for row `a`, reduced
modulo two, in `F_2^omega`, where `omega` is the number of rational split
primes dividing `N`.  The squareclass of a pair norm is

`sigma_ab = rho_a + rho_b`.

This is just the triangle identity modulo squares.  Since the pair norms
in the short-arc reduction are nonsquares, every `sigma_ab` is nonzero.  In
particular the `rho_a` are distinct, giving the familiar `M <= 2^omega`.

There is one small refinement that is sometimes useful.  If two edges share
an endpoint and have the same label, then

`sigma_ab = sigma_ac  =>  sigma_bc=0`,

which is forbidden.  Thus every nonzero label occurs on a matching of the
complete graph.  Counting edges by labels gives

`M(M-1)/2 <= (2^omega-1) floor(M/2)`.

For even `M` this is exactly `M <= 2^omega`; for odd `M` it gives
`M <= 2^omega-1`.  It does not produce a meaningful asymptotic improvement.

If one could prove that all edge labels were distinct, then the stronger
Sidon/cap estimate

`M(M-1)/2 <= 2^omega-1`

would follow.  But that premise is false as a consequence of ordinary norm
injectivity: repeated labels can occur only on disjoint edges, and their
ordinary norms can still be made different by even exponent data.

## An explicit nearcritical allocation with repeated labels

The following finite primitive construction works for every `r` and hence
for arbitrarily large `M=2^r`.  It shows that parity labels, exact norm
injectivity, and nearcritical pair norms are compatible at the allocation
level.

Take the rows to be the points `x` of `F_2^r`.  Use one base coordinate for
each nonzero linear form `ell` on `F_2^r`; there are `M-1` such coordinates.
Choose distinct split rational primes `p_ell` and odd exponents `s_ell` so
that

`s_ell log(p_ell) = t + O_r(1)`

as `t` tends to infinity.  (Taking the nearest odd integer to
`t/log(p_ell)` suffices.)  Give row `x` exponent `s_ell` at the coordinate
`ell` when `ell(x)=1`, and exponent zero when `ell(x)=0`.  Every base
coordinate has a row with exponent zero, so this allocation is already
primitive.  Let `N_base` be the product of these prime powers.

For a pair `x,y`, the difference `c=x+y` is nonzero and exactly
`2^(r-1)` nonzero forms satisfy `ell(c)=1`.  Consequently

`log Norm(A_xy) = 2^(r-1)t + O_r(1)`,

while `log N_base = (2^r-1)t + O_r(1)`.  Before the even tag is added,
the pair exponent is therefore one half of the base prime count plus the
fixed `1/2*t` excess caused by having `M-1` rather than `M` base forms.
The squareclass label is the vector `(1_{ell(c)=1})_ell`, so it depends
only on `c`; each label is repeated by the disjoint matching of pairs with
the same difference.

To make all ordinary pair norms distinct without changing any squareclass,
add one new split prime `q` and give row `x_k` the even exponent
`2*(2^k-1)` for an arbitrary ordering `x_0,...,x_{M-1}`.  The minimum
exponent is zero, so this tag preserves primitivity.  Its contribution to
the pair norm is

`q^(2 |2^i-2^j|)`.

All these exponents are distinct for unordered pairs: the 2-adic valuation
of `2^j-2^i` recovers `i`, and the remaining odd factor recovers `j-i`.
Thus the tagged pair norms are all distinct.  Their squareclasses are
unchanged because the added exponents are even.  With
`N=N_base*q^(2*(2^(M-1)-1))`, for each fixed `r` and large `t` one has

`log Norm(A_xy) = (M/2)t + O_r(1)` and
`log N = (M-1)t + O_r(1)`.

Equivalently,

`log Norm(A_xy) = (1/2 + 1/(2(M-1))) log N + O_r(1)`.

For `M>=4`, the `t/2` lower-margin and the
`(M-2)t/(2M)` upper-margin dominate the fixed tag and rounding errors, so
all tagged norms lie in
`[sqrt(N), N^(1/2+1/M)]`.  Every pair remains a nonsquare because its base
label is nonzero and the base exponents are odd.

All primes can be chosen congruent to one modulo four, so this is compatible
with the split-prime allocation bookkeeping and gives ordinary divisors of
the resulting odd `N`.  It is an exact integer construction; no assumption
of equal prime weights is being made.

## Scope of the construction

The rows above specify exponents, hence pair norms and their squareclasses,
but they do not control the arguments of the Gaussian primes.  Choosing
Gaussian factors above the split primes gives actual complex pair factors,
yet their arguments need not lie in one short arc.  Therefore the
construction does not challenge the geometric theorem and does not show
that a large actual cluster exists.

It does show that squareclass identities plus distinct ordinary norms alone
cannot justify a Sidon condition or a forbidden four-row additive relation.
Any such improvement must use extra angular information.  The conclusion from
parity information alone here is the matching multiplicity bound,
which is essentially the existing `M <= 2^omega` estimate.

The [exact allocation checker](check_pair_norm_squareclass_code.py) verifies
the M=4 and M=8 constructions using distinct fixed split primes and large
odd exponents. It checks primitivity, all norm inequalities with integer
powers, distinct nonsquare norms, and the repeated perfect-matching labels.
It does not test or assert the endpoint angular condition.
