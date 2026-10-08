# A half-mass barrier to Walsh extraction from Paley cut profiles

This note gives a quantitative obstruction to extracting a repeated-Walsh
core from a linear-size, non-Walsh cut system. It concerns literal primitive
Gaussian circle points with actual distinct split-prime weights, but makes
no claim that their arguments lie on a fixed endpoint arc. The profiles
pass every subset determinant inequality required by such an arc and the
pair and balanced four-row phase-height inequalities. Those necessary
conditions therefore do not by themselves imply a Walsh extraction with
sublinear exceptional prime-cut mass.

## Profile and exact determinant margin

Let `q` run through primes congruent to `3 mod 4`, set `M=q+1`, and let
`H` be the normalized Paley Hadamard matrix described in
[the quartic-defect note](general_four_row_quartic_slack_obstruction.md).
Delete its all-ones column. Take a fixed `b>=5` physical copies of each
of the remaining `M-1` columns and flip one physical entry in each row,
using distinct physical columns. There are `r=b(M-1)` columns. Assign
them distinct split rational primes with

```text
tau <= w_j=log p_j < tau+1/r,       tau=4log M+O(1).
```

The prime selection is given in the quartic-defect note. Choose one
Gaussian prime `pi_j` over each `p_j` and set

`z_x=product_(j:s_xj=+1) pi_j product_(j:s_xj=-1) bar(pi_j)`.
Every `z_x` has squared modulus `R^2=product_j p_j`. Every physical
column is nonconstant even after its single flip, so the complete
Gaussian gcd of the rows is one. The rows are distinct because their
sign words are distinct. Their arguments remain uncontrolled. Put
`W=sum_j w_j=log R^2`.

For a selected set `I` of `k>=2` rows define its inherited determinant
defect by

```text
D_I=(1/4)sum_j w_j (sum_(i in I) s_ij)^2.
```

Before the flips, Hadamard orthogonality and deletion of the constant
column give, **for every** `I`,

```text
D_I^0 = b tau k(M-k)/4,
kW^0/4-D_I^0 = b tau k(k-1)/4.                 (1)
```

Only flips assigned to rows in `I` affect `D_I`; there are exactly `k`
of them. A flip at row `i` changes the relevant squared row sum by
`-4 sum_(l in I, l!=i) s_ij s_lj`, hence changes `D_I` by at most
`(k-1)tau`. Consequently, at equal weights,

```text
kW^0/4-D_I >= (b/4-1)tau k(k-1).             (2)
```

The distinct-prime perturbations have total weight less than one. Since
`|k-(sum_(i in I)s_ij)^2|/4 <= k(k-1)/4`, they lower the margin by at
most `k(k-1)/4`. Thus the actual integer-prime profile satisfies

```text
kW/4-D_I >= [(b/4-1)tau-1/4] k(k-1)
           >= [(tau-1)/4] k(k-1).              (3)
```

In particular, for all sufficiently large `M`, (3) is stronger for
every subset than the fixed `C0=1/2` endpoint determinant requirement
`D_I <= kW/4-k(k-1)log 2`. The pair and balanced four-row phase-height
inequalities also hold for sufficiently large `M`, by the estimates in
the quartic-defect note. Any other fixed positive arc constant changes
only their bounded logarithmic thresholds.

## A fixed fraction of mass is exceptional on every Walsh cube

For four distinct retained rows `Q`, let

```text
T_Q=sum_j w_j product_(i in Q)s_ij.
```

The Paley character-sum estimate, the four one-entry flips incident to
`Q`, and the total prime-weight perturbation give uniformly in `Q`

```text
|T_Q| <= b(3sqrt(q)+4)tau+8tau+1 = o(W).       (4)
```

Now select any retained rows and label a subset of them by an affine
cube over `F_2`; it suffices that the labels contain a parallelogram
`x_1+x_2+x_3+x_4=0` on four distinct rows `Q`. Every Walsh character,
with either orientation, has four-sign product `+1` on `Q`.

Suppose physical columns of total weight `W-E` agree on the retained
rows with Walsh characters, allowing the usual one distinct entry flip
per retained row across all good columns. The four rows of `Q` touch at
most four such flips, of weight at most `tau+1/r` each. Every other good
column contributes `+w_j` to `T_Q`, while exceptional columns contribute
at worst `-w_j`. Therefore

```text
T_Q >= W-2E-8(tau+1/r),
E >= (W-T_Q)/2-4(tau+1/r)
  >= W/2-O(b sqrt(M)tau).                        (5)
```

For fixed `b`, (5) says `E/W>=1/2-O(M^(-1/2))`. The same conclusion
holds after arbitrary row deletion, including a restriction to just
four rows intended as a complete Walsh square. Constant restricted
columns have four-sign product `+1`, so retaining them as common
Gaussian content does not evade (5). Column reorientation also leaves
the four-sign product unchanged. The inequality charges actual
logarithmic prime-cut mass and includes all physical columns; it does
not assume equal formal weights.
Even if one permits a fixed sign twist at each retained row, every good
column has the same four-sign product `sigma_Q` on `Q`. Replacing
`T_Q` by `sigma_Q T_Q` in (5) and using `|T_Q|=o(W)` gives the same
half-mass lower bound.

There is also an entrywise version, independent of the number or balance
of proposed Walsh labels. Let the retained rows be a complete affine
cube `A` of size `n>=4`. For each physical column `j`, choose any Walsh
character on `A` and any orientation, and let `B` be the total weighted
number of disagreements over all `n*r` row-column entries:

```text
B=sum_j w_j #{x in A: s_xj differs from the chosen character}.
```

For every affine parallelogram `Q` in `A`, a column with four-sign
product `-1` must disagree at an odd number of its four entries.
Equations (4) and the exact negative-mass identity give

```text
(W-T_Q)/2 <= sum_(x in Q) sum_(j: mismatch at x,j) w_j.
```

Average this over all affine parallelograms in `A`. Every row occurs
with probability `4/n`, while `(W-T_Q)/2 >= W/2-O(b sqrt(M)tau)`
holds for each `Q`. Hence

```text
B >= nW/8-O(n b sqrt(M)tau).                    (6)
```

The relative error term tends to zero for fixed `b`. Thus even a
weighted `o(nW)`-entry approximation by Walsh characters is impossible
on a retained complete cube. The usual one-flip-per-row corrections
have entry cost at most `n(tau+1/r)`, far below (6).
The argument also permits fixed row sign twists: on each parallelogram
the target product is a common sign, and its correlation with the
Paley profile still has absolute value at most `|T_Q|`.

The obstruction is specific to a full Walsh cube or any proposed
approximation whose retained row labels contain a parallelogram. An
affinely independent labelling can avoid all such quadruples, but it is
not a complete repeated-Walsh row family. The argument does not exclude
other structures or simultaneous source-phase constraints. In
particular, the uncontrolled arguments of the
literal Gaussian points prevent interpreting this profile as an
endpoint counterexample.

With a capacity-two flip assignment, the subsequent
[prime-field theorem](paley_primefield_polynomial_character_gap.md)
shows that these same unbounded profiles also satisfy every individual
saturated primitive half-height inequality. Thus passing all those
inequalities and all subset determinant tests still does not justify
the proposed Walsh extraction. The joint source phases remain an
additional hypothesis that has not been established for these profiles.

The [exact checker](check_paley_walsh_mass_extraction_barrier.py) tests
the determinant margin on every nontrivial row subset at `q=11` and on
1,000 sampled subsets at `q=19`. It checks every quadruple at both
orders, including the rational-weight exceptional-mass budget with up
to four good-column flips, and tests the entrywise inequality on an
eight-row labelled cube. The analytic prime clustering and the
uniform Paley character-sum estimate are supplied by the linked
quartic-defect note, not by these finite checks.
