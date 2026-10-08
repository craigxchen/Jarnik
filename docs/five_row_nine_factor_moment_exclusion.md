# Five equal first-and-third moments require at least ten factors

**Theorem (including repeated factors).** Let `alpha_1,...,alpha_k` be
distinct positive real numbers. If five distinct integer rows `x_i`
have both sums

```text
sum_j x_ij alpha_j,       sum_j x_ij alpha_j^3
```

independent of the row `i`, then their total coordinate width satisfies

```text
sum_j (max_i x_ij-min_i x_ij) >= 10.
```

In particular, five distinct binary rows with equal first and third
moments for nonzero real coefficients of distinct absolute values require
at least ten columns. With repeated coefficients, rows must first be
identified when they give the same aggregate allocation at each absolute
value; exchanging identical copies does not produce distinct rows for
the theorem.

Sections 1--3 prove the distinct-absolute-value case by extending the
real `K5` exclusion to one zero edge. Section 4 removes the distinctness
restriction on factor copies using an additional exact polynomial
certificate. The previously certified
[ten-factor real witness](tenfactor_real_outercut_witness.md) shows that
the threshold is attained over the reals, with a different code.

This is an exact moment theorem. No replacement of arbitrary nearby
Gaussian factors by these fixed coefficient moments is asserted, and no
general lattice-circle growth bound follows from it.

## 1. Pairwise distance and the parity extension

First take the binary case with nonzero coefficients `a_j` of distinct
absolute values and `n` columns.
Subtract two row equations. On their differing coordinates this gives
nonzero numbers `c_j=+a_j` or `-a_j` of distinct absolute values with
`sum c_j=sum c_j^3=0`. There must be at least five terms. For one or two
terms this is immediate. For three terms Newton's identity gives
`3c_1 c_2 c_3=sum c_j^3=0`, a contradiction. For four terms the first
and third elementary symmetric functions vanish, so their monic root
polynomial is even. Its nonzero roots occur in opposite pairs, again
contradicting the distinct absolute values.

Thus the five binary rows have pairwise Hamming distance at least five.
Summing the ten pair distances, each coordinate contributes `k(5-k)<=6`,
where `k` is its number of ones. Hence

```text
50 <= sum_(i<j) d_H(b_i,b_j) <= 6n.                   (1)
```

This excludes every `n<=8`.

For `n=9`, append the parity of each row as its tenth bit. Every new pair
distance is even and is at least six. The lower and upper bounds for the
total distance are now both 60. Consequently every pair distance is six
and every column splits the five rows two against three.

Complement any column whose smaller side consists of zeros, so each
column is the incidence vector of an unordered pair of rows. Replace its
coefficient by its negative when complementing: this changes every first
moment by the same constant and likewise every third moment, preserving
the two equalities. Give the appended coordinate coefficient zero; it
does not change either moment equation.

Let `n_ij` count columns whose two-row side is `{i,j}`, and let
`d_i=sum_(j!=i) n_ij`. There are ten columns, so `sum_i d_i=20`, and

```text
d_i+d_j-2n_ij=6.
```

Summing over the four `j!=i` gives `d_i+20=24`, hence `d_i=4` and
`n_ij=1`. The ten columns are exactly the ten edges of `K5`, with one
edge coefficient zero and the other nine having distinct nonzero
absolute values. It remains to exclude precisely this boundary system.

## 2. The zero edge still passes only the old four order cones

Label the zero edge `01`, place it first in increasing absolute-value
order, and assign it an arbitrary formal sign `+`. This is a convention
for the order test, not a definition of `sign(0)`.

A difference of two vertex equations uses six edges. If it omits `01`,
the six nonzero terms have one of the four successor sign patterns used
in the [real K5 proof](k5_real_first_third_moment_exclusion.md). If it
includes `01`, its five nonzero terms have, up to negation, the critical
sign pattern `++--+` by the
[critical odd-moment theorem](critical_odd_moment_root_order.md).
The formal first sign of the zero term is its assigned sign times its
incidence difference, so it can be either sign. All four possibilities
are covered by the same successor patterns:

| Formal zero sign | Five-term sign pattern | Six-term order word |
| --- | --- | --- |
| `+` | `++--+` | `+++--+` |
| `+` | `--++-` | `+--++-` |
| `-` | `++--+` | `-++--+` |
| `-` | `--++-` | `---++-` |

Therefore every putative boundary solution survives the same necessary
finite order test. Its 48 words still lie in the same four orbits. The
enumeration is a necessary superset; no assertion that every word has
feasible magnitudes is needed.

## 3. The exact certificates stay positive at the boundary

Write the ordered absolute values as
`0=r_0<r_1<...<r_9`, with gaps `g_0=0` and
`g_i=r_i-r_(i-1)>0` for `i>=1`. The first-moment kernel and its exact
six-parameter parametrizations in the K5 proof do not change. In each
orbit, one parameter is literally `g_0`; the other five are literal
positive gap coordinates.

Let `P(lambda)` be the integer polynomial combination of the four
cube-sum equations used in that orbit. Restrict it by setting the zero-gap
parameter to zero. The existing certificates give:

| Orbit | Parameter set to zero | Surviving monomials | Least coefficient |
| ---: | --- | ---: | ---: |
| 0 | `lambda_0` | 31 | 3 |
| 1 | `lambda_1` | 48 | 108 |
| 2 | `lambda_2` | 33 | 3 |
| 3 | `lambda_2` | 45 | 3870 |

Every listed coefficient is positive. Since all other parameters are
positive, each restricted polynomial is strictly positive. It is still
a polynomial combination of the four equations that would vanish in a
moment solution. This contradiction excludes the zero-edge K5 system,
finishes the `n=9` case, and proves the theorem.

The [standard-library checker](check_five_row_nine_factor_moment_exclusion.py)
reconstructs all four certificates from the original edge incidence and
kernel columns, checks the complete order enumeration, and performs the
zero-parameter restriction exactly. It also checks the parity-extension
example and the unique linear cut-incidence system. Root supplied the
parity reduction and boundary argument; Luna and Sol independently
reconstructed the boundary certificates and audited the proof.

## 4. Repeated coefficients and arbitrary integer allocations

Translate each coordinate of the five integer rows so that its minimum
is zero, and put `m_j=max_i x_ij`. This changes both moments by constants
independent of the row. Remove coordinates with `m_j=0`. It suffices to
exclude `n=sum_j m_j<=9`.

For a difference of two distinct aggregate rows, replace an integer
coefficient `c_j=x_ij-x_lj` by `|c_j|` copies of the signed number
`sign(c_j) alpha_j`. There are no opposite pairs: each magnitude occurs
with only one sign. The same Newton identities as in Section 1 show
that a nonempty list of at most four such numbers cannot have both its
first and third power sums zero. Thus every aggregate pair has
`sum_j |x_ij-x_lj|>=5`, even when its terms have repeated magnitudes.

Embed each coordinate in its `m_j` threshold bits
`b_(i,j,r)=1` if `x_ij>=r`, for `1<=r<=m_j`, giving each bit coefficient
`alpha_j`. This preserves the two moments and makes Hamming distance
equal to aggregate integer distance. The distance and parity argument
of Section 1 applies unchanged. For `n<=8` it gives a contradiction;
for `n=9` it gives exactly the ten `K5` edge columns, with only the
appended parity coefficient zero.

Within each coefficient group the threshold supports are nested.
Every support has size two or three. Two thresholds cannot have the
same support, because the ten minority edge columns are all distinct.
Consequently each group has at most two thresholds. If there are two,
their supports are a two-set `A` and a three-set `B`, with `A` contained
in `B`. Minority orientation turns them into the disjoint edges
`A` and the complement of `B`, with coefficients `+alpha_j` and
`-alpha_j` respectively.

The [opposite disjoint edge lemma](repeated_k5_zero_edge_exclusion.md)
proves that a `K5` solution with a zero edge and this opposite nonzero
pair must have another zero edge. Here all nine original coefficients
are nonzero, so this is impossible. If there is no repeated group, the
nine coefficients have distinct absolute values and Sections 2--3 give
the contradiction. This proves the full integer-allocation theorem.

The added lemma has two symmetry cases. Its exact integer identities
express `3200 u` and `9600` as polynomial combinations of the four
cube-sum differences, where `u` is a second edge in the first case.
The [separate checker](check_repeated_k5_zero_edge.py) verifies both
identities and the complete linear parametrizations. Astra supplied
the certificates; root and Sol independently audited the reduction,
and Sol independently reconstructed the two algebraic cases.

## 5. Limited consequence for fixed factor constructions

Consider fixed rational coefficients, now allowing repetitions, opposite
values, zero values and inactive columns, in the signed construction

```text
w_i(T)=(-1)^|S_i| product_(j in S_i)(a_j+iT)
                    product_(j notin S_i)(a_j-iT).
```

Let `T` tend to positive infinity through integers, and use one fixed
denominator clearing for the rational coefficients.

Write the two oriented factors as `l_a=-(a+iT)` and `r_a=a-iT`.
Changing `a` to `-a` exchanges `l_a,r_a` exactly. At `a=0` they are
identical. Therefore remove zero factors and group all nonzero factors
by their positive absolute value `alpha_j`. Let `x_ij` count the
`l_(alpha_j)` factors in row `i`, after these exchanges, and identify
rows with the same aggregate vector.

The exact common polynomial factor at group `j` is
`l_(alpha_j)^(min_i x_ij) r_(alpha_j)^(m_j-max_i x_ij)`, where `m_j`
is the number of copies in that group. Removing these factors leaves
a primitive tuple over `Q(i)[T]` of degree

```text
n_eff=sum_j(max_i x_ij-min_i x_ij).
```

There are no other common polynomial factors: the distinct positive
`alpha_j` give distinct Gaussian linear factors. Polynomial Bezout
applied to the whole reduced tuple, with fixed denominators cleared,
gives polynomial coefficients whose linear combination is a fixed
nonzero Gaussian integer. Evaluating at integral `T` shows that the
remaining common Gaussian integer gcd divides that fixed integer.
Thus its modulus is bounded and the primitive radius is of order
`T^n_eff`. Cancelling the common polynomial factor changes all phases
equally and preserves their angular span.

The successive row-dependent angular terms are proportional to the
aggregate first and third moments at orders `T^(-1)` and `T^(-3)`.
For `7<=n_eff<=9`, bounded normalized endpoint span would force both
moments to agree, since `n_eff/2>3`. The theorem excludes five distinct
aggregate rows in this effective-width range.

This is a consequence for fixed rational coefficients after exact common
factor removal. Raw factor count cannot replace `n_eff`. The separate
[repeated six-factor classification](sixfactor_repeated_magnitude_classification.md)
settles effective width six at `C<=1/2`; together with the elementary
smaller-width cases it excludes five persistent rows of every effective
width at most nine. Moving coefficients and arbitrary source tuples
still require separate arguments. There is no general growth improvement
or uniform circle-point bound here.
