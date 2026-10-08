# An endpoint bound in a cancellation-free block-phase cone

There is an elementary uniform endpoint theorem for the cone defined below.
It allows arbitrary moving sign matrices and Gaussian block sizes. Its
substantive hypothesis is that, after orientation at one actual anchor,
the block phases have one sign and do not wind. This hypothesis is not
known to follow from a general endpoint arc.

The complementary exact realization lemma shows why replacing the endpoint
scale by the assertion that the angular width tends to zero cannot supply
that missing reduction. Every finite sign matrix admits primitive actual
circle realizations with angular width tending to zero and any prescribed
positive rational limiting block weights. In particular the Paley
four-row obstruction survives actual simultaneous angular collapse.

## 1. The cone and its uniform theorem

Suppose the source points have the literal factorization

```text
z_i = d u product_(j=1)^r Gamma_j^((1+s_ij)/2)
                            bar(Gamma_j)^((1-s_ij)/2),
R = |d| product_j |Gamma_j|.
```

Here `d` is a nonzero Gaussian integer, `u` is one common Gaussian unit,
and the `Gamma_j` are nonunit Gaussian integers coprime to their own
conjugates. Every column is nonconstant. A whole cut class can be grouped
into one block; `r` counts only the nonconstant nonunit blocks actually
used. Constant blocks belong in `d`, and unit blocks are removed while
retaining the resulting row units. Thus the common-unit hypothesis must
hold **after** this removal. Distinct blocks need not have disjoint prime
support for this theorem.

Choose an actual anchor row `0` and conjugate individual blocks, with
the corresponding column reversal, so that `s_0j=+1` for every `j`.
Write `A_i={j:s_ij=-1}`. Assume that the blocks have literal argument
lifts `phi_j=arg(Gamma_j)` satisfying

```text
phi_j >= 0,
sum_(j in A_i) phi_j < pi/2       for every row i.       (1)
```

One may simultaneously reverse all phase signs. Condition (1) includes
the phase branches: moving each block independently to a convenient
sector by a Gaussian unit can change the row units and is not permitted
without checking the displayed literal factorization again. The simpler
condition `phi_j>=0, sum_j phi_j<pi/2` suffices.

If r>=1 and the points lie in any arc of angular width `Delta`, then

```text
Delta >= 2 (|d|/R)^(1/r).                              (2)
```

Indeed, their exact relative arguments are

```text
arg(z_i/z_0) = -2 sum_(j in A_i) phi_j in (-pi,0].      (3)
```

The geodesic distance of this row from the anchor is therefore the
absolute value of the right side of (3), namely
`2 sum_(j in A_i) phi_j`, and is at most `Delta`.
Every block occurs in at least one `A_i`, so `2phi_j<=Delta` for all j.
A conjugate-primitive nonunit cannot be real. Its nonzero integral
imaginary coordinate consequently gives

```text
1 <= |Im Gamma_j| = |Gamma_j| sin(phi_j)
  <= |Gamma_j| phi_j <= |Gamma_j| Delta/2.
```

Multiplying over the actual blocks proves (2). This is a simultaneous
phase budget: its use of the nonnegative separating sums in (3) is
unavailable from separate scalar character heights alone.

At endpoint scale `Delta<=C/sqrt(R)`, (2) gives, for `r>=3`,

```text
R <= [(C/2)^r/|d|]^(2/(r-2)).                         (4)
```

If `C<=2`, no such `r>=3` configuration exists: (4) would give `R<=1`,
whereas the product contains three nonunits. If `C>2`, uniformly over
all `r>=3`,

```text
R <= (C/2)^6.                                         (5)
```

The cases `r=0,1,2` are separate. With a common row unit they give at
most `1,2,4` distinct points, respectively. Therefore a valid completely
explicit count for this cone is

```text
M(C) = 4                                      if C<=2,
M(C) = max(4, 4 floor((C/2)^6)+2)               if C>2. (6)
```

For (6), a radius-R integer circle has at most two points at each
integer x-coordinate, hence at most `4floor(R)+2` points. The constants
are deliberately coarse; the result is a uniform count for the stated
cone, not a new general growth estimate. Arbitrary independent row units
are outside this theorem. The existing
[finite-sector certificate](finite_sector_certificate.md) instead uses
a specific signed four-block character and allows cancellation; (2)
uses an anchor and every nonconstant block, with no specified sign matrix.

## 2. Every sign profile has actual primitive angular-collapse realizations

Let `S` be any finite sign matrix with distinct rows and nonconstant
columns, and prescribe positive integers `d_1,...,d_r`. Choose distinct
positive even integers `a_1,...,a_r`, and a positive integer `L` divisible
by every `|a_j^2-a_k^2|` for `j!=k` (take `L=1` when `r=1`). For
positive multiples T of L put

```text
H_j(T)=a_j T+i,                  Gamma_j(T)=H_j(T)^d_j,
z_i(T)=product_j Gamma_j(T)^((1+s_ij)/2)
                    bar(Gamma_j(T))^((1-s_ij)/2),
E=sum_j d_j,                    R_T=product_j |H_j(T)|^d_j.
```

These are distinct, Gaussian-primitive integer-circle points with the
**exact** intended block profile. Here is the complete arithmetic check.

* Each `H_j` is coprime to its conjugate: a common Gaussian prime would
  divide `2i`, but `N(H_j)=a_j^2T^2+1` is odd.
* The rational norms are pairwise coprime. A prime dividing the norms at
  j and k divides
  `a_k^2 N(H_j)-a_j^2 N(H_k)=a_k^2-a_j^2`. Such a prime divides T,
  making both norms congruent to one modulo that prime, a contradiction.
* At a prime belonging to block j, some row chooses each orientation
  because that column is nonconstant. Thus the whole tuple has Gaussian
  gcd one. If two rows differ, a prime in a differing block distinguishes
  their valuations. They therefore give distinct actual points.

There is no primality assumption and no unproved assertion about primes
in narrow sectors. Individual blocks may have several split prime
factors, but every one belongs to precisely the prescribed column.

Their lifted arguments and norm weights satisfy

```text
theta_i(T)=sum_j s_ij d_j atan(1/(a_j T)),
Delta_T <= 2 sum_j d_j/(a_j T),
log N(Gamma_j(T)) = 2d_j log T + 2d_j log a_j + O(T^-2),
W_T=log R_T^2 = 2E log T + O(1).
```

Consequently `Delta_T -> 0`, while the normalized block weights tend to
`d_j/E`. Positive rational weight vectors are all obtained, and are dense
in the positive weight simplex. This disproves any additional closed
condition on limiting normalized cut weights inferred solely from actual
primitivity and `Delta_T -> 0`, unless that condition already holds for
every positive cut-weight vector in the chosen sign matrix.

If S has an all-plus row, this realization lies in the cone (1) for all
large T. A different row has a nonempty negative set A, so

```text
theta_0(T)-theta_i(T)
  = (2/T) sum_(j in A) d_j/a_j + O(T^-3).
```

It follows that `Delta_T` is asymptotic to a positive constant times
`T^-1`, and

```text
Delta_T sqrt(R_T) is asymptotic to a positive constant
times T^(E/2-1).                                      (7)
```

Thus for E>2 these actual angular-collapse families escape every fixed
endpoint constant. Equation (2) detects the escape without taking a
limit or keeping S fixed.

## 3. An actual Paley counterexample to a scale-free extraction rule

Use the normalized order `M=q+1` Paley matrix from
[the quartic-defect note](general_four_row_quartic_slack_obstruction.md),
where `q=3 mod 4` is prime, and omit its constant column. Take
`r=q`, `d_j=1` in Section 2. Denote its column signs by `s_ij` and put

```text
w_j(T)=log N(H_j(T)),
G_xy(T)=sum_j w_j(T) s_xj s_yj,
T_Q(T)=sum_j w_j(T) product_(x in Q) s_xj,
W_T=sum_j w_j(T).
```

All these are the actual prime-layer weights and correlations of the
primitive integer-circle tuple, since different block norms are coprime.
Hadamard orthogonality gives, for x!=y,

```text
G_xy(T)=-2log T+O_q(1).
```

In particular every pair-Gram necessary inequality `G_xy<=-log 4` for
the formal choice `C=1,D=0` holds for large T. That does not assert the
corresponding actual endpoint condition.

For a quadruple Q, let `t_Q=sum_j product_(x in Q) s_xj`. The character
estimate already proved and sourced in the quartic-defect note gives
`|t_Q|<=3sqrt(q)+4`. Therefore, uniformly over **all** quadruples at
this fixed q,

```text
(W_T-T_Q(T))/W_T -> 1-t_Q/q
                  >= 1-3/sqrt(q)-4/q.                  (8)
```

The right side tends to one with q. Choosing T increasingly large at
each q produces actual primitive tuples with `M->infinity`,
`Delta_T->0`, strict pair obtuseness, and quartic defect
`(1-o(1))W_T` at every quadruple. Every conic compatibility identity and
every actual phase identity is automatically satisfied. Row deletion
and nonnegative averaging over quadruples cannot remove this defect.

This is an exact counterexample to the candidate implication

```text
actual primitive circle points + Delta->0 + strict pair obtuseness
  => some quadruple has quartic defect o(W).
```

It is not a counterexample with bounded `Delta sqrt(R)`. In fact (7)
has exponent `q/2-1>0`. No `W=O(Mlog M)` claim is made. The earlier
prime-box note achieves that radius cost without controlling angles;
the present construction controls the actual simultaneous angles and
retains every intended cut, but its endpoint constant diverges. The
fixed-Hadamard Roth theorem already excludes endpoint realizations of
this unmutated profile; Section 1 provides a direct, uniform elementary
exclusion in its no-cancellation cone.

## Verification and remaining requirement

Run [the companion checker](check_one_sided_block_phase_cone.py). It
checks Gaussian norms, all pairwise block-norm gcds, distinctness and the
whole tuple gcd exactly; verifies the Paley correlations and every
quadruple at four orders; and checks the cone/product inequalities on
several literal Gaussian fixtures. Its floating-point phase evaluations
are diagnostics only. The proofs above establish the infinite claims.

A forward general reduction would need to control oppositely signed
anchored block phases or their winding while retaining the factor
`sqrt(R)`. Compactness of unscaled source angles discards precisely the
information used in (2)--(5). Neither such control nor a general growth
improvement is established here.

Root and a separate Astra instance independently audited the proof,
including the literal unit and phase-branch assumptions and every
primitivity and limiting-weight claim.
