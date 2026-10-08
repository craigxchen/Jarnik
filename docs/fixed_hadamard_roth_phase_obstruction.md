# Roth separation for fixed Hadamard block profiles

For each fixed normalized Hadamard profile, five nonunit Gaussian blocks
already prevent an unbounded family on arcs of length `C sqrt(R)`.
Independent Gaussian units on the rows are allowed. The argument is an
ineffective bounded-radius theorem, not a uniform endpoint theorem.

The Walsh case without independent row units is substantially covered by
existing, stronger elementary results: Sections 6.1--6.2 of
[the research notes](codex_uniform_bound_research.md) exclude every
common-unit block-simplex using seven rows and discriminant-code
collisions. The [order-eight theorem](semibent_eight_point_phase_obstruction.md)
also handles all normalized order-eight Hadamard matrices at `C=1/2`,
including their nonbalanced columns and unit blocks. The addition here is
the direct Roth argument with arbitrary row units and the fixed-matrix
extension in Section 5. No novelty claim beyond this repository comparison
is intended.

## 1. The fixed Hadamard theorem

Fix `C>0` and an `m x m` Hadamard matrix whose first column is all ones.
Delete that column to obtain `H`, an `m x (m-1)` sign matrix. Thus

```text
H^T H = m I,                  H^T 1 = 0.                (1)
```

Choose nonzero Gaussian integer blocks `kappa_v`, each satisfying
`gcd(kappa_v,bar(kappa_v))=1` up to units. Unit blocks are allowed. Let
`d` be any nonzero Gaussian integer, and let `u_i` be arbitrary Gaussian
units, independently chosen for each row. Put

```text
z_i = d u_i product_v kappa_v^((1+H_iv)/2)
                         bar(kappa_v)^((1-H_iv)/2),
n_v = |kappa_v|^2,
R   = |d| product_v |kappa_v|,
product_v n_v = R^2/|d|^2 <= R^2.                       (2)
```

**Theorem.** There is a finite constant `B(H,C)` such that, if at least
five blocks are nonunits and all the `z_i` lie on an arc of length at
most `C sqrt(R)`, then `R<=B(H,C)`.

Distinctness and pairwise disjoint prime support are not needed for this
theorem. They may be imposed to recover the intended primitive independent
block model. The constant can be chosen uniformly over the finitely many
Hadamard matrices of a fixed order `m`; it is generally ineffective.

In particular, take `m=2^r>=8`, rows `i` in `F_2^r`, columns
`v` in `F_2^r \ {0}`, and `H_iv=(-1)^(i dot v)`. If all `m-1` blocks
are nonunits, the proposed full Walsh profile has bounded radius for every
fixed `m,C`, even with independent row units.

## 2. Phase inversion retains all torsion

Write `theta_v=arg(kappa_v)` using any real representatives and
`u_i=exp(pi i t_i/2)`, with `t_i` integral. A containing arc has angular
width at most

```text
Delta = C R^(-1/2).
```

Choose midpoint lifts of its point arguments. After absorbing `arg(d)`
into a common real `c`, there are integers `l_i` and errors `e_i` with

```text
H theta + (pi/2)t = c 1 + 2 pi l + e,
|e_i| <= Delta/2.                                       (3)
```

Multiplication by `H^T/m` gives the exact equality

```text
theta_v = (2 pi/m) sum_i H_iv l_i
          - (pi/(2m)) sum_i H_iv t_i
          + (1/m) sum_i H_iv e_i.                        (4)
```

Consequently every block argument is within `Delta/2`, modulo `2 pi`,
of a member of the fixed finite set

```text
T_m = { pi j/(2m) mod 2 pi : 0 <= j < 4m }.              (5)
```

These targets depend on the row lifts and units, but belong to one finite
set independent of `R` and the blocks. The unknown center `c` disappears
because the columns are balanced. The integer lifts and row-unit term in
(4) are essential; dropping either can incorrectly shrink the target set.

## 3. Gaussian arguments near a fixed algebraic direction

The input is Roth's theorem: for a fixed real algebraic irrational
`alpha` and every `epsilon>0`, there exists `c(alpha,epsilon)>0` such that

```text
|alpha-p/q| >= c(alpha,epsilon) q^(-2-epsilon)            (6)
```

for all reduced rational `p/q`, `q>=1`. Its usual finiteness formulation
implies this all-denominator version by decreasing the positive constant
to cover the finite exceptional set. See K. F. Roth,
[Rational approximations to algebraic numbers, Mathematika 2 (1955), 1--20](https://www.cambridge.org/core/journals/mathematika/article/rational-approximations-to-algebraic-numbers/EFB89E2873019B246F004FAA06E05A7F),
Theorem on printed page 2. This is a fixed-target application.

**Phase lemma.** Fix a finite set `T` of rational multiples of `pi`.
For each `epsilon>0` there is `c(T,epsilon)>0` such that every
conjugate-coprime nonunit `kappa=a+bi` satisfies

```text
dist(arg(kappa),T) >= c(T,epsilon) |kappa|^(-2-epsilon).
                                                               (7)
```

Distances are taken modulo `2 pi`. To prove this, choose a nearest
`beta in T`. If the angular distance is at least `1/8`, (7) is immediate
after reducing the constant. Otherwise choose a slope chart as follows:

* If `|cos(beta)|>=1/sqrt(2)`, use `x=b/a`, `alpha=tan(beta)`.
* Otherwise use `x=a/b`, `alpha=cot(beta)`.

Along the short interval joining the two arguments the chosen sine or
cosine has absolute value greater than `1/2`. The denominator is nonzero,
and the derivative of the chosen chart has absolute value at most four.
Thus, if the angular distance is `eta`,

```text
|x-alpha| <= 4 eta.                                     (8)
```

Both charts cover vertical directions without a limiting argument. The
reduced denominator `q` of `x` satisfies `q<=|kappa|`. The target `alpha`
is algebraic because `exp(2 i beta)` is a root of unity. If it is
irrational, (6)--(8) give (7).

If `alpha=A/B` is rational in lowest terms, a nonzero difference satisfies
`|x-alpha|>=1/(Bq)`, giving the stronger estimate

```text
eta >= c(T) |kappa|^(-1).                               (9)
```

The difference cannot vanish. Indeed, exact equality would make
`kappa/bar(kappa)=exp(2 i beta)` a root of unity in `Q(i)`. Such an element
is an algebraic integer in `Z[i]` of modulus one, hence is one of
`1,-1,i,-i`. Therefore `kappa` lies on a coordinate axis or a diagonal.
A conjugate-coprime Gaussian integer on an axis is a unit; on a diagonal
it is divisible, together with its conjugate, by `1+i`. Both contradict
the hypotheses. This also explains why unit blocks must be excluded from
the phase lemma. Taking minima over the finite target set completes its
proof.

## 4. Multiplying five block lower bounds

Apply the phase lemma to (4)--(5). For each nonunit block,

```text
(C/2) R^(-1/2) >= c(m,epsilon) n_v^(-1-epsilon/2),
n_v >= a(m,C,epsilon) R^(1/(2+epsilon)),                 (10)
```

where `a>0` is fixed. This is the exponent for the **norm** `n_v`,
not for the modulus `|kappa_v|`. If its target has rational slope in the
chosen chart, (9) instead gives `n_v >= a'(m,C) R`.

Select any five nonunit blocks and multiply (10). All unselected block
norms are at least one, so (2) gives

```text
R^2 >= a^5 R^(5/(2+epsilon)).                            (11)
```

Choose any `0<epsilon<1/2`; the exponent on the right exceeds two.
For example `epsilon=1/4` gives `R^(2/9)<=a^(-5)`, hence a finite
bound `R<=max(1,a^(-45/2))`. This proves the theorem. The general
`s`-nonunit comparison is `s/(2+epsilon)>2`; five is the first integer
for which Roth alone supplies a strict exponent gap. No assertion about
the existence of four-block endpoint families follows from equality at
the limiting exponent.

## 5. Fixed matrices beyond Hadamard

Let `S` be any fixed `m x k` sign matrix such that the `k+1` columns
of `[1 S]` are linearly independent over `Q`. There is then a fixed
rational matrix `L` satisfying

```text
L S=I_k,                       L 1=0.                   (12)
```

Use the block model (2) with `S` replacing `H`. Applying `L` to the lifted
phase equation gives

```text
theta = 2 pi L l - (pi/2)L t + L e.                     (13)
```

If `D` is a common denominator of the entries of `L`, the targets belong
to `{pi j/(2D) mod 2 pi}`, and the error in coordinate `v` is at most
`(Delta/2) sum_i |L_vi|`. The same finite-target lemma and product
comparison show that five nonunit blocks force `R<=B(S,C)`.

Full column independence is sufficient, not necessary: the same proof
works whenever five chosen nonunit coordinates have rational isolating
rows `L` satisfying `L1=0` and `LS` equal to the corresponding five
coordinate vectors. The total product budget still includes every block.

Precisely, fix an anchor row `0`, and let `V` be the rational row span of
the vectors `S_i-S_0`, `i!=0`. Call coordinate `v` isolable if its
standard coordinate vector `e_v` belongs to `V`. This is equivalent to
the existence of a rational row vector `lambda_v` with
`lambda_v 1=0` and `lambda_v S=e_v`: write a linear combination of row
differences and give the anchor coefficient the negative sum of all the
other coefficients. The converse follows by reversing that operation.

**Isolable-coordinate theorem.** For a fixed sign matrix `S` and `C>0`,
there is `B(S,C)<infinity` such that every model (2) on an arc of length
at most `C sqrt(R)` with `R>B(S,C)` has at most four isolable nonunit
coordinates. Choose isolating rows once for all isolable coordinates,
apply (13) coordinatewise, and minimize the finitely many positive
constants in (10). Five such coordinates would contradict (11).

There is no pairwise support restriction in this statement. Blocks may
share Gaussian primes, including blocks arising from nested threshold
layers, provided each individual block is conjugate-coprime and the
displayed sign factorization and radius product are valid. Eliminating
the constant row direction is essential: membership of `e_v` merely in
the row span of `S` would leave the unknown center `c` in (13), and does
not justify a fixed algebraic target.

## 6. Scope

The matrix and its rational inverse data are fixed before the radius
grows. Roth's constants depend on the resulting algebraic targets.
No effective or uniform bound as `m` grows is obtained. More fundamentally,
a large family of balanced cut columns need not have even five rationally
isolable coordinates: there may be many more columns than row differences.
That prevents using (4) or (13) on general conductor profiles.

Independent row units are harmless here because they only enlarge a fixed
finite target set. Independent moving Gaussian twists would produce moving
targets, so the proof does not cover them. Conjugate-coprimality of each block and
the actual shrinking phase condition are both indispensable assumptions.
The domain-capacity examples therefore remain useful incidence obstructions;
this result excludes the specified fixed full Walsh realizations, without
closing the uniform lattice-circle endpoint problem.

The exponent comparison and the phase, torsion, rational-target, vertical,
and unit-block cases have been checked directly above. This is a prose
theorem using Roth; it has not been formalized in Lean.

An independent Sol audit and the root proof review checked the row-unit
terms, cancellation of the common phase, both slope charts, the exact-target
exclusion, the norm exponent, and the fixed-matrix extension. Exact rational
checks also verify the Walsh inverse identities through order 64 and the
exponent difference `5/(2+1/4)-2=2/9`; these do not replace Roth's theorem.
