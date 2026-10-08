# An integer-lattice minimum for endpoint sign patterns

This note treats exact Gaussian sign-pattern rows at fixed `M`, with
arbitrary moving split primes, nested powers, complete common Gaussian
content, and the row-unit ranges of
[isolable_pattern_bessel_phase_gap.md](isolable_pattern_bessel_phase_gap.md).
The integer row-difference lattice cannot contain a short nonzero
vector at unbounded endpoint radius. The proof combines the existing
whole-pattern Bessel mass bound, a fixed rational phase target, and
the all-order multiplicative independence from
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md).
The weighted rank-one theorem in
[weighted_rank_one_power_class_packing.md](weighted_rank_one_power_class_packing.md)
packs many two-coordinate certificates; here one short integral
certificate suffices, in any corank. This is a local structural
restriction, not a uniform point-count bound.

Use the exact threshold-layer factorization from the isolable-pattern note:

```text
z_i=d epsilon_i product_(j=1)^T Gamma_j^((1+s_ij)/2)
                              bar(Gamma_j)^((1-s_ij)/2),
W=log R^2,       V_j=log Norm(Gamma_j)<=2W/M.             (1)
```

All `Gamma_j` are conjugate-primitive nonunits. The columns `s_j` represent
distinct unoriented nonconstant cuts. Define the rational row-difference
space and its integer lattice by

```text
V_Q={lambda^t S:lambda in Q^M, sum_i lambda_i=0},
L_S=V_Q intersect Z^T.                                  (2)
```

**Integer pattern-lattice theorem.** Fix `M>=9` and `C` in either
endpoint range: one literal row-unit class with `0<C<=2`, or arbitrary
row units with `0<C<=sqrt(2)`. There is a finite, generally ineffective
`B(M,C)` such that every exact endpoint cluster with `R>B(M,C)` obeys

```text
8 ||v||_1 >= M       for every nonzero v in L_S.        (3)
```

The cutoff is strict: a vector with `8||v||_1<M` bounds the radius.
The constant is uniform over changing rational primes, nested widths,
common Gaussian content, and all possible sign-class matrices at fixed
`M`; it is not uniform as `M` grows. A standard coordinate `e_j` has
norm one, recovering the isolable-pattern exclusion for `M>=9`.

For the corank-one application, suppose the row-difference matrix has
rank `T-1`. On an actual endpoint arc, affine allocation rigidity makes
this rank `M-1`, because every allocation column is a sum of its sign
threshold columns modulo a constant. Hence `T=M`. If the primitive
integer kernel generator `r` had a zero coordinate, its standard
coordinate `e_j` would lie in `L_S` and violate (3) for `M>=9` at
large radius. Thus all kernel coordinates are nonzero there.
Conjugate `Gamma_j` and reverse `s_j` when `r_j<0`; the source rows
are unchanged. Write `a_j=|r_j|>0`. For each pair
`j!=h`, the integral vector

```text
g=gcd(a_j,a_h),
v_jh=(a_h/g)e_j-(a_j/g)e_h in L_S,
||v_jh||_1=(a_j+a_h)/g.                                  (4)
```

Thus every pair at unbounded radius satisfies
`8(a_j+a_h)/gcd(a_j,a_h)>=M`, regardless of whether its cuts are
comparable. In particular, at `M>=17` all positive primitive kernel
weights must be pairwise distinct, forcing
`max_j a_j>=M` and `sum_j a_j>=M(M+1)/2`. These are necessary profile
conditions, not sufficient endpoint constructions.

## 1. Every integral relation gives an exact Gaussian direction

Take a nonzero `v in L_S`. Clear denominators and parity to choose an
integer row vector `lambda`, with `sum_i lambda_i=0`, and an integer
`t>=1` such that

```text
lambda^t S=2t v.                                      (5)
```

The signed row product has the *literal* identity

```text
Q=product_i z_i^lambda_i
 =u product_j (Gamma_j/bar(Gamma_j))^(t v_j),
u=product_i epsilon_i^lambda_i in mu_4.               (6)
```

The complete common factor cancels because `sum lambda_i=0`. At every
moving prime, all threshold layers contribute according to their exact
class in (5), including nested layers at the same prime. No
independent-prime approximation is being made.

Form the raw Gaussian integer

```text
B_raw=product_(v_j>0) Gamma_j^v_j
      product_(v_j<0) bar(Gamma_j)^(-v_j).              (7)
```

Divide out its positive rational split-prime content to obtain the
conjugate-primitive `Beta_v`. Then
`B_raw/bar(B_raw)=Beta_v/bar(Beta_v)` and (6) becomes
`Q=u(Beta_v/bar(Beta_v))^t`. Exact valuations and (1) give

```text
V(Beta_v)=log Norm(Beta_v)
 <=sum_j |v_j| V_j
 <=2 ||v||_1 W/M.                                      (8)
```

## 2. Endpoint affine rigidity makes the reduction nonunit

The all-order independence statement in
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
applies under the same two arc/unit ranges for `R>1`: a nonzero
integer `lambda` with `sum lambda_i=0` cannot make
`product_i z_i^lambda_i` a Gaussian unit. Here `lambda` is nonzero
because (5) has nonzero right side. To see the exact mechanism, if
`Q` were a unit, its valuation at each chosen split Gaussian prime
would give `sum_i lambda_i a_i(p)=0`; these equalities and
`sum lambda_i=0` would be a forbidden affine relation among the
allocation rows. The row-unit factor and complete common gcd have no
effect on those valuations.

If `Beta_v` were a unit, (6) would make `Q` a unit. Therefore every
nonzero `v in L_S` reduces to a conjugate-primitive nonunit on an
actual endpoint arc. This step is stronger than a purely combinatorial
test for two cuts. In corank one, proportional-power collision of any
two blocks is impossible on the actual endpoint arc whenever their
pair vector (4) lies in `L_S`; nested cuts do not create an exception.

## 3. Phase height and the lattice cutoff

Lift point arguments in the containing arc of angular width
`Delta<=C/sqrt(R)`. Equation (5) and `sum lambda_i=0` give

```text
|sum_i lambda_i theta_i|<=||lambda||_1 Delta/2.
```

Writing `phi=arg(Beta_v)`, equation (6) isolates it to the fixed
rational-`pi` grid

```text
T_(lambda,u)={(2pi m-arg u)/(2t):m in Z} mod 2pi,
dist(phi,T_(lambda,u))
    <=||lambda||_1 Delta/(4t).                         (9)
```

For the nonunit `Beta_v`, the fixed-target Gaussian phase lemma and
Roth give, for every `epsilon>0`,

```text
V(Beta_v)>=W/(4+2epsilon)-K_(S,v,C,epsilon).         (10)
```

When `8||v||_1<M`, choose `epsilon>0` with
`1/(4+2epsilon)>2||v||_1/M`. Equations (8) and (10) bound `W`, hence
`R`. For fixed `M`, only finitely many sign-class matrices and integral
vectors with `8||v||_1<M` occur; each has one fixed integer certificate
(5), and the row-unit shifts are finite. Taking the maximum of their
bounds makes `B(M,C)` support-uniform.

The theorem does not claim that Bessel obtuseness by itself produces a
short integral vector. The earlier weighted transportation criterion
can still exclude profiles whose integer lattice passes this local
minimum test.

## 4. A full-support corank-one witness to the norm-only gap

Take `M=T=64`. Index rows by `x in F_2^6`. Use the 63 nonconstant
Walsh columns `chi_v(x)=(-1)^(v dot x)`, and one extra column

```text
f(x)=(-1)^(x_1 x_2+x_3 x_4+x_5 x_6)
     with its sign flipped at x=0.                        (11)
```

The unflipped bent function has every Walsh coefficient `+8` or `-8`.
The flip at zero subtracts two, so every nonconstant coefficient of
`f` is `+6` or `-10`. The primitive row-difference kernel therefore has
old-column weights `3` or `5` and extra-column weight `32`, all
nonzero. Orient the old columns to make these weights positive.
At least two old columns share a weight; all distinct balanced Walsh
cuts are incomparable. Their pair vector has `||v_jh||_1=2<64/8`, so the theorem
rules out a large-radius endpoint realization of this exact profile.

Yet with abstract layer weights `1` on every old Walsh column and
`1/2` on the extra column, the row Gram has diagonal `63+1/2` and
every off-diagonal entry `-1+f(x)f(y)/2<0`. Thus the norm/obtuse test
alone accepts this sign matrix. The checker also realizes it with
literal distinct split primes, arbitrary row units, and a nonunit
common factor; its pair norms satisfy the endpoint-derived `C=1/2`
inequality. Their actual arguments are uncontrolled. The local phase
comparison, rather than the norm test, supplies the obstruction.

## 5. The full fair matrix passes this test by an exponential margin

For the matrix containing all `2^(M-1)-1` nonconstant unoriented cuts,
the lattice (2) has the exact minimum

```text
min_(0!=v in L_S) ||v||_1 = 2^(M-2).                   (12)
```

Indeed, put `c_i=2lambda_i`. Since `sum lambda_i=0`, the value at a
cut `T` is `v(T)=sum_(i in T)c_i`. The singleton cut classes force
every `c_i` to be an integer, and `sum c_i=0`. Conversely every such
integer vector gives a member of `L_S`. Reversing a representative cut
only changes the sign of that coordinate.

Choose a nonzero coordinate `c_j`. Pair the subsets not containing `j`
with those containing it. The triangle inequality gives

```text
sum_(T subset [M]) |sum_(i in T)c_i|
 =sum_(U subset [M]\{j}) (|sum_U c_i|+|sum_U c_i+c_j|)
 >=2^(M-1)|c_j| >=2^(M-1).
```

The empty and full subset values vanish. Complementary subsets have
opposite values, so dividing by two gives the lower bound in (12).
Equality holds for `c=e_i-e_j`, since precisely half of all subsets
separate these two labels. Thus even the shortest full-fair character
is much longer than the forbidden threshold `M/8`. The lattice theorem
does not exclude the general extracted fair profile.

Run [check_corank_one_local_bessel_pair_filter.py](check_corank_one_local_bessel_pair_filter.py)
for exact Walsh/kernel, comparability, certificate, Gaussian source,
pair-norm identities, and finite checks of (12). Roth and the Bessel inequality remain proof
inputs, not finite-check conclusions.
