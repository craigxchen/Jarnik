# The exact Gaussian Fourier lattice of a Walsh completion

Walsh completion has integral Fourier coefficients, exact local Gaussian
congruences, and an explicit orthogonal lattice basis. Its minimum-norm
vectors can be classified completely. For a fixed core this recovers only
the canonical completion, up to a common unit and Walsh sign twists.
It supplies no coefficient-height gain or implication that unobserved
completed rows have small phases.

The distinction matters for the full-fair representation: the known rows
have labels `0,e_1,...,e_t` in an ambient group of size `2^t`. They are
sparse, not the positive-density set required in
[the Walsh restriction theorem](walsh_dense_restriction_uniform_bound.md).

## 1. Integral Fourier coefficients and exact energy

Let `G=F_2^t`, `N=2^t`, and write `chi_v(x)=(-1)^(v dot x)`.
For every nonzero `v`, take an odd Gaussian block
`Gamma_v=a_v+i b_v`, and put `p_v=Norm(Gamma_v)` and
`A=product_v p_v`. The full-fair application has conjugate-primitive
nonunits with disjoint rational prime supports; the operator identities
below require only odd nonzero norms. Keep a common Gaussian factor `d`
explicit and assume `d!=0`. The stripped completion is

```text
Z_x=d product_(v!=0) (a_v+i chi_v(x)b_v),
|Z_x|^2=Norm(d) A=:R_0^2.                              (1)
```

Here "stripped" means that independent row units and correcting factors
have not been inserted into the completed rows. Let `H` be the Walsh
matrix, with `H^2=N Id`, and define

```text
C_y=(1/N) sum_x chi_y(x) Z_x.
```

These coefficients are Gaussian integers, despite the division by `N`:

```text
C_y=d sum_(E subset G\{0}, xor(E)=y)
       i^|E| product_(v in E)b_v product_(v notin E)a_v. (2)
```

Expansion and Walsh orthogonality prove (2). Parseval and the constant
modulus of (1) give

```text
sum_y |C_y|^2=R_0^2,
sum_y C_(y+h) bar(C_y)=0       for h!=0.                 (3)
```

The partition sum is `sum_x Z_x=N C_0`. For the known sparse labels,

```text
sum_y C_y=Z_0,
sum_(y_i=1) C_y=(Z_0-Z_(e_i))/2.                       (4)
```

Thus small source chords constrain selected sums of the integral Fourier
coefficients. They do not bound the individual coefficients in those sums.

## 2. Explicit orthogonal operators

Let `P_v` translate a coefficient vector by `v`, and set

```text
T_v=a_v Id+i b_v P_v,       T=product_(v!=0) T_v.
```

The translations commute and are self-adjoint involutions. Consequently

```text
T_v^* T_v=p_v Id,
T^* T=A Id,
C=d T delta_0.                                         (5)
```

In physical coordinates `H T H/N` is diagonal, with entries `Z_x/d`.
Every `T_v` splits into `N/2` two-by-two blocks of determinant `p_v`;
therefore `det T=A^(N/2)`, an odd integer. The Gaussian lattice

```text
L_d=d T Z[i]^N                                         (6)
```

has the explicit orthogonal basis `d T delta_h`, `h in G`, each of
squared length `R_0^2`. Its index in `Z[i]^N` is `(Norm(d) A)^N`.
The minimum calculation uses the orthogonal basis, not merely that index.

## 3. Local congruences and the full product criterion

First suppress `d` and consider one block. Its image is exactly

```text
im(T_v)={B in Z[i]^N:
 Gamma_v | B_y+B_(y+v),
 bar(Gamma_v) | B_y-B_(y+v)       for all y}.             (7)
```

Necessity follows by adding and subtracting the two coordinates of
`T_v U`. For sufficiency, the inverse numerator obeys

```text
2(a_v B_y-i b_v B_(y+v))
 =bar(Gamma_v)(B_y+B_(y+v))
  +Gamma_v(B_y-B_(y+v)).
```

The right side is divisible by `p_v`. Since `p_v` is odd, so is the
undoubled numerator. Hence `T_v^* B/p_v` is Gaussian-integral, proving
(7). Conjugate primitivity is not needed for this single-block assertion.

If the ordinary norms `p_v` are pairwise coprime, then

```text
T Z[i]^N = intersection_v im(T_v).                     (8)
```

Indeed, for a vector in every image, write `B=T_v U_v`. Commutativity
gives `T^* B=p_v (product_(h!=v)T_h^*)U_v`, so each `p_v` divides
every coordinate of `T^* B`. Pairwise coprimality makes their product
`A` divide it, and `T^{-1}B=T^*B/A` is integral. This proves exact
saturation of the local congruence system in the disjoint-core setting.
With common content, require `B in d Z[i]^N` and apply (7)--(8) to `B/d`.

There is also an exact criterion allowing overlapping odd blocks. Define
the physical core rows

```text
A_x=product_(v!=0)(a_v+i chi_v(x)b_v).
```

Then, without any prime-disjointness assumption,

```text
B in T Z[i]^N  iff  A_x | (H B)_x for every x.           (9)
```

Necessity follows by diagonalizing `T`. For the converse the row
divisibilities make `U=H[(H B)_x/A_x]_x/N` an element of
`(1/N) Z[i]^N`. The identity `U=T^{-1}B=T^*B/A` also makes it an
element of `(1/A) Z[i]^N`. Since `N` is a power of two and `A` is odd,
Bezout gives `U in Z[i]^N`. For arbitrary common `d`, again first
require `B/d` Gaussian-integral and apply (9) to it. Physical-row
divisibilities alone must not silently discard the ramified part of `d`.

The shorter local system (8) does require its stated disjointness.
For example, take `G=F_2^2`, `Gamma_1=Gamma_3=2+i`, and
`Gamma_2=3+2i`. Put `p=5` and `B=p T_2 delta_0`. This satisfies every
single-block congruence (7) and has squared norm `A=5^2*13`, but

```text
T^{-1}B=(4 delta_0-2i delta_1-delta_2-2i delta_3)/5
```

is not integral. Thus intersecting the single-block images can create
spurious vectors even on the purported minimum shell. The full criterion
(9) rejects this example. The shared prime is compatible with actual nested
thresholds on the sparse source labels `0,e_1,e_2`: its allocations
are `2,0,1`, with threshold cuts `{0,e_2}` and `{0}`. Thus this warning
is not limited to arbitrary formal blocks. No endpoint realization
is asserted for the example.

## 4. The entire minimum shell and its exact scope

For any nonzero `B=d T U in L_d`, (5) gives

```text
||B||_2^2=R_0^2 ||U||_2^2 >= R_0^2.
```

Equality holds precisely when `U=epsilon delta_h`, with
`epsilon in {1,-1,i,-i}`. Since the operators commute with translations,
the `4N` minimum vectors are exactly

```text
B=epsilon P_h C.
```

Their physical rows are

```text
(H B)_x=epsilon chi_h(x) Z_x.                          (10)
```

The values at `0,e_1,...,e_t` already distinguish these choices: the
anchor fixes `epsilon`, and the `t` other labels fix every bit of `h`.
Thus exact source values give fixed-core uniqueness. If one instead
imposes small angular windows, (10) simply tests the original source
rows under these recorded sign and unit choices. It introduces no
new auxiliary near-phase rows.

In particular, Parseval plus the exact congruence lattice is a complete
reformulation of this fixed-core minimum problem. Every choice of blocks
already has its canonical vector (2) on that shell. This classification
does not establish that a moving core has no shell vector meeting the
source phase windows, nor does it prove a universal negative statement
about every possible use of higher Fourier moments. An arithmetic gain
would have to use how the actual block parameters place the known source
values in those windows, beyond the identities (2)--(10).

## 5. Independent row units and small corrections

For a full row-unit function `epsilon_x in mu_4`, the Fourier operator is

```text
U_epsilon=(1/N) H diag(epsilon_x) H.
```

It is unitary and commutes with every `T_v`. It preserves `Z[i]^N`
(equivalently, preserves `L_d`) **if and only if**

```text
epsilon_x=epsilon_0 chi_h(x)                           (11)
```

for one unit and one character. To prove necessity, a Gaussian-integral
unitary matrix has a single unit entry in each column. The columns of
`U_epsilon` are translates of its first column, so it must be a unit
times one translation. Walsh diagonalization gives (11). Sufficiency
is immediate. Independent source units therefore cannot be silently
included in the same integral Fourier lattice.

For example, changing just one physical row unit from `1` to `i`
changes each coefficient by `(i-1)Z_x chi_y(x)/N`. For odd core rows
and `d=1`, this is not Gaussian-integral when `N>=4`. Common ramified
content can hide that denominator in a particular vector, but does
not change the lattice-preservation criterion (11).

Now suppose the actual known source rows have

```text
z_i=epsilon_i K_i Z_i,       Norm(K_i)=K^2,
|z_i|=R=K R_0.                                         (12)
```

The equal correction norms follow from the common source and core norms.
The exact form of the small-chord information in (4) is

```text
sum_(y_i=1) C_y=(z_0/(epsilon_0 K_0)-z_i/(epsilon_i K_i))/2,

|sum_(y_i=1) C_y
  -(z_0/2)(1/(epsilon_0 K_0)-1/(epsilon_i K_i))|
 =|z_0-z_i|/(2K) <= C_arc sqrt(R)/(2K).                 (13)
```

This retains the actual error and the correction phase center. When all
`epsilon_i K_i` agree, the center vanishes; generally it can have size
comparable to `R_0`, even if the correcting heights are small compared
with every core height.

A prescribed full equal-norm correction function `K_x` acts by
`H diag(K_x) H/N`, a scaled unitary operator. It need not be
Gaussian-integral: at order four, take `K_x=2+i` except for one
value `2-i`. Its Fourier coefficients have a half-integral change.
Nor is its extension from the sparse known labels canonical. Bounds
on `log|K_i|` alone do not determine those missing phases or restore
the lattice (6). For the uniform problem, selecting a common source
unit class is legitimate, but the actual fair-core corrections in
(12)--(13) must still be retained.

## Verification

[check_walsh_completion_fourier_ideal_lattice.py](check_walsh_completion_fourier_ideal_lattice.py)
checks integral Fourier expansion, energy and autocorrelation, explicit
orthogonal bases, the local and physical-row membership tests, the overlap
counterexample, all order-four row-unit patterns, and correction phases.
The literal Gaussian examples are not claimed to lie on endpoint arcs.
