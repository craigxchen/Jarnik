# Endpoint-map existence reduces to small source parameters and a square test

For a fixed primitive endpoint source, the search for an arbitrary rational
nonconformal endpoint image can be reduced to three small arithmetic
parameters. They are a divisor of the **source** squared radius, a Gaussian
integer, and one positive integer. The potentially large shear entry is
then forced by a linear formula; the remaining diagonal is recovered by an
integer square test. The target radius is computed only after this recovery.

This is a complete finite criterion, not an existence theorem. It makes no
claim that a candidate passes the square test or the final endpoint test.
Its new step is to reverse the physical projection-content identity from
the [reciprocal stretch theorem](mobius_reciprocal_stretch_grid.md), instead
of treating its integral projection as a consequence of an already chosen
map. In particular the search has no independent shear parameter of size
`N^(1/4+O(1/m))`.

## 1. Source data and explicit bounds

Let `z_0,...,z_(m-1)` be a primitive Gaussian integral circle tuple of
common norm `N`, with `m>=8`, in an arc of phase width `Theta` satisfying

```text
Theta N^(1/4) <= 2.
```

Primitive means the complete common Gaussian gcd is a unit. All target
radii below are also the least primitive squared radii of their full
phase tuples, with every source label retained.

Use the constants from the
[critical triangular reduction](mobius_finite_sample_critical_form.md):

```text
k = floor(m/4),           F = floor((m-1)^2/4),
Qm = floor(m^2/4),
q = m(m-1) / (2k((m-1)(m-2k-1)+1)),
g = m/(4F),
ell = 1/(4(m-1)) + Qm q/(m(m-1)),
b_m = ell+g,             h = ell+2g,
D_m = 4 * 2^(m(m-1)/F + m b_m),
u = min_(2<=n<=m-1) [h_n+h_(m-n+1)],
h_n = floor(n^2/4)/(n(n-1)).
```

Here `0<q<1`. Set

```text
rho = 1-4/m,
X = (16N)^(1/rho),
Emax = 2^(m q) N^q,
Amax = D_m X^h,
Zmax = 5 * 2^(u+m q) D_m^2 X^(u+2h) N^(q-1/2).
```

The letters `E,a,Z` below are arithmetic search parameters; `Z` is a
Gaussian integer, distinct from a physical source point `z_j`.

## 2. The criterion

For each actual source anchor `j`, make the following finite search:

1. Choose a positive divisor `E|N` with `E<=Emax`.
2. Choose `Z in Z[i]` with `0<|Z|<=Zmax`.
3. Choose an integer `1<=a<=Amax`.
4. Form the exact rational Gaussian number

   ```text
   W = Z conjugate(z_j)/E.
   ```

   Require `W` to be Gaussian integral. Write `W=r+i s`. Require

   ```text
   b = s/(2a) in Z,
   d^2 = a^2-b^2-r
   ```

   to be a positive integer square, with its positive root `d<=Amax`.
   Require `gcd(a,b,d)=1`.
5. In the relative source frame at `z_j`, apply

   ```text
   M = ((a,b),(0,d)).
   ```

   Compute the full primitive image squared radius `N'` and its shortest
   containing phase-arc width `Theta'`. Accept precisely when

   ```text
   Theta' (N')^(1/4) <= 2.
   ```

**Criterion.** A rational nonconformal projective half-angle map from the
full source to an endpoint image of constant at most two exists if and
only if this search accepts at least one candidate.

The inequality in the last step is an actual arc inequality, not a test
on a hypothetical radius or on an upper bound for the radius. It may be
evaluated with certified angular intervals. A coarser chord test alone
is not substituted for it.

There is no need to assert that `N/E` equals the clipped conductor of
every candidate. It equals that conductor for the representative used in
the necessity proof. Candidates that have a different clipped conductor
are harmless: acceptance is checked on their literal full images.

## 3. Why every endpoint map occurs in the search

Suppose an endpoint image exists, with least squared radius `N'`.
Conjugate the target if needed to make the map orientation-preserving;
this preserves its radius and phase-arc width. The
[radius-unordered retention bound](mobius_arbitrary_unit_retention_addendum.md)
applied in both directions gives

```text
N' >= N^rho/16,       N >= (N')^rho/16.
```

Put `P=max(N,N')`; thus `P<=X`. Apply the critical triangular reduction
with the larger realization as its source. If that is the original
target, take the triangular adjugate afterward. An integral triangular
matrix and its adjugate exchange the diagonal entries and negate the
shear; their determinants, traces, condition numbers, and discriminants
are the same. Therefore, in either radius ordering, at a corresponding
actual **original source** anchor there is a primitive representative

```text
M=((a,b),(0,d)),       1<=a,d<=D_m P^h.
```

This uses the same physical map, or its inverse while applying the
reduction. It is not the operation of reanchoring the source and then
applying a fixed diagonal map.

Write `delta=ad` and `T=a^2+b^2+d^2`. The finite-sample bound gives
`kappa<=4(2P)^u`, so

```text
T = delta(kappa+kappa^(-1))
  <= 5 * 2^u D_m^2 P^(u+2h).                         (1)
```

Let `Q` be the clipped conductor. The reciprocal stretch theorem proves
that it is the same for the actual-anchored inverse and that

```text
Q | gcd(N,N'),
Q >= 2^(-m q) P^(1-q).
```

Consequently `E=N/Q` is a source divisor with

```text
E <= 2^(m q) N/P^(1-q) <= 2^(m q) N^q = Emax.       (2)
```

The Gram anisotropy and physical projection are

```text
Wgram = a^2-b^2-d^2+2iab,
C = Wgram z_j.
```

The exact projection-content theorem says that `Q` divides both integer
coordinates of `C`. Hence the following is a Gaussian integer:

```text
Z = C/Q.
```

Since the map is nonconformal, `Wgram!=0`, so `Z!=0`. Moreover
`|Wgram|=sqrt(T^2-4delta^2)<T`, whence (1)--(2) give

```text
|Z| = |Wgram| sqrt(N)/Q
    < T E/sqrt(N)
    <= 5 * 2^(u+m q) D_m^2 X^(u+2h) N^(q-1/2)
    = Zmax.                                         (3)
```

Finally `Norm(z_j)=N=QE` proves the exact reverse identity

```text
Z conjugate(z_j)/E = Wgram.                          (4)
```

Its imaginary coordinate recovers `b`, and its real coordinate gives
the displayed square equation for `d`. Thus this representative appears
in the finite search and passes its final endpoint test.

Conversely every accepted candidate is a nonsingular rational map of
the required source: `a,d>0`. Its Gram anisotropy is exactly the nonzero
Gaussian integer `W` recovered from `Z`, so it is nonconformal. Its
actual image satisfies the requested endpoint inequality by the final
test. This proves sufficiency without any extra conductor assumption.

## 4. Exact target reconstruction

For each source phase `z_i/z_j`, choose an ordinary coordinate-primitive
integer half-angle row `H_i`, including `H_j=(1,0)`. Odd-odd rows are
allowed. Put `Y_i=M H_i` without requiring these raw rows to be primitive.
For every pair set

```text
D_ik = Y_i dot Y_k,
F_ik = det(Y_i,Y_k),
c_ik = gcd(D_ik,F_ik),
epsilon_ik = 2 if D_ik/c_ik and F_ik/c_ik are both odd, else 1,
n_ik = (D_ik^2+F_ik^2)/(c_ik^2 epsilon_ik).
```

Then the exact least squared radius is

```text
N' = lcm_(i<k) n_ik.
```

All pairs are needed; an anchor-star lcm can lose opposite Gaussian
orientations. The final phase-arc computation uses `Y_i/conjugate(Y_i)`.
Common row contents cancel from these phases and the pair formula, so no
ordinary or Gaussian cancellation is omitted.

## 5. A source-radius prefilter by close complementary factors

The search has a necessary scalar test which does not use any source
anchor. Taking norms in (4) and using the discriminant identity gives

```text
E^2 Norm(Wgram) = N Norm(Z),
(ET)^2-N Norm(Z) = 4(Ead)^2.
```

Therefore, with `S=Norm(Z)`, the two positive integers

```text
Fminus = E(T-2ad),       Fplus = E(T+2ad)
```

satisfy

```text
Fminus Fplus = N S,
0 < Fplus-Fminus = 4Ead <= 4 Emax Amax^2,
1 <= S <= Zmax^2,        S is a sum of two integer squares.   (5)
```

Positivity of `Fminus` follows from nonconformality and `a,d>0`:
`T-2ad=(a-d)^2+b^2>0`. Thus if no positive sum-of-two-squares
multiplier in this range gives such a close complementary factor pair,
no full endpoint image exists. This is a necessary prefilter only; it
does not replace the Gaussian integrality, square, and endpoint tests.

The multiplier and gap powers of `N` in (5) are respectively

```text
2q-1+2(u+2h)/rho = 27/m+O(1/m^2),
q+2h/rho = 21/(2m)+O(1/m^2).
```

It is the small residual `Z` that makes (5) a restriction on a small
multiple of the original radius `N`. Divisibility of the large
triangular discriminant alone does not directly give this statement.

## 6. The number of small parameter choices

There are at most

```text
m floor(Emax) floor(Amax) (2floor(Zmax)+1)^2             (6)
```

raw choices: the divisor and disk conditions only reduce this box count.
This is a count of candidate arithmetic data, not a claim about factoring
or bit-operation complexity. For fixed `m`, its power of `N` is

```text
v_m = 3q-1 + (2u+5h)/rho
    = 137/(4m) + O(1/m^2).                            (7)
```

In particular all three searched parameters have size `N^O(1/m)`.
For an elementary bound with absolute constants, when `m>=256` the
established inequalities `q<=8/m`, `h<=7/m`, `u<=1/2+2/m`, and
`D_m<=2^14` give

```text
Emax <= 2^8 N^(8/m),
Amax <= 2^15 N^(8/m),
Zmax <= 2^43 N^(27/m).
```

For the last inequality, the exponent of `N` is at most
`8/m+18/(m-4)<=27/m`; the power of two is less than 43.
Substitution into (6), using `2floor(Zmax)+1<=3Zmax`, bounds the
number of raw choices by

```text
m * 2^113 * N^(70/m).                                (8)
```

The constants deliberately leave room for rounding. Formula (6) with
the exact constants is substantially smaller. The scalar prefilter has
the corresponding coarse bounds `S<=2^86 N^(54/m)` and
`Fplus-Fminus<=2^40 N^(24/m)`.

The criterion still does not construct an accepted candidate for an
arbitrary source. The remaining existence question is whether the small
source divisors and Gaussian integers admit the simultaneous integrality,
square, and actual endpoint conditions. Near-conformal approximation or
minimization of the unnormalized row-norm product does not supply those
conditions. The criterion isolates that question without independently
searching a shear of size `N^(1/4+O(1/m))`, an unknown target radius, or
three unrestricted target directions.

The [exact checker](check_endpoint_map_small_parameter_reconstruction.py)
checks the recovery at every actual source anchor of its fixtures, both
radius orders, mixed half-angle parities, nontrivial clipping, ordinary
contents, the full target-radius formula, and the rational exponents.
