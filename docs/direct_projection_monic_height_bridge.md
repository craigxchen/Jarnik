# Direct source projections give a monic square model, but not a small kernel

Let `z_0,...,z_(m-1)` be distinct Gaussian integers of norm `N`, with
Gaussian gcd a unit. Fix the actual anchor `z_0`. No second configuration
or projective map is assumed. Gaussian primitivity implies that `N` is
odd and has only split prime factors.

The direct source already supplies a common integer parameter for the
[even-product square lemma](runge_even_product_square_bound.md). What is
missing is the small root and arithmetic support height used by the
[conditional growth theorem](reciprocal_squareclass_runge_growth.md).

## 1. Exact primitive pair dictionary

For every nonanchor row choose ordinary coprime integers `a_i,b_i` with

```text
z_i/z_0=(a_i+ib_i)/(a_i-ib_i),
f_i=a_i^2+b_i^2,       d_i=2N/f_i.
```

The antipodal case is represented by `(a_i,b_i)=(0,1)`. The denominator
`f_i` divides `2N`, so `d_i` is a positive integer divisor of `2N`.
Indeed, setting `X_i=Re(conjugate(z_0)z_i)` gives

```text
N-X_i=2N b_i^2/f_i in Z,
gcd(f_i,b_i^2)=1.
```

Define the radial deficit and global parameter

```text
c_i=N-X_i=|z_i-z_0|^2/2>0,       Z=2N.
```

Exactly,

```text
c_i=d_i b_i^2,       Z-c_i=d_i a_i^2.               (1)
```

Thus an even set of labels whose product `product d_i` is a square gives

```text
product (Z-c_i)=(sqrt(product d_i) product a_i)^2.  (2)
```

Every root has at most two physical rows, since it specifies one linear
projection on a circle. Select distinct roots before applying the
even-product lemma. The parity vectors of the `d_i` are supported on
the primes of `2N`; equivalently, `c_i` and `d_i` have the same
squareclass. All half-angle parity factors remain in the integer `2N`.

If the source lies in an arc of length at most `C N^(1/4)`, then

```text
0<c_i<=C^2 sqrt(N)/2.                              (3)
```

In particular, the bounded primitive imaginary residue `b_i` does not
itself improve the polynomial root height: its complementary factor
`d_i` turns `d_i b_i^2` into exactly the physical radial deficit.

## 2. Exact common divisor with the global parameter

Include the anchor deficit `c_0=0`, and set

```text
g_0=gcd_Z(Re z_0,Im z_0),
G=gcd_Z(2N,c_0,...,c_(m-1)).
```

Then

```text
G is either g_0 or 2g_0.                            (4)
```

Here is a proof at arbitrary odd prime powers. For `p^e||N`, write
`e_i=v_pi(z_i)` at one Gaussian orientation. Primitivity of the full
tuple means that the allocations attain both `0` and `e`. The reduced
primitive pair norm has valuation `|e_i-e_0|`, so

```text
v_p(d_i)=e-|e_i-e_0|.
```

If `e_i!=e_0`, then `p|f_i` and ordinary primitivity forces `p` not to
divide `b_i`. Hence `v_p(c_i)=e-|e_i-e_0|`. If `e_i=e_0`, the valuation
of `c_i` is at least `e` (the anchor zero has infinite valuation).
Taking the minimum over all rows gives exactly

```text
min_i v_p(c_i)=min(e_0,e-e_0)=v_p(g_0).
```

No odd prime outside `N` enters `G`; since `N` is odd its two-part is
at most two. This proves (4), with both Gaussian orientations and every
prime-power exponent retained.

Consequently the largest common integer by which one can divide the
parameter `Z` and all these unshifted roots is at most `2g_0`. If the
actual anchor has ordinary content one, that division saves at most a
factor of two. Equation (4) does not rule out other affine changes of
parameter or a different construction. After division by `G`, (2)
remains an integer square identity: an even product divides its square
by `G^(2k)` and rational squares which are integers are integer squares.

## 3. Direct squareclass multiplicity and the remaining bridge

There is a useful direct spacing consequence. In a fixed squareclass
write `d_i=s v_i^2` with `s` squarefree, and put `t_i=|v_i a_i|`.
Then

```text
2N-c_i=s t_i^2.
```

If `D=C^2 sqrt(N)/2<2N`, all the possible nonnegative integers `t_i`
lie in an interval of length

```text
D/[sqrt(s)(sqrt(2N)+sqrt(2N-D))]
 <= C^2/(2sqrt(2)).                                (5)
```

Thus a squareclass has at most
`2(floor(C^2/(2sqrt(2)))+1)` physical rows. This is a bound on
multiplicity in one squareclass, not a bound on how many squareclasses
occur. For `C<=2` and `N>=4`, the sharper denominator using `D<=N`
makes (5) at most `2/(sqrt(2)+1)<1`; each squareclass then has at most
two physical rows.

This is a direct radial proof of bounded squareclass multiplicity. The
earlier [squareclass packing baseline](squareclass_baseline_and_power_groups.md)
already gives a stronger pairwise separation for `C<2`; no new uniform
count for bounded squareclass rank is claimed here.

The direct construction has parameter `Z=2N`, roots bounded by
`O_C(sqrt(N))`, and labels supported on `2N`, whose log height is
`log N+O(1)`. The conditional argument instead has roots bounded by
`H^2`, with `log H=O(log N/m)`, while its parameter still has a fixed
positive power of `N`. That separation excludes squareclass relations
of length proportional to `m`; (3) alone does not provide it.

The [cotangent offset dictionary](integer_cotangent_offset_height_cut_residue_dictionary.md)
gives a related precise obstruction: a common residue modulus with
`log L=o(log N)` forces offset-kernel log height at least
`(1/2-o(1))log N` as the point count grows, rather than
`O(log N/m)`. The current direct bridge would need a new global
normalization or identity supplying both small root height and small
parity-support height. Neither bounded primitive pair residues nor (4)
provides that bridge. No new general endpoint count is asserted.

## 4. An actual root-height floor at every anchor for short arcs

The root cost in (1) has a uniform **lower** bound on almost every row
when the actual primitive source arc has `C<=sqrt(2)` and `N>1`.
This range controls all literal Gaussian units at once. Fix any actual
anchor `z_0` and retain the exact half-angle parity

```text
f_i=epsilon_i n_(0i),       epsilon_i in {1,2},
d_i=2N/f_i,               c_i=d_i b_i^2.
```

The chord inequality in (3) and `|b_i|>=1` give
`f_i>=4sqrt(N)/C^2>=2sqrt(N)`, hence
`n_(0i)>=sqrt(N)`. The same argument applies to **every** pair
`i,j`, so the threshold Gram matrix

```text
mu_ij=1-2 log n_ij/log N,       mu_ii=1
```

is positive semidefinite, has diagonal one, and has nonpositive
off-diagonal entries without selecting a Gaussian-unit cohort. Its
positive semidefiniteness follows from the weighted allocation-sign
Gram identity in
[projective orbit radius rigidity](projective_orbit_radius_rigidity.md):
`mu_ij=E_layers(F_i F_j)` with each `F_i` a unit sign vector.
Write `mu=I-B`, where `B` is symmetric with nonnegative entries and
zero diagonal. Positivity of `mu` gives `lambda_max(B)<=1`.
For any real vector `v`, entrywise nonnegativity gives
`|v^t Bv|<=|v|^t B|v|<=||v||^2`; thus every eigenvalue of `B`
lies in `[-1,1]`, and `||mu||_op<=2`. Consequently

```text
sum_(i!=0) mu_0i^2=(mu^2)_00-1<=2mu_00-1=1.         (6)
```

For any `epsilon>0`, if `n_(0i)>N^(1/2+epsilon)`, then
`mu_0i<-2epsilon`. Equation (6) shows that the number of such rows
is **strictly less than** `1/(4epsilon^2)`. For every other row,
the exact parity factor gives

```text
c_i=(2N/(epsilon_i n_(0i)))b_i^2
    >=N/n_(0i)>=N^(1/2-epsilon).                    (7)
```

The proof works after fixing **any** actual anchor, and an arbitrary
subselection cannot turn more than `1/(4epsilon^2)` rows into
small-root exceptions for that anchor. This is a height floor for
the specific monic square model (1), not a no-go theorem for a
different projection, normalization, or global identity. In
particular it does not construct the small arithmetic kernel required
by the conditional squareclass theorem.

The [checker](check_direct_projection_monic_height_bridge.py) verifies
primitive half-angle divisibility, (1), (4), even-product identities,
and direct squareclass multiplicity on literal Gaussian tuples.
