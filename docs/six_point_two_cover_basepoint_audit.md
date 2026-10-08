# The multiplication-by-two cover: base points and controlled twists

Changing the Abel--Jacobi base point changes the descent of the cover by
an exact Kummer class. A half of the difference of the base points gives
the isomorphism; it need not be rational. Nevertheless, large coordinates
of that half do **not** imply an uncontrolled field or uncontrolled twist
coefficients. Over a controlled branch-splitting field one can replace the
target-dependent radical values by squareclass representatives of height
`U^O(1)`. The field, pole, and place costs of this replacement remain.

Throughout, degrees and genus are fixed, and constants in `O(1)` may be
very large but are independent of the point height `H`.

## 1. Exact change of base point

Let `C/k` be a smooth genus-five curve, `J=Jac(C)`, and `A,B in C(k)`.
Define the covers as fiber products

```text
X_A = {(P,Q) in C x J : 2Q=[P-A]},
X_B = {(P,Q) in C x J : 2Q=[P-B]}.
```

They are geometrically connected etale covers of degree `2^10=1024`.
Their deck group scheme is `J[2]`, which need not be constant over `k`.
This is the isogeny-pullback construction in Section 7.2.6 of
[Poonen, Lectures on rational points on curves](https://math.mit.edu/~poonen/papers/curves.pdf).

Write `D=[B-A]`, and choose `T in J(kbar)` with `2T=D`. Then

```text
(P,Q) |-> (P,Q-T)                                      (1)
```

is an isomorphism `X_A -> X_B` over `k(T)`. Its field degree is at most
1024. For `sigma in Gal(kbar/k)`, the difference between its conjugate
and itself is translation by `-(sigma T-T)`. Thus its descent obstruction
is precisely the Kummer cocycle of `D`.

As covers over `C`, they are isomorphic over `k` exactly when

```text
[B-A] belongs to 2J(k).                                (2)
```

Indeed every geometric isomorphism over `C` differs from (1) by a deck
translation, so it is translation by a different half of `D`; descent
means that half is rational. This also follows from the injective Kummer
map `J(k)/2J(k) -> H^1(k,J[2])` (Poonen, Section 4.5).

For a target point `P`, the fiber of `X_A` is the finite etale scheme
`[2]^{-1}([P-A])`. The fiber of `X_P` above `P` contains `(P,0)`.
Identifying this target-based cover with `X_A` requires a half of
`[P-A]`, exactly the same division problem as lifting `P` to `X_A`.
Choosing `A=P` therefore supplies no rational identification with a
previously fixed cover.

When `J[2]` is nonconstant, individual residue fields of this fiber need
not have power-of-two degrees. The bound 1024 is the degree of the
whole etale fiber, not a claim that its residue fields form a
multiquadratic extension of `k`.
Under the additional hypothesis that the linear Galois image is the full
`S_12`, `s12_even_subset_h1.md` gives the possible orbit lists
`(1,66,495,462)`, `(12,220,792)`, and `(1024)`. The last case requires
the full translation kernel. That hypothesis must be checked over the
actual base field: adjoining a branch point changes the linear image.

## 2. A controlled Weierstrass base

In the controlled frame write the affine equation as

```text
y^2=f(x),       deg f=12,       f in Z[x],
height(f) <= U^c,       disc(f)!=0.
```

A fixed small projective change of parameter can ensure that infinity
is not a branch point, if needed. For any root `alpha`, the Weierstrass
point `A=(alpha,0)` is defined over `k=Q(alpha)`, of degree at most twelve
and absolute multiplicative height `U^O(1)`. Consequently `X_A` is a
controlled construction over this field; no rational base point of
controlled height has been assumed.

Let `a` be the leading coefficient. The root `a alpha` is integral,
annihilated by the monic polynomial `a^11 f(X/a)` of controlled height.
Root bounds and discriminants show that both this root field and the
full splitting field `K` have discriminant `U^O(1)`, with
`[K:Q]<=12!`. For the splitting field, one can alternatively use the
bounded-degree ramification argument in Section 3 below. The extension
is unramified outside primes dividing `a disc(f)`.

Over `K`, label the roots `alpha_1,...,alpha_12` and take
`A=(alpha_12,0)`. There is an explicit normalized radical model of `X_A`:

```text
c_i=(alpha_12-alpha_i) f'(alpha_12),
phi_i=(x-alpha_i)/((x-alpha_12)c_i),       1<=i<=10,
K(X_A)=K(C)(sqrt(phi_1),...,sqrt(phi_10)).          (3)
```

Each divisor is `2(W_i-A)`. These ten classes form a basis of geometric
two-torsion. At `A`, use the local parameter `y`; since

```text
x-alpha_12 = y^2/f'(alpha_12) + O(y^4),
phi_i = y^(-2)(1+O(y^2)),
```

the fiber above `A` splits completely over `K`. This normalization
identifies (3) with the based maximal elementary abelian unramified
cover: geometrically it has all ten independent quadratic subcovers,
and the split fiber at `A` removes the constant Kummer twist.
All coefficients in (3) have height `U^O(1)`.

For an affine target away from the branch points, the lift field over
`K` is exactly

```text
K(sqrt(phi_1(P)),...,sqrt(phi_10(P))).              (4)
```

Its degree is `2^r`, where `r` is the rank of these ten **evaluated**
squareclasses in `K*/K*2`. No lower bound on `r` follows merely from
the geometric independence in (3).

## 3. The field cost is independent of the target height

There is a nonzero integer `E=U^O(1)`, independent of `P`, such that the
Jacobian has good reduction away from `E`; one may include
`2a disc(f)` and the frame denominators in `E`.

If `v` is a finite place of a base field away from `E`, its Jacobian
extends to an abelian scheme over the valuation ring. Every rational
Jacobian point extends to a section by properness. Multiplication by
two is finite etale there, so its fiber above this section is a finite
etale scheme over the valuation ring. Therefore every half-point field
is unramified at `v`. This proof imposes no bound on the coordinates of
the section and no affine-integrality assumption on `P`.

In particular, for the controlled root field `k` above and any rational
target `P`, a field of one half has

```text
[k(T):Q] <= 12*1024,
k(T)/Q unramified outside primes dividing E,
|disc(k(T))| <= U^O(1).                            (5)
```

Here the second assertion includes the controlled ramification of `k`.
The same statements, with a different fixed degree bound, hold after
adjoining all branch roots or all conjugate halves.

For completeness, the last implication has an explicit uniform bound:

```text
v_p(disc(L)) <= n(1+floor(log_p n)),
|disc(L)| <= lcm(1,...,n)^n rad(E)^n                (5a)
```

whenever `[L:Q]<=n` and `L` is unramified outside the divisors of `E`.
To prove it, a completion of ramification degree `e` and residue degree
`f` is totally ramified over its maximal unramified subfield. A uniformizer
has an Eisenstein polynomial of degree `e`. The nonzero terms of its
derivative have valuations distinct modulo `e`, so cancellation is
impossible; the leading term bounds the different exponent by
`e-1+e v_p(e)`. Multiplying by `f` and summing over completions bounds
the discriminant exponent by `n(1+floor(log_p n))`. Multiplication over
the ramified primes proves (5a), since
`product_(p<=n) p^floor(log_p n)=lcm(1,...,n)`.
This keeps wild primes explicitly. The local different bound is recorded
in equation (4.12), following Example 4.25, of
[Conrad, The different ideal](https://kconrad.math.uconn.edu/blurbs/gradnumthy/different.pdf).

The geometry-of-numbers proof of Hermite's theorem also supplies a
primitive integral generator with defining polynomial of height a fixed
power of the discriminant. One way to see the polynomial dependence is
to choose `n` independent integral elements with embeddings bounded by
`C_n sqrt(|disc|)` using successive minima. A bounded integral linear
combination separates the `n` embeddings: avoid their at most
`n(n-1)/2` equality hyperplanes. Its conjugates, hence its polynomial
coefficients, obey a fixed power bound. Thus (5) gives a defining
polynomial of height `U^O(1)` and at most `U^O(1)` possible isomorphism
classes of such bounded-degree fields.

This is a field bound, not a bound for the coordinates of `T` as a point
on `J`, and not a bound for the target parameter `H`.

## 4. Controlled representatives of target-dependent twists

There is a further useful consequence over the splitting field `K`.
Let `b_i=phi_i(P)`. Each quadratic extension `K(sqrt(b_i))/K` is
unramified away from the fixed bad set by Section 3. It therefore has
bounded degree and absolute discriminant `U^O(1)`.

If the extension is nontrivial, choose a small integral primitive
generator `beta` for it over `Q`, and let `sigma` be its nontrivial
automorphism over `K`. Then

```text
delta=beta-sigma(beta) != 0,
a_i=delta^2 in K*,
a_i/b_i in K*2,       height(a_i)<=U^O(1).          (6)
```

The nonvanishing holds because a primitive generator cannot be fixed
by `sigma`. Its square belongs to `K`, and generates the same quadratic
extension; this gives the squareclass equality. Height subadditivity
and the small-generator bound prove the final estimate. If `b_i` is a
square already, take `a_i=1`.

The normalized twist

```text
z_i^2=phi_i/a_i,       1<=i<=10,                    (7)
```

now has a `K`-rational point over the target `P` and coefficients of
height `U^O(1)`. Choosing a bounded-height basis of `K`, coefficients
can also be written in rational coordinates of height `U^O(1)` by
linear algebra and the discriminant bound. This gives a finite family
of controlled presentations over `K`. The actual square roots
`sqrt(b_i/a_i)` at the target may still have large height.

Thus raw factors `phi_i(P)`, whose displayed heights involve `H`, are
not an intrinsic coefficient-height obstruction to normalized twists.
Unlike the triangle radicals, these functions define covers unramified
over the complete curve; that is what confines their quadratic fields
to the fixed bad-prime set.

This conclusion takes place over `K`. It neither asserts a rational
identification in (1) over `Q`, nor supplies small rational equations
with all descent and marking data for every `Q`-twist. Nor does it
remove the archimedean places, bad places, or pole fields introduced
by `K`. A Runge argument on (7) must use those actual data.

## 5. What base-point choice does and does not establish

Changing controlled base points translates the Kummer class by the
fixed class of their difference. Forcing a large specialization degree
would require an additional assertion about these translated classes
and their Galois orbits. Neither the Abel--Jacobi construction nor its
geometric degree supplies such an assertion. At its own base point a
cover even has a rational lift, though a branch base is outside the
nondegenerate endpoint locus and is not by itself a counterexample to
a possible arithmetic theorem restricted to that locus.

Conversely, choosing the target as base point is not ruled out simply
because its naive equations contain `H`: Sections 3--4 provide a
controlled-field and controlled-twist replacement after bounded-degree
extension. Whether its pole orbits and required places yield a strict
Runge surplus is a separate question. The large constant field already
falls within the scope of the constant-torsion archimedean audit when
that audit's Galois and triangle-orbit hypotheses hold.

The proved outputs are the translation identity (1), its rationality
criterion (2), the division-field bound (5), and the normalized-twist
construction (6)--(7). They do not establish high specialization degree
or a uniform endpoint height bound, and do not exclude all nonconstant
group covering arguments.
