# Extension auxiliary circle: exact identity and similarity audit

## Status

The extension discriminant really does define an auxiliary integral circle, and
the relevant algebraic identities are formalized in
`GaussianChain/ExtensionAuxCircle.lean`.

However, the earlier interpretation as a primitive-norm descent was incorrect.
After the affine map is written in the original displacement coordinates, it is
exactly a dot-cross similarity.  It scales the entire original circle and arc by
the same linear factor and therefore creates no endpoint gain.

The similarity identity is formalized in
`GaussianChain/AuxiliarySimilarity.lean`.

---

## 1. The discriminant is an auxiliary circle

For a geometric certificate `C`, the extension discriminant is a quadratic in
`W`:

```text
Delta(W) = -A^2 W^2 + L W + C0,
A = 2 y z V.
```

The negative-square leading coefficient is the formalized identity

```text
discLeading = -4 y^2 z^2 V^2.
```

If `Delta(W)=Y^2`, completing the square gives

```text
(2 A Y)^2 + (2 A^2 W - L)^2 = L^2 + 4 A^2 C0.       (1.1)
```

Thus every possible extension maps to an integral point on one fixed auxiliary
circle.  This identity is exact.

---

## 2. Reduced coordinates for an integral slope

Suppose the third rooted point has rooted coordinates

```text
(k,s) = (z,z r).
```

The certificate identities imply

```text
z(r^2+1)=2hN.
```

For an extension with rooted coordinate `t` and adjacent determinant index `W`,
define

```text
X = W-hzr,
Y = h(2t-z)-rW.
```

Then

```text
X^2+Y^2=2h^3Nz.                                     (2.1)
```

This is the theorem

```text
GeomCertificate.extension_reducedCircle
```

in `GaussianChain/ExtensionAuxCircle.lean`.

---

## 3. Displacement-coordinate interpretation

Let the primitive root vector be `(A,B)`, and let the third chord displacement
be

```text
d=(m,n).
```

For a later point, write its displacement from the third point as

```text
Delta=(dm,dn).
```

The rooted-coordinate difference is

```text
t-z = -(A dm+B dn),
```

and the determinant index is

```text
W = m dn-n dm.
```

The third-point relations imply

```text
m = -2hA+r n,
n = -2hB-r m.
```

Therefore

```text
2h(t-z)-rW = m dm+n dn.                              (3.1)
```

After subtracting the image of the third point, the auxiliary map is exactly

```text
Delta |-> (cross(d,Delta), dot(d,Delta)).             (3.2)
```

The Brahmagupta identity gives

```text
cross(d,Delta)^2+dot(d,Delta)^2
  = |d|^2 |Delta|^2.                                 (3.3)
```

Thus the auxiliary transformation is a Euclidean similarity with scale

```text
|d| = sqrt(2hz).
```

Theorems `dot_cross_norm_identity`, `dot_cross_sqDist_identity`, and
`extensionDifference_is_similarity` formalize this calculation in
`GaussianChain/AuxiliarySimilarity.lean`.

---

## 4. Why there is no descent

Equation (2.1) factors numerically as

```text
2h^3Nz=(hz)^2(r^2+1).
```

The image of the third point itself is divisible by `hz`, but later image
points need not share that factor.  Divisibility of the norm by `(hz)^2` does
not imply coordinatewise divisibility at split Gaussian primes.

For example, on the circle of norm `5`, one can choose rooted data for which the
third point maps to `(-5,5)` while a later point maps to `(-1,7)`.  Both lie on

```text
X^2+Y^2=50,
```

but the latter has no common factor `5`.

Hence one cannot divide the whole transformed family by the scale of the base
point.

More conceptually, (3.2)--(3.3) show that no such division should be expected:
the transformed configuration is simply a scaled and rotated copy of the
original configuration.  Angular separations are unchanged.  If the original
arc has endpoint angular width, so does the image after normalization by its
new radius.

---

## 5. The parameter lambda remains useful only pointwise

For one rooted point, writing

```text
g=gcd(k,s),
k=gq,
s=gr,
gcd(q,r)=1
```

gives `q | g`; writing `g=lambda q` yields

```text
k=lambda q^2,
s=lambda qr,
q^2+r^2=2hN/lambda.                                 (5.1)
```

For a fixed `lambda`, the small coordinate `q` lies in an axis-aligned endpoint
range and only `O_C(1)` integer values of the complementary coordinate `r` are
possible.  Thus a fixed `lambda` fiber is bounded.

But the auxiliary-circle map does not globally reduce all different lambda
fibers to a smaller common circle.  It merely repackages them by a similarity.
Any unbounded family would still have to use many distinct lambda values.

---

## 6. Exact rational parametrization and the integral remainder

The discriminant calculation can be sharpened to an exact birational
identification.  Put

```text
r = n3num/J
```

when the third rooted slope is integral.  The geometric certificate identities
give

```text
N V^2-2 h y z q = -2 y z V r,
z(r^2+1)=2hN.
```

After dividing the extension polynomial by `yV`, its equation is

```text
2hzt^2-2z(rW+hz)t+NW^2=0.                            (6.1)
```

For `W != 0`, set

```text
m = 2h(t-z)/W-r.
```

Then (6.1) is parametrized by

```text
t = 2hN/(m^2+1),
W = 2hz(r-m)/(m^2+1).                                (6.2)
```

Conversely every rational `m` gives a rational point of (6.1).  Thus the
extension conic has no rational rigidity beyond the original rooted circle.
For positive geometric data, the rational sign and strict-admissibility
conditions reduce to `0<m<r` together with the endpoint bound on `t`.

Write `m=a/b` in lowest terms and put `e=a^2+b^2`.  Formula (6.2) gives

```text
t = 2hN b^2/e.
```

Since `gcd(e,b)=1`,

```text
t is integral  iff  e divides 2hN.                  (6.3)
```

The remaining reconstruction conditions are explicit linear divisibilities
for `W,w',c',J'`.  In particular, a fixed certificate has at most

```text
sum_(e|2hN) r_2(e) <= 4 tau(2hN)^2
```

raw rational-slope candidates.  This is a useful exact finite search, but its
size is not uniform in the conductor.

For the known certificate

```text
(h,N,z,r)=(5,22698161,1597,377),
```

an exact enumeration of the primitive sum-of-two-squares divisors of `2hN`
finds no candidate with

```text
1597<t,   4t^2<=22698161.
```

So the known four-point cluster has no fifth admissible extension.  This is a
property of that finite divisor list, not a general four-to-five obstruction.

Equivalently, if one instead regards the extension polynomial as a quadratic
in `W`, its discriminant is

```text
4y^2z^2V^2 t(2hN-t).
```

The rooted circle equation makes this a square automatically.  Only the
integrality congruences remain.

---

## 7. Verdict

The auxiliary-circle identities are correct and useful for exact reconstruction
and finite search.  They do not provide a descent mechanism for the uniform
arc theorem.

The precise outcome is:

```text
fixed lambda fiber  -> bounded,
auxiliary circle     -> exact similarity,
full family          -> no reduction in endpoint scale.
```

This route is therefore exhausted as a standalone proof strategy.  The valid
Lean lemmas are retained because they expose the geometry cleanly and prevent
the same false descent interpretation from recurring.
