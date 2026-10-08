# Square descent of the degree-twelve branch form: an applicability audit

The square-value condition has an exact finite descent, with polynomially
many twists in the coefficient height. It does **not**, by itself, give a
height bound on an isolated endpoint solution. Primitive Gaussian realization
adds no extra existence condition to rationality on the nondegenerate cover.
This note does not use Runge or assert that square descent cannot ultimately
be strengthened.

Write `F=Delta_u`, with integer coefficients of maximum absolute value `B>=2`,
where the controlled frame gives `B<=c U^b` for absolute constants. Assume F
is squarefree of degree 12 as a binary form. Fix primitive `(s,t)` with

```
F(s,t)=d^2 !=0,       H=max(|s|,|t|),
0<|d|<=C_end U^a H.
```

All constants below depend only on degree unless explicitly indexed by F.

## 1. The integral etale algebra and every bad prime

We may arrange that `a0=F(1,0)!=0` by an integral unimodular substitution
of bounded height: one of `(k,1)`, `0<=k<=12`, is not a zero. Complete that
column by `(1,0)`. This changes B and H by absolute factors only and
preserves primitivity. Put

```
f(X)=F(X,1),                    a0=leading coefficient(f),
g(X)=a0^11 f(X/a0),             A=Q[X]/(g),
theta=X mod g,                 L=a0 s-theta t.
```

Here g is monic integral, squarefree, and degree 12. The algebra A is a
product of number fields; its maximal order O is the product of their rings
of integers. No irreducibility or total reality is assumed. Exactly,

```
Norm_A/Q(L)=a0^11 F(s,t)=a0^11 d^2.                (1)
```

Let `D=|disc(g)|>=1` and `Q=2|a0|D`. Include in S **all** prime ideals of O
above rational primes dividing Q. Then

```
D=|a0|^110 |disc(f)| <= c B^132,
Q<=c B^133,                  #S<=12 omega(Q).      (2)
```

The discriminant identity follows by scaling the twelve roots; the bound
uses the degree-22 coefficient homogeneity of the discriminant of f.
If `g=product g_j` is its factorization into monic irreducibles, the identity

```
disc(g)=product disc(g_j) product_(j<k) Res(g_j,g_k)^2
```

shows explicitly that S includes the cross-factor resultants, all order
indices, and all field discriminants, in addition to the leading coefficient
and 2. The factor 2 is conservative and harmless.

For `p` outside Q, if a prime over p divides L then `p` cannot divide t:
otherwise primitivity makes `a0 s` nonzero modulo p. Since the order is
etale at p and theta must reduce to `a0 s/t` in F_p, there is exactly one
prime of the entire product O that can divide L. Its residue degree is one.
Taking valuations in (1) therefore proves

```
v_P(L)=2v_p(d) is even                     (P outside S).   (3)
```

This argument accounts for reducible F as well: treating its factors as
independently square-valued without the resultant primes would be incorrect.

## 2. Exact ideal classes and unit factors

In each field component, and hence in the product, there is a unique
factorization

```
(L)=a b^2,      a=product_(P in S) P^(v_P(L) mod 2),           (4)
```

with a squarefree integral S-supported ideal and b integral. For each of
at most `2^#S` choices of a, retain precisely the classes c satisfying
`c^2=[a]^-1` in `Cl(O)`. There are either zero or `|Cl(O)[2]|` such classes.
Choose an integral ideal b0 representing each retained class and choose a
generator gamma of the principal ideal `a b0^2`. Since `[b]=[b0]`, write
`b=(beta)b0`, with `beta in b0^-1`. Equality of principal ideals gives

```
L=epsilon gamma beta^2.
```

Write the unit epsilon as `epsilon0 v^2`, absorbing v into beta. If the
field components have signatures `(r_j,s_j)`, their unit squareclass group
has size `2^(sum_j(r_j+s_j))<=2^12`. Thus the final finite list is exactly

```
a0 s-theta t = eta beta^2,       beta in b0^-1,                (5)
eta=epsilon0 gamma,
number of possible twists <=2^(#S+12) |Cl(O)[2]|.              (6)
```

Equation (5) must retain (1), primitivity, and the rational-square norm
condition: (3) alone is only a necessary condition. Conversely, retaining
these conditions gives exactly the original nonzero square-value solutions.
In particular, replacing the unit epsilon by a fixed bounded unit rather
than absorbing its square into beta is invalid.

The only external facts used here are the standard ideal-class and unit
theorems, checked in [Milne, Algebraic Number Theory, Theorems 4.3, 4.4 and
5.1](https://www.jmilne.org/math/CourseNotes/ANT.pdf): integral class
representatives have norm at most the Minkowski constant times the square
root of the field discriminant, the class group is finite, and the unit
rank is `r+s-1`.

For clarity, the list size is `B^O(1)=U^O(1)`, independent of H. Indeed
`2^#S<=Q^12`. Minkowski bounds the norm of b0 in a degree-n component by
`c_n sqrt(|D_j|)`. Counting all sublattices of that index in `Z^n` by
Hermite normal form gives a fixed-power bound for its class number.
Since `product |D_j|<=D` and the total degree is 12, multiplying gives
`|Cl(O)|<=c D^C` for an absolute C. This deliberately crude estimate is
sufficient for (6); it supplies no saving on H.

## 3. What the archimedean and norm conditions actually say

Fix the representatives in (5) and define their explicit archimedean cost

```
E_F=max_(eta,sigma) max(1, |sigma(eta)|, |sigma(eta)|^-1).
```

For the initially chosen generators and unit representatives, polynomial
list size alone does not bound E_F. However, the later
[polynomial squareclass representative theorem](six_point_polynomial_squareclass_representatives.md)
repairs this issue: replace eta and its accompanying ideal together using
a small anti-invariant element of the quadratic extension. It proves that
one can choose `E_F=B^O(1)=U^O(1)` and a polynomially bounded integral
denominator for beta. No regulator bound is needed. The formulas below
apply after this replacement as well. For every embedding sigma, Cauchy's
root bound gives `|sigma(theta)|<=2B`,
so (5) gives the exact formula and upper estimate

```
|sigma(beta)|^2=|a0 s-sigma(theta)t|/|sigma(eta)|,
|sigma(beta)| <= (3 B E_F H)^(1/2).                           (7)
```

Its product is exactly

```
|Norm(beta)|=|a0|^(11/2) |d| / |Norm(eta)|^(1/2).              (8)
```

Thus norm discreteness offers a lower bound on |d| independent of H; the
original integer condition already gives `|d|>=1`. It does not turn the
upper bound `|d|<=C_end U^a H` into an upper bound on H.

For a fixed F, an endpoint solution approaching a simple real root has
one factor `|s-alpha t|=O_F(C_end^2 U^(2a) H^-9)` and the other eleven
factors comparable to H, once H is sufficiently large relative to the
root-separation constants and `C_end U^a`. In (5) this becomes one small
embedding of beta of order at most `O_(F,C_end,U)(H^-9/2)`, while the
other embeddings have order `H^1/2`. Their product still has size O(H),
exactly as (8) permits. Indeed the fixed fractional-ideal denominator and
the other eleven upper bounds only force the small embedding to be
`Omega_F(H^-11/2)`, which is compatible with `O_F(H^-9/2)`. This is an
unbalanced algebraic integer problem,
not an integrality contradiction. The constraint that `eta beta^2` lies
in the rational two-plane `span_Q(a0,theta)` retains the difficult
Diophantine content; discarding it loses the problem, and retaining it
has not supplied a moving-weight approximation bound.

## 4. Gaussian primitivity is already encoded

On the nonsingular distinct-point locus, a rational solution of the cover
gives the rational Gaussian tuple Z_i by the explicit formula in
[six_point_isotropic_circle_cover.md](six_point_isotropic_circle_cover.md),
with equal positive norms and the original central weights. Clear one
common denominator D0 and divide all six coordinates by their common
Gaussian gcd G. Then

```
z_i=D0 Z_i/G in Z[i],        gcd_G(z_0,...,z_5)=1,
|z_i|^2=D0^2 |Z_i|^2/Norm(G).
```

This construction always exists. It preserves angular positions and the
central moment identities. Gaussian division is a common similarity, not
an additional equation on `(s,t,d)`. It realizes the least radius for the
Gaussian rational projective tuple by the Bezout argument already proved
in the cover note. Particular fair-core conductor or valuation hypotheses
may select a smaller subset, but **primitive Gaussian realizability alone**
does not. Also, the necessary window with a positive power of U must not
be claimed equivalent to a fixed endpoint bound; the moving-radius note
records a polynomial gap between its necessary and sufficient windows.

The remaining possible gain would have to use the two-plane equations,
additional proven conductor restrictions, or a new uniform approximation
estimate. The finite square descent above supplies none of these estimates
and leaves the uniform circle-arc bound unresolved.
