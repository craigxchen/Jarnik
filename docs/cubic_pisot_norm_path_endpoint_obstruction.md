# An endpoint obstruction for cubic Pisot norm paths

This note treats a fixed real cubic field with two nonreal embeddings.
It strengthens the rational-offset result in
[the phase note](norm_descent_rational_phase_gap.md) to every fixed
algebraic offset in that field. It uses a published consequence of the
Subspace Theorem, and does not cover arbitrary extracted integer profiles.

## 1. Statement

Let `K` have degree three and signature `(1,1)`, with its real embedding
fixed. Suppose

```text
P(t)=X(t)+iY is monic of degree h,
X in K[t],       Y in K\{0},
P(r)=0 for a fixed r in Q(i),
Norm_(K(i)/Q(i))(P(t)) has constant imaginary part.
```

Let `epsilon in K` be a Pisot algebraic unit, positive at the real
embedding, and fix any `xi in K`. Put

```text
T_n=epsilon^n,
A_n=Norm_(K(i)/Q(i))(P(xi+epsilon^n)),
eta_n=A_n/bar(A_n).
```

The norm is taken after evaluating the element. For a Gaussian rational
`z=p/q` in coprime Gaussian integers, put `H(z)=max(|p|,|q|)`.

**Theorem.** Fix `C>0` and `0<kappa<1/6`. For every sufficiently large
`n`, if a rational unit-circle rotation `u` makes `u eta_n` the ratio
of two nonzero Gaussian integers of common modulus `R` in an origin-centered
arc of length at most `C sqrt(R)`, then

```text
H(u)>T_n^kappa.                                         (1)
```

The threshold depends on the fixed data and is not asserted effective.
Thus rotations of height `exp(o(n))` cannot produce endpoint pairs even
on an infinite sparse subsequence. Common Gaussian content and the actual
circle radius are retained. No linear-splitting hypothesis is needed.

## 2. Published inputs and exceptional sets

We use [Corvaja--Zannier, *A lower bound for the height of a rational
function at S-unit points*, Theorem 1 and formula (1.2)](https://arxiv.org/pdf/math/0311030).
For fixed coprime polynomials in two variables, not both vanishing at
the origin, their full-place logarithmic gcd on a finitely generated
multiplicative group is sublinear in input height outside finitely many
translated proper subtori. Proposition 1, formula (2.2), similarly
bounds negative local logarithms of a fixed Laurent polynomial containing
a nonzero constant monomial. We apply that to `1+c(y/x)^ell`.

The exceptional sets are fixed once a positive error exponent is fixed.
Each will meet our power sequence finitely, so the resulting bounds hold
for every sufficiently large index, not just an infinite subsequence.

Work in a Galois field `E` containing all conjugates and `i`. Write
`epsilon_1=epsilon`, `epsilon_3=bar(epsilon_2)`. The norm is positive,
so

```text
epsilon_1 epsilon_2 epsilon_3=1,
|epsilon_2|=|epsilon_3|=epsilon^(-1/2).
```

If `epsilon_1^a epsilon_2^b epsilon_3^c` is a root of unity, take
absolute values after automorphisms sending each conjugate to the real
one. This gives `2a=b+c`, `2b=a+c`, `2c=a+b`, hence `a=b=c`.
Thus any two conjugates are multiplicatively independent even modulo
roots of unity. A translate of a proper subtorus satisfies some
`x^a y^b=c_0`, `(a,b)!=(0,0)`. Two indices of
`(epsilon_j^n,epsilon_k^n)`, `j!=k`, in that translate would contradict
this independence. Each such translate therefore contains at most one
index; a finite exceptional set has the same harmless effect.

## 3. Common Gaussian content is subexponential

For the embeddings `sigma_j` of `K`, put

```text
xi_j=sigma_j(xi),
P_j^+(t)=sigma_j(X)(t)+i sigma_j(Y),
P_j^-(t)=sigma_j(X)(t)-i sigma_j(Y),
F^+(x,y)=P_1^+(xi_1+x) P_2^+(xi_2+y)
          (xy)^h P_3^+(xi_3+(xy)^(-1)),
F^-(x,y)=P_1^-(xi_1+x) P_2^-(xi_2+y)
          (xy)^h P_3^-(xi_3+(xy)^(-1)).                (2)
```

These are ordinary polynomials with nonzero constant terms. Indeed,
`P(xi)!=0` because `X(xi)` and `Y` are real and `Y!=0`; every conjugate
is nonzero. The reversed third factor has constant term one by monicity.
The minus signs give the same conclusion.

The polynomials are coprime over an algebraic closure. Their irreducible
factors have forms `x-a`, `y-b`, and `xy-c`, all with nonzero constants.
Different forms cannot coincide. Within one form, a common factor would
give a common root of `P_j^+` and `P_j^-`, whose difference is the
nonzero constant `2i sigma_j(Y)`. Neither polynomial has a factor `x`
or `y`.

At `(x,y)=(epsilon_1^n,epsilon_2^n)`, the two values are

```text
F^+(x,y)=(xy)^h A_n,
F^-(x,y)=(xy)^h bar(A_n).                              (3)
```

Conjugation permutes the coefficient embeddings and changes the sign of
`i`, giving the second identity. The multiplier is an algebraic unit
and has valuation zero at every finite place.

There is a fixed positive integer `D_0` with `D_0 A_n in Z[i]` for all
`n`: clear the fixed polynomial and offset denominators before taking
the norm of an algebraic integer. Theorem 1 applied to (2), with the
exceptional indices excluded as above, gives

```text
log |gcd_G(D_0 A_n,D_0 bar(A_n))|=o(n).                (4)
```

Here multiplying by `D_0` changes the estimate by a fixed constant.
The normalized finite-place sum agrees with the logarithm of the
Gaussian gcd modulus for Gaussian integers, and extending the field
does not change that sum. The full-place gcd height also has nonnegative
archimedean contributions, which may be discarded in this upper bound.
Thus no unproved bounded-content assertion is used.

The growing embedding contributes `T_n^h(1+o(1))`, and the other
factors tend to nonzero limits. It follows that

```text
A_n/T_n^h -> B=Norm(P(xi))/P(xi) !=0,
H(eta_n)=T_n^(h+o(1)).                                (5)
```

## 4. Irrational limiting phase

Set `eta=B/bar(B)`. Expansion gives `eta_n=eta+O(T_n^(-1/2))`.
The fractional-linear identity in the earlier phase note gives

```text
Q(i,eta)=Q(i,U(xi)),       U=X/Y.                       (6)
```

If `eta notin Q(i)`, its relative degree is three. Clearing its minimal
polynomial gives, for rational unit-circle `u`,

```text
|u-eta^(-1)|>=c H(u)^(-3).                             (7)
```

Suppose `H(u)<=T_n^kappa` and an endpoint pair exists. Its ratio
`gamma=u eta_n` has `H(gamma)<=R`; this uses the actual Gaussian
integer points, including any common-factor removal. Hence

```text
R>=H(eta_n)/H(u)>=T_n^(h-kappa-o(1)),
|u-eta^(-1)|<=O(T_n^(-1/2))
                +C T_n^(-(h-kappa)/2+o(1)).            (8)
```

The second inequality is the chord bound and the unit moduli.
Here `h>=2`: a degree-one monic row with a root in `Q(i)` would have
rational real and imaginary coefficients, so its cubic norm would be its
cube and have nonconstant imaginary part. For `kappa<1/6`, both powers
on the right of (8) exceed `3kappa`, since `h>7kappa`. This contradicts
(7).

## 5. Rational limiting phase and contact order

Now suppose `eta in Q(i)`. Equation (6) gives `q=U(xi) in Q`.
The ramified-prime theorem in the earlier phase note forces `xi` to be
irrational. It generates `K` and has two nonreal conjugates. Put

```text
ell=ord_(t=xi)(U(t)-q),       1<=ell<=h.                (9)
```

In fact `ell<h`. Otherwise monicity gives
`P(t)=(t-xi)^h+Y(q+i)`. The root `r` is nonreal because `Y!=0`.
Evaluation and conjugation give

```text
((r-xi)/(bar(r)-xi))^h=(q+i)/(q-i) in Q(i).            (10)
```

This quotient has modulus one at the real embedding, so every embedding
fixing `Q(i)` must give it modulus one. But if `r=u+iv`, `v!=0`, and
`xi_j=a+ib`, `b!=0`, its squared modulus is

```text
[(u-a)^2+(v-b)^2]/[(u-a)^2+(v+b)^2] !=1.
```

The denominator cannot vanish: a degree-three conjugate cannot equal
`bar(r)`, of degree at most two. This contradicts (10).

Let `a_j=sigma_j(U^(ell)(xi)/ell!)`. Then `a_2!=0` and
`a_3=bar(a_2)`. Taylor expansion at the two contracting embeddings gives

```text
eta=((q+i)/(q-i))^2,
eta^(-1) eta_n-1
 =-2i/(q^2+1) [a_2 epsilon_2^(ell n)+a_3 epsilon_3^(ell n)]
    +O(T_n^(-(ell+1)/2)).                             (11)
```

The growing factor has phase error `O(T_n^-h)`, since its imaginary
part is constant. The product of the contracting errors is `O(T_n^-ell)`.
Both fit the remainder because `1<=ell<h`.

Factor the bracket as

```text
a_2 epsilon_2^(ell n)
 [1+(a_3/a_2)(epsilon_3/epsilon_2)^(ell n)].
```

Apply Proposition 1, formula (2.2), to
`1+(a_3/a_2)(y/x)^ell` at `(epsilon_2^n,epsilon_3^n)`. The exceptional
indices are finite. At the selected complex place the last bracket is
at least `exp(-o(n))` in modulus and is bounded above by a constant.
The resulting lower bound dominates the strictly smaller remainder in
(11), proving

```text
|eta^(-1) eta_n-1|=T_n^(-ell/2+o(1)).                  (12)
```

This controls all sufficiently large indices, including possible sparse
subsequences with unusually strong cancellation.

For `u!=eta^(-1)`, Gaussian rational separation gives
`|u-eta^(-1)|>=c/H(u)`. Combining it with the first inequality of (8),
the bound `|eta_n-eta|=O(T_n^-ell/2)` and the chord bound contradicts
`H(u)<=T_n^kappa`: here `ell/2>kappa` and `h>3kappa`.
Thus any such putative pair must eventually have `u=eta^(-1)`.

For that fixed rotation, (5) gives `R>=T_n^(h+o(1))`. Equation (12)
then gives

```text
|gamma-1| sqrt(R)>=T_n^((h-ell)/2+o(1)) -> infinity,
```

contradicting the fixed constant `C`. This completes the theorem.

## 6. Application and limits

The degree-16 example in section 5 of the earlier phase note has a
rational limiting phase at its algebraic offset. The present theorem
shows that aligning that limit cannot produce an endpoint pair: the
primitive radius and finite contact order must also be counted.
Corrections of subexponential height cannot repair that path.

In that example, with `alpha^3=1/2`, exact differentiation gives

```text
U'(alpha)=-987/125+(1599/125) alpha+(4371/250) alpha^2 !=0.
```

Thus `ell=1`, the primitive ratio height is `T_n^(16+o(1))`, and
the normalized chord for the fixed aligning rotation grows as
`T_n^(15/2+o(1))` at the least possible radius. These asymptotics use
the proof above; the exact derivative is also checked in the existing
[arithmetic checker](check_norm_descent_rational_phase_gap.py).

The field, polynomial, unit and offset are fixed. Totally real cubic
fields, higher degrees, varying offsets and arbitrary extracted integer
profiles need further arguments. No improved general point-count growth
or radius-uniform theorem follows from this result.
In particular, the varying prime factors of arbitrary extracted profiles
need not lie in one fixed finitely generated multiplicative group. The
exceptional-set argument above does not provide uniformity over such groups.
