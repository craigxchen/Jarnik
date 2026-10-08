# Rational values and correction height in Pisot norm descent

This note excludes one combination of two previously investigated
mechanisms: constant-imaginary polynomial norm descent and a parameter
path with one growing embedding and all other embeddings contracting to
a rational value. It also gives a positive height cost for rational
rotations that try to repair that path. It does not give a residual-height
gap for arbitrary endpoint configurations or arbitrary algebraic-unit
paths.

## 1. A ramified-prime obstruction to rational values

Let `K` be a real number field of degree `d>1`, and suppose

```text
P(t)=X(t)+iY,       X in K[t],       Y in K\{0},
P(r_0)=0,          r_0 in Q(i),
Q(t)=Norm_(K(i)/Q(i))(P(t))=A(t)+iC,
A in Q[t],         C in Q.
```

Set `U=X/Y`. Then

> `U(c)` is not rational for any rational `c`.

The [fixed-root parity theorem](fixed_gaussian_root_norm_descent_parity.md)
first gives odd `d>1` and

```text
C=(-4)^((d-1)/2) Norm_(K/Q)(Y).                     (1)
```

If `U(c)=a in Q`, evaluation at the rational parameter commutes with
the coefficient norm. Thus (1), after dividing by the nonzero norm
of `Y`, would give

```text
(a+i)^d-(a-i)^d=(2i)^d.
```

Write `a=A_0/B_0` in lowest terms, with `B_0>0`. Equivalently,

```text
(A_0+iB_0)^d-(A_0-iB_0)^d=(2iB_0)^d.               (2)
```

Let `pi=1+i` and normalize `v_pi(pi)=1`. If two Gaussian integers
`x,y` are units at `pi`, the quotient

```text
(x^d-y^d)/(x-y)=sum_(j=0)^(d-1) x^(d-1-j)y^j
```

is a unit: every summand is one modulo `pi`, and `d` is odd.

If `A_0,B_0` have opposite parity, both bases on the left of (2)
are units. Consequently the left valuation is
`s=v_pi(2iB_0)=2+2v_2(B_0)>0`, whereas the right valuation is `ds`.
This is impossible.

If both are odd, each base has valuation one. After dividing them
by `pi`, the two resulting units have difference of valuation one.
The same quotient argument now gives left valuation `d+1`, while
the right valuation is `2d`, again impossible for `d>1`.
These exhaust coprime integer pairs, proving the assertion. This is
a direct local calculation; no general Fermat theorem is being used.

## 2. The limiting phase on a contracting parameter path

Assume `P` is monic of degree `h`. Let `epsilon in K` be an algebraic
unit with the following properties at the distinguished real embedding:

```text
epsilon>1,
rho=max_(sigma!=id)|sigma(epsilon)|<1.
```

Thus the chosen unit is a Pisot unit, with all embeddings counted;
`K` need not be Galois or totally real. Fix `c in Q`, and define

```text
t_n=c+epsilon^n,
A_n=Norm_(K(i)/Q(i))(P(t_n)) in Q(i),
T_n=epsilon^n,
delta=min(1,-log(rho)/log(epsilon))>0.
```

Here `A_n` is the norm of an evaluated algebraic element. In general
it is not `Q(t_n)`: norm and evaluation do not commute at a parameter
outside `Q(i)`.

At the growing embedding, `P(t_n)/T_n^h=1+O(T_n^-1)`. At each
other embedding the factor tends to `sigma(P(c))`, with error
`O(rho^n)`. Every limiting factor is nonzero, since `P(c)!=0` is an
algebraic element (`X(c)` and `Y` are real, with `Y!=0`). Hence

```text
A_n/T_n^h = B+O(T_n^-delta),
B=Norm(P(c))/P(c) !=0.                              (3)
```

Put

```text
eta_n=A_n/bar(A_n),       eta=B/bar(B).
```

These lie on the unit circle, and (3) gives
`|eta_n-eta|<=A T_n^-delta` for a fixed positive `A` and all large
`n`. Crucially,

```text
eta not in Q(i).                                    (4)
```

Indeed, `Norm(P(c))` is in `Q(i)`. If `B/bar(B)` were in `Q(i)`,
then `P(c)/bar(P(c))` would also be in `Q(i)`. Writing this last
ratio as `z`, we have `z!=1` and

```text
U(c)=i(z+1)/(z-1) in Q(i) intersection K=Q,
```

contrary to section 1. This proof also works with nonreal embeddings
of `K`: the full norm and its conjugate still belong to `Q(i)`.

Consequently no fixed rational rotation can make `eta_n` tend to
one. Removing common Gaussian content from circle points changes no
pair ratio and cannot repair this obstruction.

## 3. A quantitative lemma for any irrational algebraic phase limit

The following height argument does not require the norm construction.
For a Gaussian rational number `z=p/q` in coprime Gaussian integers,
define

```text
H(z)=max(|p|,|q|).
```

Associates do not change this height. It is at least one;
`H(zw)<=H(z)H(w)` and `H(z^-1)=H(z)`. On the unit circle,
`|p|=|q|=H(z)`.

Fix positive constants `A,delta,C`. Suppose a fixed algebraic `eta`
lies on the unit circle outside `Q(i)`, with

```text
D=[Q(i,eta):Q(i)]>=2,
eta_n in Q(i),       |eta_n|=1,
|eta_n-eta|<=A T_n^-delta,       T_n>=1,       T_n -> infinity.
```

If rational unit-circle rotations `u_n` make
`gamma_n=u_n eta_n` the ratio `z_1/z_0` of two nonzero Gaussian
integer points with `|z_0|=|z_1|=R_n`, contained in an arc of length
`C sqrt(R_n)` on that origin-centered circle, then

```text
H(u_n) >= c T_n^[delta/(D(2D+1))]                   (5)
```

for a fixed `c>0` depending only on `eta,A,delta,C`. The estimate
uses the actual radius after any common-content removal, in both the
chord bound and the denominator bound below.

To prove it, clear denominators in the minimal polynomial of `eta`
over `Q(i)` to obtain `f in Z[i][t]` of degree `D`. For a rational
unit-circle `z=p/q`, `q^D f(p/q)` is a
nonzero Gaussian integer. Also
`|f(z)-f(eta)|<=L|z-eta|` on the unit circle, where
`L=sum_j j|f_j|>0`. Therefore

```text
|z-eta| >= L^-1 H(z)^-D.                            (6)
```

Applying this to `eta_n` yields
`H(eta_n)>=c_0 T_n^(delta/D)`. Write `H=H(u_n)`.
Because a reduced denominator of a point ratio divides one of
the Gaussian integer points, `H(gamma_n)<=R_n`. Thus

```text
c_0 T_n^(delta/D) <= H(eta_n) <= H R_n,
R_n >= c_0 T_n^(delta/D)/H.                          (7)
```

The chord bound and unit moduli now give

```text
|u_n-eta^-1|
 <= |eta_n-eta|+|gamma_n-1|
 <= A T_n^-delta+C R_n^-1/2
 <= A T_n^-delta+C c_0^-1/2 H^1/2 T_n^-delta/(2D).
```

Use (6) for `eta^-1`, which generates the same field over `Q(i)`
as `eta` and therefore has the same degree `D`. Since
`H>=1` and `T_n>=1`, the first term is also at most
`A H^1/2 T_n^-delta/(2D)`. Hence

```text
c_1 H^-D <= c_2 H^1/2 T_n^-delta/(2D),
H^(D+1/2) >= (c_1/c_2) T_n^(delta/(2D)).
```

This is exactly (5).

## 4. Application and scope

Apply section 3 to the phase limit in (3)-(4). Its relative degree
`D` is at most `[K:Q]`. In fact `Q(i,eta)=Q(i,U(c))`, by the
same fractional-linear relation used above. Every correction that
rotates this norm row into an endpoint configuration must therefore
pay a fixed positive power of `epsilon^n` in Gaussian height.

In particular, corrections of height `exp(o(n))` cannot work, even
if the primitive circle has large common-content cancellation.
If the rotation is `u_n=K_n/bar(K_n)` for a Gaussian integer
correction factor, then `H(u_n)<=|K_n|`, so (5) also lower-bounds
that correcting factor. With fixed rotations, (7) already makes
the primitive radius grow exponentially, while the limiting phase
stays away from one.

For a concrete field, take `alpha^3=2` and
`epsilon=1+alpha+alpha^2=(alpha-1)^-1`. This is an algebraic integer
unit of norm one. Its distinguished value exceeds one; its two
other conjugates are a complex-conjugate pair, each of modulus
`epsilon^-1/2`. Thus `delta=1/2`. The phase degree is `D=3`, since
section 1 excludes degree one and the field is cubic. In this case
(5) reads

```text
H(u_n) >= c epsilon^(n/42).                         (8)
```

The hypotheses are nonempty. The earlier explicit quartic

```text
P(t)=(t^2+1)(t^2+4)
     +alpha[-t(t^2+3)+2i]+alpha^2(t^2+1)
```

is monic, vanishes at `i`, and has imaginary polynomial norm `-64`.
At `c=0`, `Norm(P(0))=68-64i`. No splitting assertion for this
quartic is needed for the phase theorem.

The constant-imaginary polynomial norm hypothesis is essential to
the present proof of (4); it is not necessary for every possible
unit-parameter construction. Algebraic offsets need a separate
condition, as section 5 shows. Several growing embeddings and parameters
that vary their field or unit with `n` are also outside the result.
It supplies no lower bound for all integer profiles extracted from
arbitrary short arcs. The uniform endpoint theorem and the general
growth improvement remain open.

## 5. Algebraic offsets can have rational limiting phase

For any fixed offset `xi in K`, the asymptotic formula (3) still
holds, now with `B=Norm(P(xi))/P(xi)`. Each contracting embedding
tends to its corresponding `sigma(xi)`. The same fractional-linear
argument gives the exact criterion

```text
B/bar(B) in Q(i)  if and only if  U(xi) in Q.         (9)
```

Unlike section 1, this criterion can hold. Here is an explicit example;
there is no assertion of linear splitting or a full Boolean profile.
Take the positive real `alpha` with `alpha^3=1/2`, and define

```text
g(t)=(4t^4-2t^3+10t^2-7t+6)/5,
h(t)=(34t^4-12t^3+20t^2-37t-14)/100,
b=2(g^2+1),
a=g+b h,
c=1+4g h+2b h^2,
U=a+b alpha+c alpha^2.
```

Direct multiplication gives `a^2+1=bc/2`. The second elementary
symmetric polynomial in the coefficient conjugates is consequently

```text
e_2(U)=3(a^2-bc/2)=-3.
```

Since `g(i)=-i` and `h(i)=-i/4`, one has `U(i)=-i`. The degrees
of `a,b,c` are `12,8,16`, with leading coefficient of `c` equal to
`4624/15625`. Set

```text
Y=15625/(4624 alpha^2),       P=Y(U+i).
```

Then `P` is monic of degree 16, has constant nonzero imaginary part
`Y`, vanishes at `i`, and satisfies
`Im Norm(P)=(e_2(U)-1) Norm(Y)=-4 Norm(Y)`.

At the algebraic offset `xi=alpha`, exact reduction by `alpha^3=1/2`
gives

```text
g(alpha)=1-alpha+2alpha^2,
g(alpha)^2+1=5alpha^2,
h(alpha)=-alpha*g(alpha)/5,
a(alpha)=0,       b(alpha)=5/alpha,
c(alpha)=1/(5alpha^2),       U(alpha)=26/5.
```

This does not contradict the coefficientwise polynomial identity
`e_2(U)=-3`: conjugating the evaluated element also moves the
algebraic argument, so evaluation no longer commutes with `e_2`.

The field has the Pisot unit
`epsilon=1+2alpha+2alpha^2=(2alpha^2-1)^-1`.
Indeed `beta=2alpha^2` has `beta^3=2`, and this unit is
`1+beta+beta^2`, the one in section 4. For the actual parameter
path `t_n=alpha+epsilon^n`, the limiting phase is rational:

```text
eta=((26+5i)/(26-5i))^2 in Q(i).
```

To check this last formula, put `q=26/5`. Since `P(alpha)=Y(q+i)`,
its norm is `Norm(Y)(q+i)^3`, so `B` is a nonzero real multiple of
`(q+i)^2`. A fixed rational rotation can therefore align this limit
to one. This does not prove that the resulting phases are close enough
relative to their primitive circle radii to meet the endpoint scale.
It does prove that the irrational-limit obstruction cannot cover all
algebraic offsets under the existing hypotheses.

The [exact checker](check_norm_descent_rational_phase_gap.py) checks
2,196 odd-power valuation instances across both parity cases, the
displayed cubic norm and unit identities, and 32 Pell fixtures for
Gaussian height, cancellation, and algebraic phase approximation.
It also verifies the degree-16 offset example by exact polynomial and
cubic-field arithmetic, including its constant imaginary norm.
These finite checks supplement the valuation and height proofs;
they do not establish their universal quantifiers by enumeration.
