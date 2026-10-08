# Fixed-root norm descent and a quadratic endpoint obstruction

Let `K` be a finite real number field, let `d=[K:Q]`, and put
`L=K(i)`. No Galois or total-reality assumption on `K` is needed.
Suppose

```text
P(t)=X(t)+iY,       X in K[t],       Y in K\{0},
P(r)=0,            r in Q(i).
```

If the relative norm has constant imaginary part,

```text
Q(t)=Norm_(L/Q(i))(P(t))=A(t)+iC,
A in Q[t],         C in Q,
```

then `d` is odd and necessarily

```text
C=(-4)^((d-1)/2) Norm_(K/Q)(Y).                     (1)
```

In particular, `C` is nonzero. The conclusion requires neither
monicity, a prescribed row degree, nor coprimality of any factors.

## 1. Proof by evaluation at the fixed root

Since `K` is real, `[L:Q(i)]=d`. Relative norm commutes with
coefficient conjugation, and evaluation at `r in Q(i)` commutes
with this norm. Consequently

```text
Q(r)=0,
bar(Q)(r)=Norm(bar(P)(r))=Norm(-2iY)=(-2i)^d Norm(Y).
```

Here `bar(Q)(r)` means coefficient conjugation followed by evaluation,
not conjugation of the value `Q(r)`. Subtracting the two polynomial
identities gives

```text
2iC=-(-2i)^d Norm(Y).                               (2)
```

For even `d`, the right side is nonzero and real, whereas the left
side is purely imaginary. For odd `d`, division by `2i` yields (1).
This proves the assertion.

For the normalized root `r=i`, the same evaluation also constrains
a nonconstant imaginary polynomial. Writing `Q=A+iB` with
`A,B in Q[t]`, it gives

```text
B(i)=(-2i)^(d-1) Norm(Y),
B mod (t^2+1) = (-4)^((d-1)/2) Norm(Y)              if d is odd,
B mod (t^2+1) = -2(-4)^((d-2)/2) Norm(Y) t          if d is even.
```

In degree three, an imaginary degree at most one therefore forces
`B` to be constant. Thus the bounded-content cubic endpoint degree
budget still requires exact constant-imaginary preservation; the
positive-degree quadratic example below has no cubic analogue of
degree one.

Every row of the
[certified four-row profile](four_row_real_polynomial_profile_certificate.md)
has the fixed factor `t-i`. Therefore any common real coefficient
field of even degree is excluded from this particular norm-descent
mechanism. The degree of that profile's coefficient field is **not**
known here. This theorem does not exclude its descent in every field,
does not construct a rational profile, and gives no lattice-point
count improvement. An extension of odd degree still has to satisfy
all the constant-imaginary polynomial identities.

## 2. The quadratic coefficient calculation

For `K=Q(sqrt(D))`, write

```text
X=A+sqrt(D) B,       Y=c+sqrt(D) e,
A,B in Q[t],        c,e in Q.
```

Direct multiplication gives

```text
Im Norm(P)=2(c A-D e B).                            (3)
```

If `X` is monic of positive degree `n`, then `A` is monic of degree
`n` and `deg B<n`. Constancy in (3) first forces `c=0`; because
`Y!=0`, it then forces `e!=0` and `B=b` constant. But the fixed root
would give

```text
A(r)+sqrt(D)(b+ie)=0.
```

Linear independence of `1,sqrt(D)` over `Q(i)` forces `b=e=0`, a
contradiction. Thus the usual trace-zero cancellation in the leading
coefficient cannot extend to an exact quadratic norm descent when
there is a fixed Gaussian root.

This excludes a **constant** imaginary part, not the endpoint arc
scale itself. For example, if `A in Q[t]` is monic, nonconstant and
divisible by `t^2+1`, and `s in Q\{0}`, the single row

```text
P=A+sqrt(D)s(t-i)
```

has constant imaginary part and the fixed root `i`, but

```text
Norm(P)=A^2-Ds^2(t-i)^2,       Im Norm(P)=2Ds^2 t.
```

Its surviving imaginary degree is exactly one. In the
[circle-content calculation](real_algebraic_norm_descent_content.md),
the bounded-content quadratic threshold is `q=d/2=1`; positive
imaginary degree can therefore still meet the endpoint degree budget.
No full Boolean profile of this form is constructed here.

## 3. A full quadratic profile cannot keep every imaginary degree at most one

Suppose a full Boolean polynomial profile over `K(i)`, where
`K=Q(sqrt(D))` and `D>0` is nonsquare, has monic rows
`P_i=X_i+iY_i` of degree `n>2`, with `X_i in K[t]` and nonzero
constants `Y_i in K`. Suppose every row has the common root `i`.
Assume the original oriented blocks and their
conjugates are pairwise coprime. As usual, the whole shared norm
product for two rows has degree `n`.

Then at least one descended row has imaginary degree at least two.
In particular this applies to the four-row, fifteen-linear-block
profile if its real coefficient field were quadratic. The theorem
does not assert that this is the coefficient field of the certified
profile.

To prove it, suppose every descended imaginary degree is at most
one. Formula (3), monicity and `n>2` force, in every row,
`c_i=0` and `deg B_i<=1`. The identity at `i` then forces

```text
P_i=A_i+s_i sqrt(D)(t-i),
A_i in Q[t] monic of degree n,     (t^2+1) divides A_i,
s_i in Q\{0}.
```

The original imaginary constants are `Y_i=-s_i sqrt(D)`. They
are pairwise distinct: otherwise a pair bracket would have degree
less than `n`, while its shared norm product of degree `n` divides
it, forcing equal rows and contradicting their distinct coprime
incident factors. Thus `s_i!=s_j`. Their shared norm product is

```text
S_ij=(s_i A_j-s_j A_i)/(s_i-s_j) in Q[t],
deg S_ij=n.
```

But the Gaussian norm of one original row is

```text
N_i=P_i bar(P_i)=F_i+sqrt(D) G_i,
F_i=A_i^2+D s_i^2(t^2+1),       G_i=2s_i t A_i.
```

Every rational polynomial divisor of `N_i` divides both `F_i` and
`G_i`. Their monic gcd is exactly `t^2+1`: one has
`gcd(A_i,F_i)=gcd(A_i,t^2+1)=t^2+1`, and
`F_i(0)=A_i(0)^2+D s_i^2>0` excludes the extra factor `t`.
Equivalently, `gcd(N_i,sigma(N_i))=t^2+1`, where `sigma` is the
nontrivial quadratic automorphism. The rational polynomial `S_ij`
divides `N_i`, but has degree `n>2`, a contradiction.

For the four-row norm circle with distinct oriented Galois orbits,
the [content theorem](real_algebraic_norm_descent_content.md) gives
bounded common content. Its normalized chord exponent is `q_i-1`
when `d=2`, so a row with `q_i>=2` rules out endpoint arcs along
the fixed-denominator limit `t -> +infinity`. If orbit collisions
create growing common content, that exponent also depends on the
content degree; the present argument alone does not exclude those
primitive circle configurations. It proves no general lattice-point
count improvement.

## 4. Whole shared factors retain every pair bracket

Suppose a Boolean profile over `L` has monic rows `P_i` of degree
`n`, and its whole shared factor `S_ij` is monic of degree `n/2`
(as the product of the monic shared Boolean factors):

```text
P_i=S_ij U_i,       P_j=S_ij U_j.
```

Assume both descended rows have constant nonzero imaginary parts,
`Q_i=A_i+iC_i` and `Q_j=A_j+iC_j`. Write
`D_ij=Norm(S_ij)`. Multiplicativity alone gives the polynomial
identity

```text
Im(bar(Q_i) Q_j)
  =D_ij bar(D_ij) Im(bar(Norm(U_i)) Norm(U_j)).       (4)
```

This factors the whole shared product at once. It does not assume
that `D_ij` and `bar(D_ij)` are coprime, or that different descended
Boolean blocks are coprime. Both `D_ij bar(D_ij)` and the possible
nonzero bracket have degree `nd`; monicity identifies the quotient:

```text
C_j A_i-C_i A_j
  =(C_j-C_i) Norm_(K/Q)(S_ij bar(S_ij)).              (5)
```

The norm on the right is the product of the coefficient conjugates
of the entire polynomial `S_ij bar(S_ij) in K[t]`, and equals
`D_ij bar(D_ij)`. If `C_i=C_j`, (4)
instead forces a zero bracket, so `Q_i=Q_j`. Thus equal norm residuals
can merge rows; they cannot be declared distinct without a separate
check.

If both original rows share a root in `Q(i)`, equation (1) substitutes
`C_k=(-4)^((d-1)/2) Norm(Y_k)` into (5), and the common nonzero
factor cancels. This bracket statement does not establish the
coprimality needed by the
[integer specialization theorem](boolean_polynomial_integer_specialization.md).

## 5. Enlarging the row's coefficient field cannot repair descent

If `F` is a proper real subfield of `K` and `P in F(i)[t]`, then,
with `e=[K:F]>1`, the tower formula gives

```text
Norm_(K(i)/Q(i))(P)=Norm_(F(i)/Q(i))(P)^e.
```

For every nonconstant polynomial `B in Q(i)[t]`, the power `B^e`
cannot have a nonzero constant imaginary part. Otherwise
`B^e-bar(B)^e` would be a nonzero constant. Factoring it over the
complex numbers gives

```text
B^e-bar(B)^e=product_(zeta^e=1)(B-zeta bar(B)).
```

Every factor would have to be constant. Two distinct roots of unity
then force both `B` and `bar(B)` to be constant, a contradiction.
The fixed Gaussian root makes the putative imaginary constant
nonzero by (2). Therefore successful fixed-root norm descent requires
the chosen real field to equal the field generated by each row's
coefficients. This statement does not assume that the norm from the
smaller field already has constant imaginary part.

There is a stronger degree statement for monic rows. Put `n=deg P>1`,
`f=[F:Q]`, and `d=ef`. Write the smaller norm as `B=A+iH`, with
`A` monic of degree `nf`. The fixed-root evaluation shows `H!=0`,
and monicity gives `h=deg H<nf`. In the imaginary part of `B^e`,
the term `e A^(e-1) H` has strictly larger degree than every term
containing three or more factors `H`. Thus

```text
deg Im Norm_(K(i)/Q(i))(P)=(e-1)nf+h >= nd/2 > d/2.
```

Consequently a proper field enlargement also fails the bounded-content
endpoint degree threshold in the
[content theorem](real_algebraic_norm_descent_content.md), even when
the smaller norm has nonconstant imaginary part. This corollary retains
that theorem's fixed-denominator parameter and bounded-content hypotheses.

## 6. Odd-degree descent can preserve a nonzero imaginary constant

The parity restriction is sharp for single rows over real fields.
Let `alpha` be the real cube root of two, put `K=Q(alpha)`, and set

```text
a=(t^2+1)(t^2+4),       U=-t(t^2+3),       c=t^2+1,
P=a+alpha(U+2i)+alpha^2 c.
```

Then `P` is monic of degree four, its imaginary part is the nonzero
constant `Y=2alpha`, and `P(i)=0`. The cubic norm formula is

```text
Norm(P)=a^3+2(U+2i)^3+4c^3-6ac(U+2i).
```

Since `ac=U^2+4`, its imaginary part is exactly

```text
12U^2-16-12ac=-64=-4 Norm(Y).                       (6)
```

This fixture has a real coefficient field with two nonreal
embeddings. It does **not** split completely into linear factors
over `K(i)`. To see this, the simple roots `alpha=3` of `x^3-2` and
`i=2` of `x^2+1` modulo five lift by Hensel's lemma and give an
embedding `K(i) -> Q_5`. In this embedding the monic polynomial has
integral coefficients and reduces to

```text
P mod 5=t(t-2)(t^2+4t+2).
```

The last quadratic has discriminant three modulo five, a nonsquare.
If `P` split over `K(i)`, all its roots in `Q_5` would be integral
because `P` is monic with integral coefficients. Reducing its linear
factors would give complete splitting over `F_5`, a contradiction.
Thus it supplies neither a full linear-factor Boolean profile over
this field nor coprime descended blocks. Enlarging the real field
to make its roots split changes its relative norm to a proper power;
section 5 proves that this destroys the nonzero constant imaginary
part. Its purpose is to show that the even-degree exclusion cannot
be asserted in every odd degree under the same single-row hypotheses.

## 7. Derivatives give further necessary conditions

All `d` conjugate factors of `Q` vanish at the fixed root `r`, so
`Q` has a zero there of multiplicity at least `d`. If `Q=A+iC`, then
`Q'=A'` is rational. Since `Y!=0` and `K` is real, `r` is nonreal.
Writing `r=u+iv` with rational `u,v` and `v!=0`, one obtains

```text
((t-u)^2+v^2)^(d-1) divides A'(t).                 (7)
```

For `d>1`, the logarithmic derivative of `bar(Q)` at `r` also gives

```text
Tr_(L/Q(i))(X'(r)/Y)=0.                            (8)
```

Indeed `bar(Q)(r)!=0`, `bar(Q)'(r)=Q'(r)=0`, and each factor of
`bar(Q)` has value `-2i sigma(Y)` at `r`. Conditions (7)--(8) do not
exclude odd-degree descent. No such conclusion is inferred from
them.

The dependency-free
[checker](check_fixed_gaussian_root_norm_descent_parity.py) verifies
the signs in (2), quadratic norm identities on fixed-root examples,
the quadratic endpoint gcd with repeated factors and zero constant terms,
the factorization (4) with deliberate shared and conjugate overlaps,
the odd-degree fixture (6) with its derivative vanishing and
modulo-five nonsplitting certificate, its first five proper powers,
the exact imaginary degrees of one hundred further polynomial powers,
and (5)
on the existing three-row rational family. These finite
checks supplement the general algebraic proofs above.
