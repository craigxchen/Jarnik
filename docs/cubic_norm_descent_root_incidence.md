# Cubic norm descent: single rows and root incidence

Let `K` be a real cubic number field, `L=K(i)`, and let `q=e_2` be
the second elementary symmetric function of the three `K/Q` embeddings.
The field need not be totally real or Galois. For

```text
P=X+iY,     U=X/Y,     X in K[t],     0 != Y in K,
```

one has

```text
Im Norm_(L/Q(i))(P) = Norm(Y) (q(U)-1).                 (1)
```

If `U(i)=-i`, the fixed-root calculation in
[the parity note](fixed_gaussian_root_norm_descent_parity.md) shows that
a descended imaginary degree at most one forces that imaginary part to
be constant. Thus the bounded-content cubic endpoint condition is
exactly `q(U)=-3` for every row. This note proves that this condition is
attainable by a single monic quartic over **every** real cubic field,
and gives additional restrictions when all input roots split over `L`.
It supplies neither a full four-row construction nor an exclusion of
all such constructions.

## 1. A quartic over every real cubic field

Write `Tr=Tr_(K/Q)`. For `z in K`,

```text
q(z) = ((Tr z)^2-Tr(z^2))/2.
```

The trace-zero plane is nondegenerate for this quadratic form. Choose
a nonzero `Y_0` on that plane with `q(Y_0)!=0`, and set

```text
x=1/Y_0,       w=3x/Tr(x)-1.
```

Indeed, `q(x)=Tr(Y_0)/Norm(Y_0)=0`, while
`Tr(x)=q(Y_0)/Norm(Y_0)!=0`. Consequently `Tr(w)=0` and `q(w)=-3`.
Choose a nonzero trace-zero `v` orthogonal to `w` for the polar form
of `q`, and put `delta=-q(v)/3`. Nondegeneracy and `q(w)!=0` give
`delta in Q\{0}`. The elements `1,w,v` form a rational basis of `K`,
and

```text
q(a+bw+cv)=3(a^2-b^2-delta c^2).                       (2)
```

Choose nonzero rational `u,v_0` satisfying
`v_0^2+delta u^2=1`, and put `k=v_0/(4 delta u^2)`. Such choices exist:
take

```text
u=2r/(1+delta r^2),       v_0=(1-delta r^2)/(1+delta r^2)
```

for any rational `r` avoiding the finitely many zero numerators and
denominators. Define rational polynomials

```text
s=t^2+1,
a=-(t^3+3t)/2+k s^2,
b=s [k t^2-t/2+k(1-2 delta u^2)],
c=s [2ku t-(delta u^2+1)/(2 delta u)].                 (3)
```

For an explicit verification put `h=delta u^2`. Polynomial expansion gives

```text
a^2+1-b^2-delta c^2
  =-s^2 (h-1)/(4h) [16k^2 h^2+h-1]=0,
```

because `16k^2 h^2=v_0^2`. In particular,

```text
a^2+1=b^2+delta c^2.
```

Hence `U=a+bw+cv` satisfies `q(U)=-3` and `U(i)=-i`. Its leading
coefficient is the nonzero element `k(1+w)`. Set

```text
Y=1/[k(1+w)],       P=Y(U+i).
```

Then `P` is monic of degree four, has nonzero constant imaginary
part `Y`, vanishes at `i`, and satisfies

```text
Im Norm(P)=-4 Norm(Y) != 0.                            (4)
```

The leading coefficient `ell=k(1+w)` obeys `q(ell)=0`, so
`Tr(Y)=q(ell)/Norm(ell)=0`. Since `Y!=0`, it cannot be rational and
therefore generates the cubic field `K`. The construction does not
hide a proper enlargement of the row's coefficient field.

The nondegeneracy used above follows directly from the trace pairing:
on the trace-zero plane, `q(z)=-Tr(z^2)/2`, and that plane is the
orthogonal complement of `1`, whose trace-pairing norm is `3`.
Thus the construction also applies to totally real cubic fields. It
does not assert that the quartic splits into linear factors over `L`.
The nonsplitting example in the parity note shows why this additional
requirement cannot be omitted.

## 2. The common rational factor has small degree

Suppose now that `P` has degree `n`, all its roots are simple,
`U(i)=-i`, and `q(U)=-3`. In a normal closure over `Q(i)`, write the three conjugates
of `p=U+i` as `p_1,p_2,p_3`. Equation (2) is not needed here. Direct
expansion gives

```text
e_2(p_1,p_2,p_3)=2i(p_1+p_2+p_3).                    (5)
```

Let `G` be the monic gcd of these three polynomials. Galois invariance
puts `G` in `Q(i)[t]`. It is squarefree because `P` has simple roots.
Equation (5) implies

```text
G^2 divides Tr(U)+3i.
```

Coefficient conjugation gives `bar(G)^2 | Tr(U)-3i`. The polynomials
`G` and `bar(G)` are coprime, since they respectively divide `U+i`
and `U-i`. Differentiating and using squarefreeness therefore gives

```text
G bar(G) divides (Tr(U))'.
```

The polynomial `Tr(U)` is nonconstant: its value at `i` is `-3i`,
whereas its coefficients are rational and real. Consequently

```text
2 deg G <= deg Tr(U)-1 <= n-1.                        (6)
```

For a degree-eight row, `deg G<=3`. Since the fixed root `i` belongs
to `G`, no complete orbit of three non-Gaussian roots can also belong
to the same row.

There is also a useful multiplicity restriction. If two conjugate
polynomials `p_j` vanish at a point, (5) forces the third to vanish
there. Thus pairwise gcds of `p_1,p_2,p_3` equal their common gcd.
If two distinct members of one nonbase `Q(i)`-orbit occur among the
roots of a row that splits over `L`, their degree-three minimal
polynomial has all three roots in `L`: a cubic extension containing
two of its roots contains the third and is normal. Applying (5) then
forces the whole orbit into that row. In degree eight, (6) and the
fixed root exclude this possibility. Each row therefore contains
at most one root from each nonbase orbit, and at most three roots
in `Q(i)`.

## 3. Consequences for the full four-row profile

Assume the full profile of
[the four-row certificate](four_row_real_polynomial_profile_certificate.md),
with all thirty oriented roots distinct, lies in `L`, and every row
obeys `q(U_j)=-3`. The following are necessary conditions.

* If `L/Q(i)` is non-Galois, all thirty oriented roots automatically
  have distinct `Q(i)`-Galois orbits. A nonbase element generates the
  cubic field; a second conjugate in that field would induce a
  nontrivial automorphism, which a non-Galois cubic does not have.
  The bounded-content conclusion of
  [the content theorem](real_algebraic_norm_descent_content.md) therefore
  applies without an additional orbit-separation hypothesis.
* In the exceptional case where `L/Q(i)` is cyclic, any two different
  input labels in the same oriented orbit must have disjoint cuts.
  Otherwise they occur together in one row, contradicting section 2.
  A collision involving opposite orientations also requires disjoint
  cuts: an intersecting row would make its descended norm and its
  conjugate share a root, although their difference is the nonzero
  constant `-8i Norm(Y_j)`.
* At most seven of the fifteen input roots lie in `Q(i)`. Each row
  contains at most three such roots, so their total row incidence is
  at most twelve. The fixed mask `15` uses four incidences. Among the
  other labels only four are singletons, and every additional label
  uses at least two incidences; the remaining budget of eight permits
  at most six labels. In particular each row has at least five roots
  outside `Q(i)`.

For clarity, normality of `L/Q(i)` occurs when `K/Q` is cyclic, or
when `K/Q` is non-Galois with quadratic resolvent `Q(i)`. In the latter
case its normal closure equals `K(i)`. Indeed, the splitting field of
a non-Galois cubic is obtained by adjoining the square root of its
discriminant, and the only quadratic subfield of `K(i)` is `Q(i)`.

Disjoint-cut collisions are not ruled out by these arguments. In
particular, two labels with opposite orientations in an orbit pair
can have bounded common content. The conclusions above must not be
strengthened to a general distinct-orbit assertion in the cyclic case.
No four-row endpoint exclusion follows from this note.

The dependency-free
[checker](check_cubic_norm_descent_root_incidence.py) verifies (3) by
exact polynomial arithmetic for sixteen signed rational choices of
`delta`, the fixed-root and leading-coefficient identities, and the
finite Boolean incidence count used above. These checks supplement
the general algebraic proofs; they do not test root splitting or
construct a full profile.
