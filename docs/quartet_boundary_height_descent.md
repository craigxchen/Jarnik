# A coefficient-height descent for quartet-supported boundary relations

A large space of numerical relations can have
[compatible boundary basepoints](coherent_boundary_basepoints_of_numeric_relations.md)
at every first smaller cut. For the subspace induced by
numerical zero relations on one fixed quartet, however, its nonzero
restriction directions cannot be short. This statement uses coefficient
heights, including cancellation against the full deeper kernel. It does
not classify arbitrary short relations or prove the uniform arc bound.

Write `m=2q>=10`, fix a quartet `A`, and put `B={1,...,m}\A`.
Use the degree-two-each invariant spaces and kernels `K_1,K_2` of
[the all-cut hierarchy](all_cut_invariant_relation_hierarchy.md).
Let `Lambda_A` be the rational space of quartet invariants that vanish
at the actual four directions. Define the rational polynomial subspace

```text
E_A = K_2(m) + span_Q { Q_A P : Q_A in Lambda_A, P in K_1(B) }.
```

Here `K_1(B)` is the balanced kernel for the `m-4` complement rows.
No bound on the coefficients of the displayed generators is assumed.
An element of `E_A` can have cancellations between many products and
an arbitrary element of `K_2(m)`.

**Height descent.** In the primitive full-profile Gaussian system, let
all core blocks satisfy `log N(H_T)>=(1-eta)w`, and let the corrections
satisfy `log|K_i|<=sigma`. If an integer polynomial `W in E_A` is not
in `K_2(m)`, its polynomial coefficient sum obeys

```text
log C_W >= 2^(m-3)(1-eta)w - 4 sum_(i in A) log|K_i|
        >= 2^(m-3)(1-eta)w - 16 sigma.                 (1)
```

Consequently all elements of this subspace below that height lie in
the deeper kernel. The assertion includes the subspace with just one
fixed quartet relation, but does not require a common polynomial factor.

## 1. Integral quartet coordinates

Temporarily label the quartet `1,2,3,4`, and set

```text
X=Delta_12 Delta_34,       Z=Delta_14 Delta_23.
```

The integral quartet invariant lattice of degree two in each row has
basis `X^2,XZ,Z^2`. This follows directly from integral bracket
straightening and `Delta_13 Delta_24=X+Z`; alternatively these are its
three standard monomials.

For an inside pair `I`, evaluate its vectors at `(0,1)` and the other
two at `(1,0)`. Complementary pairs have equal values. On the three
representatives `12,13,14`, the balanced evaluation matrix in this
basis is

```text
          X^2   XZ   Z^2
   12      0     0    1
   13      1    -1    1
   14      1     0    0 .                              (2)
```

Its determinant is `1`. Thus every integral triple `(a,b,c)` is the
balanced value vector of the unique integral polynomial

```text
Q=c X^2+(a+c-b)XZ+a Z^2.                              (3)
```

In particular, an integral vector in the rational image of `Lambda_A`
comes from an *integral* numerical-zero quartet relation. This removes
any denominator or primitive-content ambiguity in the next step.

## 2. Extracting a quartet relation from a short polynomial

First, `E_A` is contained in `K_1(m)`. Indeed, at a global balanced
cut, either at least three quartet rows are collapsed, making every
quartet invariant zero by weight balance, or at least `q-2` complement
rows are collapsed, making the factor in `K_1(B)` zero.

Now take a first smaller inside set `S`, so `|S|=q-1`, and put
`r=|S intersect A|`. If `r>=3`, the quartet factor vanishes. If `r<=1`,
at least `q-2` complement rows are collapsed and the `K_1(B)` factor
vanishes. Only `r=2` can give a nonzero restriction. Write

```text
S=I union J,       I subset A, |I|=2,
J subset B, |J|=q-3.
```

The quartet restriction is then its balanced scalar, independent of
the remaining two quartet variables. For a product the exact identity is

```text
F_(I union J,Q_A P)=c_I(Q_A) F_(J,P).                  (4)
```

Every `K_2(m)` summand has zero restriction. Fix `J` and any quadratic
monomial `z^alpha` in the outside complement variables. For an integer
`W in E_A`, form its three coefficient values

```text
a_I = [z^alpha] F_(I union J,W),       I=12,13,14.     (5)
```

They are integers, each of absolute value at most `C_W`, because the
restriction substitutes only zeros, ones and individual variables.
By (4), their vector is a rational linear combination of balanced
vectors of elements of `Lambda_A`. Equation (3) therefore reconstructs
an integer numerical-zero quartet invariant `Q_alpha` with balanced
values exactly `(a_12,a_13,a_14)`.

If `W` is not in `K_2(m)`, some first smaller restriction is nonzero.
The preceding case split then supplies a choice of `J,alpha` for which
this reconstructed quartet relation is nonzero. At least one of its
balanced values is nonzero and has absolute value at most `C_W`.
This extraction controls the coefficients after all cancellations in
`W`; it does not estimate the heights of a chosen product decomposition.

## 3. Aggregating all conductor blocks over the quartet

Use the actual primitive central numerators

```text
P_i=X_i+iY_i=K_i product_(T containing i)H_T,
gcd_G(P_i,bar(P_i))=1.
```

The blocks have disjoint odd split-prime supports, also disjoint from
their conjugates. For each pair `I subset A`, define

```text
H_I^A=product_(T intersect A=I) H_T.
```

There are exactly `2^(m-4)` factors: the membership of each complement
label is free. These global subsets are all nonempty and proper, so no
empty- or full-set convention affects this count.

Apply the exact [balanced scalar congruence](small_invariant_relations_at_balanced_cuts.md)
to any integer quartet numerical zero `Q`. At each constituent block,
the two inside imaginary coordinates are units. The outside core
factors are units, leaving only the outside corrections. Since the
constituent blocks are coprime, they combine before any height loss:

```text
H_I^A | c_I(Q) product_(j in A\I) K_j^2.              (6)
```

Thus the positive integer

```text
M_I = N(H_I^A) /
      N(gcd_G(H_I^A,product_(j in A\I) K_j^2))
```

divides `c_I(Q)`. The complementary pair has the same scalar and a
coprime norm modulus. If `c_I(Q)!=0`, this yields

```text
log|c_I(Q)| >= log N(H_I^A)+log N(H_(A\I)^A)
                    -4 sum_(j in A)log|K_j|
             >= 2^(m-3)(1-eta)w
                    -4 sum_(j in A)log|K_j|.          (7)
```

The correction loss is charged once to each aggregated block. It is
not multiplied by the number of constituent subsets. Every step keeps
the full Gaussian prime powers; no coprimality between corrections and
core blocks is presumed.

Apply (7) to the nonzero balanced value extracted in Section 2. Its
absolute value is at most `C_W`, proving (1).

## Scope and verification

For fixed `m`, the height in (1) is far above the Chow-form short-relation
threshold. Therefore a quartet-supported boundary construction can
satisfy all the geometric common-basepoint conditions while contributing
no nonzero first-smaller restriction at those short heights. This is a
positive height obstruction to that specific construction.

The argument does not show that arbitrary short numerical relations
belong to any one of the spaces `E_A`, and it does not justify replacing
one fixed quartet by a sum over all quartets. It also leaves short
relations in `K_2(m)` uncontrolled. These are substantive remaining
gaps for a uniform theorem.

The exact finite certificates are in
[check_quartet_boundary_height_descent.py](check_quartet_boundary_height_descent.py).
They check the unimodular quartet map and mixed-product coefficient
arrays at actual rational directions on `m=10,12`, including a numerical
zero term in `K_2` at twelve rows. All `210` and `792` first-smaller cuts
are checked, with `60` and `120` nonzero extracted arrays. Separate
Gaussian block examples at `m=4,6` check prime-power aggregation and
shared-prime corrections; those examples test the local divisibility
mechanism, not a common-circle endpoint configuration.
The general proof is the integral reconstruction (2)--(5) and the
aggregated prime-power divisibility (6)--(7), not an inference from
the finite examples.
