# An exact quadratic-determinant reduction for all fifteen blocks

This note concerns the full fifteen simple blocks on four nonanchor
rows, including every nonempty subset of `{1,2,3,4}` once. A row
has degree eight, each pair gcd has degree four, each triple gcd
has degree two, the common gcd has degree one, and the lcm has
degree fifteen. All roots are assumed nonreal, distinct, and
globally disjoint from their complex conjugates.

The actual constant-imaginary polynomial equations reduce to four
oriented triple roots, four real scalar parameters, and six explicit
quadratic divisibility tests. This is a necessary-and-sufficient
reduction with nondegeneracy conditions retained. It is not a proof
that the system has, or lacks, a solution, and it gives no automatic
transfer to an arbitrary integer endpoint sequence.

## 1. Normalize the common root and name the factors

Let the normalized rows be `F_i+i`, where `F_i in R[t]` have
degree eight, every pair difference has degree eight, and

```text
F_i-F_j divides F_i F_j+1.                                 (1)
```

A real affine change of variable sends the common root of all
four rows to `-i`; its slope may be negative. Then

```text
F_i(t)=(t^2+1)A_i(t)+t,       A_i in R[t], deg A_i=6,
F_i+i=(t+i)[(t-i)A_i+1].                                  (2)
```

Let `z_j` be the root in the triple block missing row `j`,
and let `e_ij=e_ji` be the root in the pair block `{i,j}`.
Define their real monic norm quadratics by

```text
Q_j=(t-z_j)(t-bar z_j),
Q_ij=(t-e_ij)(t-bar e_ij).                                 (3)
```

The eleven norm quadratics `t^2+1`, the four `Q_j`, and the
six `Q_ij` are pairwise coprime and have negative discriminant.
The pair gcd identity gives real nonzero constants `c_ij`, with
`c_ji=-c_ij`, such that

```text
A_i-A_j=E_ij product_(k notin {i,j}) Q_k,
E_ij=c_ij Q_ij.                                           (4)
```

In particular every `E_ij` has degree two. Its leading coefficient
can have either sign; `E_ij/leading(E_ij)` is the positive monic
quadratic in (3).

## 2. All six pair quadratics are two-by-two determinants

Telescoping (4) around a triple and cancelling its common fourth
quadratic gives

```text
E_ij Q_k+E_jk Q_i=E_ik Q_j.                               (5)
```

Consequently there exist real polynomials

```text
P_i(t)=c t+d_i,                                           (6)
```

with one common coefficient `c`, for which

```text
E_ij=P_j Q_i-P_i Q_j.                                     (7)
```

Here is a proof including the degree restriction. Take row 4 as
reference. Equation (5) says

```text
E_i4 Q_j-E_j4 Q_i=E_ij Q_4.
```

Because the `Q_i` are coprime, the residues `E_i4/Q_i` in
`R[t]/(Q_4)` are all equal. Represent this residue by the unique
real polynomial `P_4` of degree at most one. The polynomials

```text
P_i=(P_4 Q_i-E_i4)/Q_4
```

are then real and have degree at most one. Their coefficients of
`t` all equal that of `P_4`, because the quadratics are monic
and `deg E_i4=2`. Substitution proves (7), with its displayed
sign. This also covers `c=0`, when all the `P_i` are constant.

Put `R_i=P_i product_(k!=i) Q_k`. Equations (4) and (7) give
`A_i-A_j=R_j-R_i`. Hence there is a single real polynomial `K`
such that

```text
A_i=K-P_i product_(k!=i) Q_k.                             (8)
```

It has degree at most seven, and `[t^7]K=c`. This cancels the
common degree-seven term on the right of (8), leaving degree at
most six as required. No unsupported cancellation of lower
coefficients has been used.

## 3. The four oriented triple roots determine K uniquely

At a triple root `z_j`, every row `i!=j` vanishes in (2).
The product subtracted in (8) contains `Q_j`, so

```text
K(z_j)=-1/(z_j-i),       1<=j<=4.                         (9)
```

Conversely (9) gives the required triple-root values in every
incident row. The roots `z_j,bar z_j` are eight distinct nodes,
and none is `i` or `-i`. Interpolate (9) and its conjugate
values at these eight nodes. There is a unique real polynomial
`K` of degree at most seven with those values. Thus, once the
four oriented triple roots have been chosen, both `K` and

```text
c=[t^7]K
```

are fixed. Only the four real constants `d_i` in (6) remain
free. In particular `K` is not an additional moving polynomial
that can be tuned independently of the triple factors.

For a useful factorized form, put

```text
Z=product_j(t-z_j),
M=((t-i)K+1)/Z.                                          (10)
```

The numerator vanishes at all four `z_j`, so `M in C[t]` has
degree at most four and leading degree-four coefficient `c`.
Its conjugate satisfies the exact Bezout identity

```text
(t+i) Z M-(t-i) bar Z bar M=2i.                           (11)
```

Define the four residual polynomials

```text
W_i=(t-z_i)M-(t-i)P_i product_(j!=i)(t-bar z_j).            (12)
```

The two nominal degree-five coefficients in (12) are both `c`
and cancel. Thus `deg W_i<=4`, and direct substitution gives

```text
F_i+i=(t+i) product_(j!=i)(t-z_j) W_i.                    (13)
```

The three pair roots incident to row `i` and its private root
must be exactly the four roots of `W_i`.

## 4. Six actual-value tests, with a converse

Choose four nonreal oriented roots `z_i` and four real scalars
`d_i`. Assume the triple roots are distinct, conjugate-disjoint,
and avoid `i,-i`. Construct `Q_i,K,c,P_i,E_ij,M,W_i` by
(3), (6)--(12). The remaining tests are as follows.

1. Each `E_ij` has degree two and negative discriminant. Its
   monic normalization, all the other normalized `E_kl`, the
   four `Q_i`, and `t^2+1` are pairwise coprime.
2. Each `W_i` has degree exactly four.
3. For every pair, perform the actual polynomial division

```text
E_ij divides W_i bar W_i       in R[t].                    (14)
```

   There are six such conditions, each a remainder of degree
   at most one. One endpoint `i` per unordered pair suffices.
4. The edge and private roots recovered below are all distinct
   and globally conjugate-disjoint, including with `-i` and
   all the `z_i`.

These tests are necessary and sufficient for the desired full
fifteen-block model.

For sufficiency, form `F_i` by (2) and (8). It is real, and
(13) with test 2 makes it a quartic residual times four known
linear factors, hence a polynomial of degree eight. Its imaginary
part after adding `i` is exactly the constant one.

Fix a pair `i,j`. By (7)--(8), the real difference `F_i-F_j`
contains the factor `E_ij`. By (14), one root `e_ij` of this
irreducible real quadratic is a root of `W_i`. Thus
`F_i(e_ij)=-i`, and the difference identity forces
`F_j(e_ij)=-i` as well. Test 1 ensures this root is not one
of the known full or triple roots in (13), so it is a root of
`W_j`. The other root, `bar e_ij`, cannot also be a root of
`F_i+i`, because a real `F_i` takes the conjugate value `+i`
there. Thus exactly one orientation of this edge is selected,
coherently in both incident rows.

The three incident edge quadratics for a row are coprime, giving
three distinct roots of its degree-four `W_i`. The remaining
root is its private block. Test 4 is the explicit final check
that all fifteen roots give precisely the prescribed incidences.
The exact difference factorization and its degree imply (1).

For necessity, (13) in an existing full-profile system has its
three incident edge roots inside `W_i`; their real norms divide
`W_i bar W_i`, proving (14). All other tests follow from the
specified simple, independent fifteen-block profile.

Condition (14) is an actual small-value compatibility condition.
It does not follow merely from the determinant identities (5)
or (7), which are already built into the parameterization.

## 5. Arithmetic scope and the remaining system

If the four `z_i` have rational real and imaginary parts and the
four `d_i` are rational, then `K,P_i,E_ij,F_i` are rational
polynomials and `M,W_i` have coefficients in `Q(i)`. Whenever
the tests in Section 4 pass, every recovered edge root is also
Gaussian-rational. Indeed it is the unique complex common root
of `W_i,W_j`, so their degree-one gcd over `Q(i)` determines
it. A quadratic minimal polynomial over `Q(i)` would instead
force both conjugate roots into `F_i+i`, which is impossible.
The private roots are then Gaussian-rational as the remaining
degree-one quotients. Thus rational solutions of this reduced
system would give the rational linear-block model required for
an arithmetic construction.

The four triple roots provide eight real coordinates and the
`d_i` provide four. The six tests (14) are twelve explicit real
remainder equations, with the strict and nonvanishing conditions
listed separately. These numbers describe an exact finite system;
they are not used as an equation-count proof of impossibility or
existence. No full-fifteen-block solution, or contradiction to
this system, is asserted here.

The root orchestrator independently audited the proof and supplied
[an exact identity certificate](check_full_fifteen_reduction.py).
It checks the real interpolation, Bezout identity, degree
cancellations, all row and pair factorizations on 24 arbitrary
rational triple inputs, and two reflection slices. It deliberately
does not certify the remaining divisibility tests (14), their
solvability, or a full-fifteen-block example.

## 6. A reflection-symmetric slice retaining every support

Unlike the eight-block construction, the following symmetry does
not remove any of the fifteen supports. Reflect roots by
`J(z)=-bar z` while permuting rows by `(12)(34)`. Write the
four triple roots as

```text
z_1=a+ib, z_2=-a+ib, z_3=d+ie, z_4=-d+ie,                 (15)
```

with real coordinates satisfying the nondegeneracy conditions
above. Uniqueness in (9) makes `K` odd. Indeed the interpolation
values transform by

```text
K(-bar z_j)=-bar(K(z_j)).
```

Take

```text
P_1=ct+p, P_2=ct-p, P_3=ct+q, P_4=ct-q.                  (16)
```

Equations (2) and (8) then give

```text
F_2(t)=-F_1(-t),       F_4(t)=-F_3(-t).                   (17)
```

The invariant pair quadratics have particularly simple formulas:

```text
E_12=-2[(2ac+p)t^2+p(a^2+b^2)],
E_34=-2[(2dc+q)t^2+q(d^2+e^2)].                           (18)
```

If their leading coefficients are nonzero and

```text
rho_12=p(a^2+b^2)/(2ac+p)>0,
rho_34=q(d^2+e^2)/(2dc+q)>0,
```

their roots are `+/- i sqrt(rho_12)` and
`+/- i sqrt(rho_34)`. The even parts of `F_1` and `F_3`
vanish at the respective roots by the difference identities.
Their values there are therefore purely imaginary, and each
condition `F_i(eta)^2+1=0` is just one real equation.

For pairs `13` and `14`, choose either root `eta` of their
negative-discriminant quadratics and impose
`F_1(eta)^2+1=0`, two real equations per pair. Reflection gives
the other two mixed pairs. Thus the slice has a concrete system
of six real equations in `a,b,d,e,p,q`, followed by all the
original nondegeneracy checks. No cut has been discarded: the
full root and the two invariant pair roots lie on the imaginary
axis, and the other twelve roots form reflection pairs.

An initial bounded floating-point Newton experiment on this
system produced only near-degenerate approximations, with triple
roots approaching one another or the common root. No exact
solution or rigorous exclusion was certified. The usable result
is the exact reduction (15)--(18), not numerical evidence for
nonexistence.
