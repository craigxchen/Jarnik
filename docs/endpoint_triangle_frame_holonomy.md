# Endpoint triangle frames, holonomy, and short directed circuits

The actual four-row Boolean profile has integral triangle frames whose
pair determinants are the original small residuals. Their transition
matrices have determinant one and small denominators. The holonomy is an
explicit unipotent shear, however, and neither its trace nor its small
rational shear parameter bounds its ambient height. The triangle
relations themselves clear to the existing five-star bracket relations.

There is also a precise limitation on an additive approach using pair
quotients: two distinct directed words, and every three-word pattern
except a directed triangle, have a singleton valuation obstruction.
A primitive, unramified Pell family realizes a bounded-coefficient
directed triangle with bounded imaginary coordinates. It is a grouped
three-block example, not a realization of the full four-row profile.

## 1. The exact integral frames

Use the conventions of
[the integer five-star lift](four_row_integer_five_star_lift.md). Let
`I={1,2,3,4}` and

```text
P_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T),
G_ij=product_(T containing i,j) n_T,
t_ij=Im(bar(P_i) P_j)/G_ij,
t_ji=-t_ij.
```

The fifteen core blocks are odd, conjugate-primitive Gaussian integers
with pairwise coprime norms; the `K_i` are nonzero Gaussian integers and
may overlap core primes. Assume every `t_ij` is nonzero. Define

```text
V_i=K_i H_{i} bar(H_(I\{i})),
F_ij=H_{ij},                          1<=i<j<=4,
L_i^ell=V_i F_iell bar(F_jk),          {i,j,k,ell}=I.
```

For each omitted label `ell`, the three `L_i^ell`, `i!=ell`, are
Gaussian integers. Direct cancellation of shared core factors gives

```text
bar(L_i^ell) L_j^ell=bar(P_i) P_j/G_ij=:Z_ij.          (1)
det(L_i^ell,L_j^ell)=Im(bar(L_i^ell)L_j^ell)=t_ij.     (2)
```

Here Gaussian integers are identified with their two real coordinate
columns. In particular, the right side of (1) is independent of the
choice of the omitted label. If a pair quotient is instead defined as
`P_i bar(P_j)/G_ij`, it is `bar(Z_ij)` and its imaginary part is
`-t_ij`; no sign convention is being suppressed.

To check (1), the exclusive factors of `P_i` relative to `P_j` are
`H_i,H_ik,H_iell,H_ikell`; the exclusive factors of `P_j` are
`H_j,H_jk,H_jell,H_jkell`. These are exactly the conjugated and
unconjugated factors in `bar(L_i^ell)L_j^ell`, with the same correcting
factor `bar(K_i)K_j`.

Under the balanced endpoint hypotheses

```text
log n_T=w+o(w),
log max(1,|K_i|), log max(1,|t_ij|)=o(w),
```

each frame column has modulus `exp(2w+o(w))`, while its pairwise
determinants have height `exp(o(w))`. Its columns need not individually
lie near the real axis; their relative angles are small modulo pi.

The elementary triangle relation is

```text
t_jk L_i^ell-t_ik L_j^ell+t_ij L_k^ell=0.              (3)
```

It is not an extra independent relation on the stars. In the notation
of the five-star lift,

```text
Q_i=V_i F_ij F_ik F_iell,
L_i^ell=Q_i bar(F_jk)/(F_ij F_ik).
```

Multiplying (3) by `F_ij F_ik F_jk` gives exactly

```text
t_jk n_jk Q_i-t_ik n_ik Q_j+t_ij n_ij Q_k=0.          (4)
```

## 2. Transition matrices and exact holonomy

For two omitted labels `a,b`, define `T_ab` as the real linear map
from frame `b` to frame `a` that sends the two shared labeled columns
to their counterparts. If those labels are `i,j`, then

```text
T_ab=[L_i^a L_j^a] [L_i^b L_j^b]^(-1).
```

Equation (2) implies `det T_ab=1`. Its entries are rational with
denominators dividing the nonzero integer `t_ij`. This follows directly
from the integer adjugate formula, including when correction factors
overlap core primes.

Consider the ordered loop of omitted labels `4 -> 3 -> 2 -> 4` and put

```text
H=T_42 T_23 T_34,
v_1=L_1^4, v_2=L_2^4,
Pf(t)=t_12 t_34-t_13 t_24+t_14 t_23.
```

Then the exact holonomy formula is

```text
H v_1=v_1,
H v_2=v_2-Pf(t)/(t_13 t_14) v_1.                     (5)
```

For a proof write `v_i=L_i^4`, `u_i=L_i^3`, and `w_i=L_i^2`.
The transition maps fix the label `1` along the loop and respectively
identify the other shared labels `2`, `4`, and `3`. The triangle
determinants give

```text
v_3=-(t_23/t_12)v_1+(t_13/t_12)v_2,
u_2=(t_24/t_14)u_1+(t_12/t_14)u_4,
w_4=-(t_34/t_13)w_1+(t_14/t_13)w_3.
```

Apply the three transitions successively to the second equation and
then substitute the first and third. The coefficient of `v_2` is one;
the coefficient of `v_1` is
`(t_13 t_24-t_12 t_34-t_14 t_23)/(t_13 t_14)`, proving (5).

Thus `tr H=2` holds identically, and this loop is flat if and only if
`Pf(t)=0`. Flatness is an additional condition, not an automatic
consequence of defining the four frames. The actual weighted Pluecker
identity is

```text
n_12 n_34 t_12 t_34
 -n_13 n_24 t_13 t_24
 +n_14 n_23 t_14 t_23=0,                              (6)
```

as in [the full Boolean circuits](full_boolean_small_imaginary_circuits.md).
Its three generally distinct matching norm factors do not cancel to
the unweighted expression in (5).

The height issue can be made exact. With
`lambda=-Pf(t)/(t_13 t_14)`, equation (5) gives

```text
H-I=(lambda/t_12) v_1 (-Im(v_1),Re(v_1)),
||H-I||_op=|lambda| |v_1|^2/|t_12|.                   (7)
```

If `Pf(t)!=0` and all residuals are subpower, both upper and lower
logarithmic bounds give `log |lambda|=o(w)`. Consequently
`||H-I||_op=exp(4w+o(w))` in the balanced profile. The small rational
coefficient of the shear in the almost parallel frame therefore
corresponds to a large ambient matrix. Each individual transition also
has the elementary upper bound `exp(4w+o(w))`; a small denominator
does not give small entries. Equations (5)--(7) provide no positive
residual-height exponent.

## 3. The singleton obstruction for at most three directed words

For distinct labels `i,j`, let the bare directed word be

```text
E_ij=product_(T containing i but not j) H_T
      product_(T containing j but not i) bar(H_T).
```

Thus `bar(Z_ij)=K_i bar(K_j) E_ij`. For a nonempty proper subset
`S` of labels, `H_S` divides `E_ij` precisely when the directed edge
`i -> j` exits the cut `S`. If it does not exit, `E_ij` is coprime
to `H_S`; it may contain `bar(H_S)`, which is also coprime to `H_S`.

Suppose distinct directed words satisfy

```text
sum_(e in E) c_e E_e=0,
```

with all coefficients nonzero Gaussian integers. If some directed cut
contains all but one of these edges, reduction modulo `H_S` forces
`H_S | c_e` for the unique remaining edge. This keeps the full block
valuation depth, without assumptions about coefficient prime support.
In particular that coefficient has modulus at least `sqrt(n_S)`.
For a relation among the actual quotients the coefficient here is the
original coefficient multiplied by `K_i bar(K_j)`. Subpower correcting
factors therefore do not remove the positive block-height cost.

The following elementary classification is exact for a full Boolean
family of cut blocks:

* Any two distinct directed edges admit a cut containing exactly one.
* Three distinct directed edges admit a cut containing exactly two
  unless they form a directed three-cycle.
* A directed three-cycle has at most one edge exiting any cut, so the
  singleton test does not apply to its three-term relation.

For the second assertion, two edges can exit the same cut exactly when
neither edge's head is the other edge's tail. If some compatible pair
exists, choose such a cut. It contains two or three of the edges. In
the latter case all three edges form a bipartite graph across the cut.
A simple bipartite graph with three edges has a degree-one vertex;
moving that vertex across the cut leaves exactly two exiting edges.
If no compatible pair exists, all pairs are head-to-tail composable.
An opposite pair cannot be supplemented by a distinct third edge
incompatible with both; absent an opposite pair, pairwise
incompatibility forces the directed three-cycle. This also proves the
claimed exception. The two-edge assertion follows in the same way:
if a cut contains both edges, their two-edge bipartite graph has a
degree-one vertex, which can be moved to leave exactly one.

This classifies the singleton test for two and three words; it is not
a classification of larger relations. For example, six edges oriented
by a total ordering of four labels have no cut containing five edges,
since a cut on four vertices contains at most four edges.

## 4. A primitive bounded-coefficient directed cycle

Let

```text
u+b sqrt(13)=(649+180 sqrt(13))^n,     n>=1,
A=(u+3b)+2ib,
B=2b+i(u-3b),
C=(2u+8b)+i(u+b)=2A+B.
```

Then `u^2-13b^2=1`, `u` is odd, and `b` is divisible by six. The
recurrence preserves these facts. The coordinate columns of `A,B`
have determinant

```text
(u+3b)(u-3b)-4b^2=1.
```

Thus `A,B` are primitive integer vectors, and so is `C`, the image
of the primitive vector `(2,1)` under this unimodular matrix. Their
coordinate parities are respectively odd/even, even/odd, and even/odd.
They are consequently odd and conjugate-primitive Gaussian integers.

Their norms are

```text
n_A=1+26b^2+6ub,
n_B=1+26b^2-6ub,
n_C=5+130b^2+34ub.                                    (8)
```

These are pairwise coprime. Indeed

```text
n_A-n_B=12ub,
n_C-5n_A=4ub,
n_C-5n_B=64ub.
```

A shared prime is odd. For the first pair, the exceptional prime `3`
is excluded by `b=0 mod 3`, which gives `n_A=n_B=1 mod 3`.
Every other possible shared prime must divide `ub`. A prime dividing
`b` gives `n_A=n_B=1`; a prime dividing `u` gives
`13b^2=-1` and hence `n_A=n_B=-1` modulo that prime. These observations
also exclude all possible common primes in the second and third
pairs. The three norm logarithms are `2 log b+O(1)`, since
`u/b -> sqrt(13)` and each leading constant in (8) is positive.

The directed products are exactly

```text
W_12=A bar(B)=4ub-i,
W_23=B bar(C)=1+26b^2+2ub+2i,
W_31=C bar(A)=2+52b^2+16ub+i,
3W_12+2W_23-W_31=0.                                  (9)
```

All three have positive real parts and bounded nonzero imaginary
coordinates. They are conjugate-primitive: each is a product of two
conjugate-primitive factors with coprime norms. Thus (9) has no
ramified-prime or common-content loophole. Its three large factors
have balanced logarithmic height, while its coefficients are fixed.

The family realizes the exceptional directed-cycle pattern in a
three-factor model. It does not provide the fifteen independent
balanced blocks of an actual four-row endpoint profile, nor the
remaining row and pair compatibility conditions. Conversely, the
singleton test, trace-two holonomy, and these triangle-frame identities
alone do not establish the radius-uniform arc count.

## 5. Exact checks

The independent
[checker](check_endpoint_triangle_frame_holonomy.py) verifies 300
nonsingular fifteen-block Gaussian fixtures, including correction
factors overlapping core factors. It checks all six pair products,
all twelve rational determinant-one transitions, and the holonomy
matrix. It also checks the two-loop traces, the canonical gauges and
closed-word height bounds, the seven norm equations, and all four
additive gcd identities below. All 300 selected fixtures have nonzero
residual Pfaffian. The separate isotropic-line counterfixture is
checked modulo 17, and all 66 directed edge pairs and 220 edge triples
are enumerated, leaving exactly the eight directed three-cycles.
The checker also verifies 32 powers in the primitive Pell family. These finite
checks supplement the proofs above; the general identities do not
depend on the sampled block sizes or choices.

## 6. A second isotropic line is not automatically invariant

At a split core prime `p` not dividing `t_12 t_13 t_14`, reduce (7)
modulo `p`. For any vector `z` one obtains

```text
det(z,Hz)=-(lambda/t_12) det(v_1,z)^2.                 (10)
```

If `Norm(v_1)` is a unit and `z` is nonzero and isotropic, then
`det(v_1,z)!=0`: otherwise `v_1` and `z` would span the same isotropic
line. Consequently preservation of the line of `z` is equivalent to
`Pf(t)=0 mod p`. The assignment of an isotropic column to each frame
does not itself establish this extra preservation.

Here is a concrete counterfixture with all corrections equal to one.
Assign the following Gaussian coordinate pairs to masks `1,...,15`,
where bit `i-1` denotes label `i`:

```text
(1,2), (2,3), (1,4), (2,5), (1,6), (4,5), (2,7),
(5,6), (3,8), (5,8), (4,9), (1,10), (3,10), (7,8), (4,11).
```

Their norms are distinct split primes. The six residuals are

```text
(t_12,t_13,t_14,t_23,t_24,t_34)
 =(-2675448,-1081583,-2556377,3637775,1624445,-4799520),
Pf(t)=5298313940220.
```

At `H_12=1+4i`, with `p=17`, these reduce to

```text
(t_12,t_13,t_14,t_23,t_24,t_34)=(12,8,15,13,10,5),
Pf(t)=5,
H=[[16,2],[15,3]],
v_1=(6,6),
z=L_3^4=(2,9),             Hz=(16,6).
```

The vector `z` is isotropic, whereas `Norm(Hz)=3 mod 17`; in
particular `det(z,Hz)=4 mod 17` is nonzero. All relevant transition
denominators are units here. This refutes automatic preservation
based only on the core-prime incidence pattern. The finite fixture
does not satisfy a sequence of endpoint height bounds and does not
exclude a further argument that uses those bounds.

## 7. Traces of combined loops retain only the residual data

In the common basis `(v_1,v_2)` of frame `4`, let `H_1=H` be the
loop `4 -> 3 -> 2 -> 4`, and let `H_2` be the loop
`4 -> 3 -> 1 -> 4`. Then

```text
H_1=[[1,a],[0,1]],       a=-Pf(t)/(t_13 t_14),
H_2=[[1,0],[b,1]],       b= Pf(t)/(t_23 t_24).          (11)
```

The second formula follows from (5) by swapping labels `1,2`, which
negates the Pfaffian. Direct matrix multiplication gives

```text
tr(H_1 H_2)=2-Pf(t)^2/(t_13 t_14 t_23 t_24),
tr(H_1 H_2 H_1^(-1) H_2^(-1))
 =2+Pf(t)^4/(t_13^2 t_14^2 t_23^2 t_24^2).            (12)
```

Large matrix entries in ambient coordinates therefore do not force a
large trace of this commutator. The small pair determinant cancels
the large lengths in the exact change of basis.

There is a uniform algebraic explanation for every fixed closed word
in the transitions. Set

```text
T=max(1,max_(i<j)|t_ij|),
L=product_(i<j)|t_ij|,       D=L^3.
```

For a frame whose three labels are `a<b<c`, construct its canonical
rational columns

```text
q_a=(1,0),
q_b=(0,t_ab),
q_c=(-t_bc/t_ab,t_ac).                                (13)
```

They have exactly the same three oriented determinants as the actual
frame. There is therefore a unique rational determinant-one map
`S_ell` sending its actual columns to (13). The corresponding
canonical transition is

```text
T'_ab=S_a T_ab S_b^(-1).
```

It is determined by the residuals alone. Each canonical column
coordinate has absolute value at most `T` and denominator dividing
`L`. The adjugate formula on the two shared columns shows that each
entry of `T'_ab` has absolute value at most `2T^2` and denominator
dividing `L^2 |t_ij|`, hence dividing `D`. Reverse transitions obey
the same bound; no separate bound on inversion is needed.

It follows that the trace of any closed transition word of length
`k>=1` can be written

```text
tr(word)=N/D^k,      N in Z,
|N| <= (4T^2 D)^k <= 4^k T^(20k).                     (14)
```

Indeed conjugation by the starting frame's `S` preserves the trace,
and the trace expansion of a product of `k` two-by-two matrices has
`2^k` scalar terms. This proof neither assumes nor removes
coprimality between residuals and core factors, and therefore retains
all correction-factor overlaps.

For every fixed `k`, the numerator and denominator in (14) have
logarithmic height `o(w)` when the residuals do. If the word length
grows, the actual estimate is `O(k(1+log T))`; a subpower conclusion
then requires this quantity to be `o(w)`. There is no implicit
uniformity in arbitrary growing word length.

Thus the transition connection, up to a separate rational change of
coordinates at each frame, already depends only on the six residuals.
A new trace or commutator does not impose an additional compatibility
identity on those residuals merely by being formed. A positive gap
from this route would still require a new arithmetic link between a
residual expression such as (12) and the large core factors; that
link is not supplied by the connection or its trace identities.

## 8. Seven norm variables and an exact additive gcd

The frame identities also have a useful integer norm form. Define

```text
a_i=Norm(V_i)=Norm(K_i)n_i n_(I\{i}),
b_1=n_12 n_34,   b_2=n_13 n_24,   b_3=n_14 n_23,
B=b_1 b_2 b_3,
Z_ij=x_ij+i t_ij=bar(P_i)P_j/G_ij.                  (15)
```

The three `b`'s are pairwise coprime. The core parts of the four
`a_i` are also pairwise coprime, but the `a_i` themselves need not be
coprime to one another or to the `b`'s: a correcting factor `K_i` may
contain any core prime. No such coprimality is used below. Extend the
orientation convention by `Z_ji=bar(Z_ij)`, `x_ji=x_ij`, and
`t_ji=-t_ij`.

Let matching `r=1,2,3` respectively mean the edge pairs
`(12,34)`, `(13,24)`, and `(14,23)`. For an edge `ij` in matching `r`,
the exact norm equation is

```text
x_ij^2+t_ij^2=Norm(Z_ij)=a_i a_j B/b_r.           (16)
```

Indeed, if `ell` is a label outside `ij` and the remaining label is
`k`, then `Norm(L_i^ell)=a_i n_iell n_jk`; the corresponding formula
for `j` has the other two matching factors. Each directed triangle
also satisfies

```text
Z_ij Z_jk Z_ki=a_i a_j a_k B>0.                   (17)
```

This is the product of the three column norms of its triangle frame.
For example, when `i<j<k`, its imaginary part gives the entirely
integer cubic identity

```text
t_ij x_jk x_ik+x_ij t_jk x_ik-x_ij x_jk t_ik
 +t_ij t_jk t_ik=0.                                (18)
```

The seven norm variables also recast (6) as

```text
b_1 t_12 t_34-b_2 t_13 t_24+b_3 t_14 t_23=0.     (19)
```

There is an additional exact *additive* way to recover each `a_i`.
For distinct `i,j,k`, put `ell=I\{i,j,k}`. Multiply the triangle
relation (3) by `bar(L_i^ell)` and use (1). Its imaginary part
cancels, leaving

```text
t_ik Z_ij-t_ij Z_ik
 =t_jk Norm(L_i^ell)
 =t_jk a_i n_iell n_jk.                              (20)
```

Thus every quotient in the following gcd is an integer, and its
absolute value is one of `a_i b_1,a_i b_2,a_i b_3`:

```text
gcd_({j,k} subset I\{i})
 |(t_ik x_ij-t_ij x_ik)/t_jk|=a_i.                   (21)
```

Pairwise coprimality of the three `b`'s is the only coprimality used
in (21). This recovers the complete correction norm `Norm(K_i)` along
with the two private core norms; it does not identify their factors.

Equations (16)--(21) are necessary arithmetic constraints on a
four-row realization. They do not by themselves construct the fifteen
oriented Gaussian blocks, split each `b_r` into its two edge norms,
or impose the original small imaginary coordinates of the `P_i`.
The cycle and norm equations specify integral Hermitian Gram data for
each triangle; the existing [norm-metric realization criterion](boolean_norm_metric_realization_criterion.md)
addresses the corresponding Gaussian-integer Gram lift. A common
four-row realization additionally has to make those four triangle
frames agree through the same seven grouped Gaussian factors and,
ultimately, the same fifteen blocks.
