# Euclidean realization of the full-Boolean norm and bracket data

The actual Gaussian rows carry more information than their rank-two bracket
array. This note gives an exact criterion for the norm data to come from one
Euclidean plane, proves an integral lift when one prescribed norm is a
Gaussian norm, and clears every forced Boolean factor from its first triple
equation. The
criterion is a different coordinate form of the positive-quadratic condition
in [the many-row norm interpolation note](many_row_norm_barycentric_system.md);
it does not by itself prove a radius-uniform arc bound.

## 1. A necessary and sufficient metric criterion

Fix `m>=3`. Let `N_i` be positive rational numbers and let
`Delta_ij=-Delta_ji` be nonzero rational numbers satisfying all rank-two
Pluecker relations

```text
Delta_ij Delta_kl-Delta_ik Delta_jl+Delta_il Delta_jk=0.
```

Put `D=Delta_12` and represent this determinant array by the rational
vectors

```text
v_1=(1,0),       v_2=(0,D),
v_k=(-Delta_2k/D, Delta_1k),       k>=3.             (1)
```

The Pluecker relations say exactly that `det(v_i,v_j)=Delta_ij` for
every pair. For `k>=3` define

```text
C_k=N_1 Delta_2k^2+N_2 Delta_1k^2-N_k D^2,
S=C_3/(2 Delta_13 Delta_23).                       (2)
```

There is a unique rational symmetric form matching the first three norms:

```text
Q=[[N_1, S/D], [S/D, N_2/D^2]].                    (3)
```

It matches every `N_k` if and only if the `m-3` cleared moment identities

```text
C_k Delta_13 Delta_23=C_3 Delta_1k Delta_2k,
                                             4<=k<=m, (4)
```

hold. Its determinant equals one if and only if

```text
C_3^2=4 Delta_13^2 Delta_23^2 (N_1 N_2-D^2).     (5)
```

Since `N_1>0`, (5) makes `Q` positive definite. Thus (4)-(5) are
necessary and sufficient for **real** vectors `P_i=(X_i,Y_i)` with

```text
X_i^2+Y_i^2=N_i,       det(P_i,P_j)=Delta_ij.
```

Indeed, choose a real matrix `A` with `A^t A=Q` and `det A=1`, then
take `P_i=A v_i`. Conversely any such points pull the ordinary dot
product back to (3). Equation (5) is the scale-one Euclidean condition;
mere positivity of `Q` would permit an arbitrary area multiplier.

If all `N_i,Delta_ij` are integers and (4)-(5) hold, every dot product
`S_ij=v_i^t Q v_j` is an **integer**. It is rational by (1),(3), and

```text
S_ij^2=N_i N_j-Delta_ij^2 in Z;                 (6)
```

a rational number whose square is an integer is an integer.

## 2. Integer realization follows from one Gaussian norm

Under the integral hypotheses just stated, assume additionally that
`N_1` is the norm of a Gaussian integer. Then (4)-(5) are sufficient
for **integer** vectors `P_i in Z[i]` with all the prescribed norms and
determinants. No separate rational-rotation obstruction remains.

Set `Z_ij=S_ij+i Delta_ij in Z[i]`. Real Euclidean realization gives
the exact rank-one identities

```text
Z_ij Z_jk=N_j Z_ik,
bar(Z_1j) Z_1k=N_1 Z_jk.                            (7)
```

Include `Z_11=N_1` and let `H=gcd_(1<=j<=m) Z_1j` in `Z[i]`.
We claim that `H` has a divisor
`G` of norm `N_1`. At a split rational prime `p`, write
`p=pi bar(pi)` and `p^e || N_1`. Put

```text
a=min_(1<=j<=m) v_pi(Z_1j),
b=min_(1<=j<=m) v_bar(pi)(Z_1j).
```

The second identity in (7), with `j,k` allowed to equal 1, implies
`a+b>=e`. Choose an integer
`r` with `0<=r<=a` and `0<=e-r<=b`, and put
`pi^r bar(pi)^(e-r)` in `G`. At an inert prime `p`, the norm
assumption makes `v_p(N_1)=2e`. Since
`Norm(Z_1j)=N_1 N_j`, the Gaussian `p`-valuation of `Z_1j` is
`e+v_p(N_j)/2>=e`; thus `p^e | Z_1j` for all `j`.
At the ramified prime `rho=1+i`, if `2^e || N_1`, then
`v_rho(Z_1j)=e+v_2(N_j)>=e`, so `rho^e | Z_1j` for all `j`.
These choices give `G|H` and
`Norm(G)=N_1`.

Now define

```text
P_1=bar(G),       P_j=Z_1j/G  (j>=2).               (8)
```

They are Gaussian integers. Their norms are `N_j`, and (7) gives
`bar(P_j)P_k=Z_jk`, so every prescribed determinant is attained.
Conversely every integer realization appears in (8) with
`G=bar(P_1)` (up to a common Gaussian unit). Therefore the finite set
of divisors `G|H` with norm `N_1` is also an exact test for *extra*
requirements: conjugate-primitivity, a chosen orientation of each core
factor, positivity and order of real coordinates, and small `|Y_i|`.
Those extra requirements do not follow from (4)-(5). A Gaussian integer
realizing an odd split-prime norm can contain both primes above that
rational prime, so integral realization alone does not enforce
conjugate-primitivity.

For disjoint odd split-core profiles, these primitive-orientation
requirements can nevertheless be recovered after only subpower
factor losses; see [the primitive completion theorem](boolean_metric_primitive_completion.md).
All fixed-metric integral rotations are classified in
[the common Gaussian rotation note](gaussian_common_rotation_classification.md).
Neither result supplies the simultaneous small imaginary coordinates.

The norm assumption is necessary. For example, `N_1=3`, `N_2=3`,
`N_3=6` and all three positive brackets equal to 3 satisfy (4)-(5),
with `Q=diag(3,1/3)`, but 3 is not a Gaussian norm. This example has
a real realization and no Gaussian integer realization.

Even odd norms supported on split primes do not make every integral
realization primitive. Take the actual integer rows
`P_1=5`, `P_2=3+4i`, `P_3=2+i`. Their norms are `(25,25,5)` and
their pair determinants are `(Delta_12,Delta_13,Delta_23)=(20,5,-5)`.
They satisfy (4)-(5), but no realization with these data can make
both first rows conjugate-primitive. Indeed
`Z_12=15+20i=(2+i)^3(2-i)`. A primitive Gaussian integer of norm
25 is associated to `(2+i)^2` or `(2-i)^2`; the product
`bar(P_1)P_2` would then have `(2+i)`-valuation `0`, `2`, or `4`,
never `3`. The divisor test in (8) detects this obstruction exactly.
The [primitive-completion theorem](boolean_metric_primitive_completion.md)
shows that, for the disjoint odd split-core profile, a common rational
rotation followed by row-content removal can recover a primitive oriented
profile after only subpower core loss. That theorem does not restore the
small imaginary coordinates relative to the integer real axis.

## 3. The first triple equation after all forced core cancellation

Now use the full-Boolean arithmetic data from actual primitive rows:

```text
P_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T),       k_i=Norm(K_i),
G_ij=product_(T containing i,j) n_T,
Delta_ij=t_ij G_ij,       N_i=k_i product_(T containing i) n_T.
```

For the triple `{1,2,3}`, partition the blocks by their intersection
with this triple:

```text
U_i=product_(T: T intersect {1,2,3}={i}) n_T,
V_ij=product_(T: T intersect {1,2,3}={i,j}) n_T,
F=V_12 V_13 V_23
  (product_(T: {1,2,3} subset T) n_T)^3.           (9)
```

Empty products are one. Checking the eight intersection patterns
individually gives

```text
C_3=F (k_1 t_23^2 U_1 V_23
       +k_2 t_13^2 U_2 V_13
       -k_3 t_12^2 U_3 V_12),
G_12 G_13 G_23=F,
N_1 N_2/G_12^2=k_1 k_2 U_1 U_2 V_13 V_23.        (10)
```

By (6), `S_12` is an integer. Since `N_1N_2` is divisible by
`G_12^2` and `Delta_12` by `G_12`, equation (6) shows that
`G_12 | S_12`. Put `s_12=S_12/G_12 in Z`. The metric conditions
for this triple become the exact *integer* equations

```text
k_1 t_23^2 U_1 V_23+k_2 t_13^2 U_2 V_13
    -k_3 t_12^2 U_3 V_12=2 t_13 t_23 s_12,       (11)
s_12^2+t_12^2=k_1 k_2 U_1 U_2 V_13 V_23.       (12)
```

Eliminating `s_12` is the fully cleared determinant-one equation:

```text
(k_1 t_23^2 U_1 V_23+k_2 t_13^2 U_2 V_13
       -k_3 t_12^2 U_3 V_12)^2
 +4 t_12^2 t_13^2 t_23^2
 =4 k_1 k_2 t_13^2 t_23^2 U_1 U_2 V_13 V_23.  (13)
```

The integrality of `s_12` in (11) is already supplied by the actual
metric and does not have to be postulated as an independent condition.
For each other `k`, the analog of (11) must give the **same** `s_12`;
these are the `m-3` moment compatibilities (4) in normalized form.
Once (4) and the single determinant-one equation (5) hold, the same
positive form `Q` gives every other triple's Heron equation
automatically. Those equations are not independent extra constraints.

When every block has `log n_T=w+o(w)` and
`log k_i,log|t_ij|=o(w)`, each of the six private factors `U_i,V_ij`
has logarithmic size `2^(m-3)w+o(w)`. Each of the three large terms
in (11), and `|s_12|`, has size `exp(2^(m-2)w+o(w))` on the positive
near-real arc. The main terms in (13) have exponent
`2^(m-1)w+o(w)`, whereas its additive correction
`4t_12^2t_13^2t_23^2` is `exp(o(w))`. The forced common factor `F`
removed from `C_3` has exponent `3*2^(m-2)w+o(w)`.

Equation (13) is a substantial necessary near-square constraint beyond
Pluecker rank two or the isolated lower norm circuit. It is also exactly
the determinant-one part of the full positive-quadratic/barycentric
system already identified in the many-row norm interpolation note, in
different coordinates. Its very small additive correction does not yet
force a strict height gap: the giant square terms can cancel because
they come from one genuine Euclidean triple. A radius-uniform argument
would need to use the moment compatibilities (4) with the one
determinant-one equation (13), the allowed core orientations in (8),
and the small imaginary parts. Repeating (13) for more triples adds no
independent condition after (4)-(5).

There is a direct Gaussian lift explaining that cancellation. Let `A_i`
be the product of the actual `H_T` with
`T intersect {1,2,3}={i}`, and let `B_ij` be the product with
intersection `{i,j}`. Then `Norm(A_i)=U_i` and `Norm(B_ij)=V_ij`.
The ordinary rank-two circuit
`Delta_23 P_1+Delta_31 P_2+Delta_12 P_3=0`, divided by its common
Gaussian factors, becomes

```text
W_1+W_2+W_3=0,
W_1=t_23 K_1 A_1 bar(B_23),
W_2=t_31 K_2 A_2 bar(B_13),
W_3=t_12 K_3 A_3 bar(B_12).                         (14)
```

Their norms are the three large summands in (11), and their oriented
area is `det(W_1,W_2)=t_12 t_23 t_31`. To see the area, write each
`W_i=Delta_jk P_i/D_0` with the same Gaussian denominator
`D_0=B_12 B_13 B_23 C Norm(C)`, where `C` is the product of the
triple-intersection blocks. Then `Norm(D_0)=F`, and
`Delta_12 Delta_23 Delta_31/F=t_12 t_23 t_31`.
The law of cosines and the area identity for this Gaussian triangle
are exactly (13). Thus (13) exposes a narrow integral triangle with
large, nearly balanced, private-factor side norms; by itself it is a
repackaging of the actual three-row Gaussian circuit, not a new
height descent.

The [exact checker](check_boolean_norm_metric_realization_criterion.py)
verifies the metric conditions and constructive integer lift on rational
forms involving split, inert, and ramified primes; the primitive
obstruction fixture; and the full fifteen-block normalization and area
identity for four rows. The arguments above do not rely on the finite
fixtures.
