# The first Gram lift and its exact dual quotient

The at-most-two theorem in
[coupled_orientation_frame_packing.md](coupled_orientation_frame_packing.md)
counts metrics for fixed rational source vectors and contact norms.
It does not exclude the first metric. This note gives an exact
first-lift-dependent description of the whole integral Gram lattice,
and audits one natural construction from the source vectors. The
construction does not currently produce an existence obstruction.

## 1. The congruence lattice after one lift is known

Fix primitive integer contact vectors `lambda_e` and pairwise coprime
odd split rational contact norms `d_e`, and put `D=product d_e`. Set

```text
Gamma={G=[[a,b],[b,c]] in Sym_2(Z):
       lambda_e^t G lambda_e=0 mod d_e for every e}.
```

Each evaluation functional is primitive, so `Gamma` has index `D` in
`Sym_2(Z)`. Suppose one determinant-one frame `U=(z w) in SL_2(Z)`
exists, and write `F=U^t U`. Choose each Gaussian prime-power
orientation from this frame and let their product be `Delta`, so
`Norm(Delta)=D`. No prime-power exponents are discarded.

Define the integer linear map

```text
Psi_F(G)=a w^2-2bzw+c z^2 in Z[i],
T_F(G)=tr(F^(-1)G) in Z.                              (1)
```

Then the following statements are exact:

```text
G in Gamma  <=> Delta divides Psi_F(G),               (2)
ker(Psi_F:Sym_2(Z)->Z[i])=Z F,                         (3)
Psi_F(Sym_2(Z))=Z+2iZ,                                (4)
T_F(G)^2-Norm(Psi_F(G))=4 det(G).                      (5)
```

For (2), work at a chosen `pi^e` over `p^e || d_e`.
Write `lambda_e=(r,s)`. If `s` is a unit modulo `p`, then
`rz+sw=0 mod pi^e` gives

```text
Psi_F(G)=z^2 s^(-2)(a r^2+2brs+c s^2) mod pi^e.
```

Here `z` is a unit modulo `pi`: otherwise both frame coordinates
would vanish, contrary to determinant one. If `r` is a unit, use the
analogous formula with `w`. The multiplier is invertible. An ordinary
integer is divisible by `pi^e` exactly when it is divisible by `p^e`,
which proves both directions of (2).

For the other claims, the integral unimodular congruence

```text
H=U^(-t) G U^(-1)=[[alpha,beta],[beta,gamma]]
```

is an automorphism of `Sym_2(Z)`. Direct expansion gives

```text
Psi_F(G)=-(alpha-gamma+2i beta),
T_F(G)=alpha+gamma.
```

This proves (3)--(5). In particular there is an exact sequence

```text
0 -> Z F -> Gamma -> (Delta) intersect (Z+2iZ) -> 0.   (6)
```

More explicitly, `Gamma` is in bijection with pairs `(T,xi)` satisfying

```text
xi in (Delta) intersect (Z+2iZ),
T in Z,    T=Re(xi) mod 2.                            (7)
```

The reconstruction is

```text
H=[[ (T-Re xi)/2, -Im xi/2 ],
   [ -Im xi/2, (T+Re xi)/2 ]],
G=U^t H U.                                           (8)
```

Thus its positive determinant-one forms are exactly the pairs in (7)
with `T>0` and `T^2-Norm(xi)=4`. Every positive integral binary form
of determinant one is an integral Gram form: choose a primitive
shortest vector, complete it to an integral unimodular basis, and shear
the second basis vector to obtain `|b|<=a/2` and `c>=a`. Then `1=ac-b^2>=3a^2/4` forces `a=1`, `b=0`,
`c=1`. Consequently these forms correspond to all admissible frames
up to common Gaussian units, not merely to a larger relaxation. Indeed
a row of an integral unimodular frame is primitive; if its norm is
divisible by `p^e`, exactly one Gaussian orientation carries the whole
exponent. Thus a Gram solution of the rational contact congruences
supplies an admissible orientation at every prescribed prime power.

The dependence on the first frame is essential. Its unknown matrix
`U` occurs in both the coordinate change and the physical height
constraint `tr(G)<=H_bound^2`. The abstract description (6) therefore
cannot be used to assert that a short first lift exists or fails to
exist.

## 2. A precise sufficient obstruction from a second integral form

If `G in Gamma` is not a rational multiple of `F`, then (2)--(3) imply

```text
|Psi_F(G)|>=sqrt(D).
```

Cauchy--Schwarz applied to the coefficient vectors
`(a,sqrt(2)b,c)` and `(w^2,-sqrt(2)zw,z^2)` gives

```text
|Psi_F(G)| <= ||G||_F (|z|^2+|w|^2)
             =||G||_F ||U||_F^2.
```

Hence every independent integral form in the same congruence lattice
satisfies the deterministic bound

```text
||G||_F >= sqrt(D)/||U||_F^2.                          (9)
```

This is an ideal-divisibility bound, not a covolume estimate. To exclude
a first frame with `||U||_F<=H`, it suffices to construct a nonzero
`G in Gamma`, known not to be proportional to any positive definite
form, with

```text
||G||_F < sqrt(D)/H^2.                                (10)
```

A rank-one or indefinite `G` would meet the independence condition.
For the full paired cut data, including private norms, and
`A=2^(m-3)`, the profile values are

```text
log D=(8A-2)w+o(w),    log H=Aw+o(w).
```

The threshold in (10) is therefore `exp((2A-1)w+o(w))`.
A valid use of this criterion needs uniform error bounds and a strict
height margin; matching leading exponents alone is not enough.

## 3. What products of source linear forms actually cost

Let `v_i=(r_i,s_i)` be the common primitive source vectors, and put

```text
L_i(X,Y)=det(v_i,(X,Y))=r_i Y-s_i X.
```

For each paired nonconstant cut, let its median collision cluster be
`C_e`. Its representative `lambda_e` is projectively equal modulo
`d_e` to every source row in that cluster. Thus, whenever `i in C_e`,

```text
d_e divides L_i(lambda_e).
```

For distinct labels `i,j`, define the complete missing-contact factor

```text
M_ij=product_(C_e disjoint from {i,j}) d_e.
```

The integral symmetric form corresponding to

```text
G_ij(X,Y)=2 M_ij L_i(X,Y)L_j(X,Y)                     (11)
```

belongs to `Gamma`: every contact either divides `M_ij` or divides
one of its linear factors on evaluation. The factor two makes the
middle matrix coefficient integral. This construction retains every
prescribed prime-power depth and needs no clean-outsider assumption.
It gives an indefinite form for distinct projective source rows, with

```text
det(G_ij)=-M_ij^2 det(v_i,v_j)^2,
Psi_F(G_ij)=2 M_ij Q_i Q_j,    Q_i=Uv_i.               (12)
```

There are `4A-1` paired nonconstant cuts. An anchor belongs to `A`
median clusters and an outside label to `2A`. Two anchors belong to
no common cluster; an anchor and an outside label belong to `A/2`
common clusters; two outside labels belong to `A` common clusters.
For the last case there must be two outside labels, so `m>=5`.
Each prescribed rational contact norm has logarithm `2w+o(w)`.
The resulting height audit is

| Labels | Missing clusters | Source linear-form height sum | `log ||G_ij||_F` |
|---|---:|---:|---:|
| Two anchors | `2A-1` | `o(w)` | `(4A-2)w+o(w)` |
| Anchor and outside | `3A/2-1` | `Aw+o(w)` | `(4A-2)w+o(w)` |
| Two outside | `A-1` | `2Aw+o(w)` | `(4A-2)w+o(w)` |

Here the first three source vectors have subexponential height and
each outside vector has height `exp(Aw+o(w))`. The norm of the
symmetric product of two nonzero real linear forms is within absolute
constant factors of the product of their norms, so coefficient
cancellation does not invalidate the displayed height exponents.

Every choice in (11) is therefore larger than the sufficient threshold
(10) by `exp((2A-1)w+o(w))`. In the coordinates (1), formula (12)
likewise has logarithmic modulus `(6A-2)w+o(w)`, whereas the first
nonzero ideal element only has the lower bound `(4A-1)w+o(w)`.
The gap is the same. Moving the selected labels from anchors to
outside rows saves precisely as many missing contact factors as it
costs in source-vector height.

If only nonprivate contacts are fixed, omit the `m` singleton clusters.
Two selected labels cover exactly two of those private contacts, so the
missing multiplier loses `m-2` factors in every row of the table. Thus

```text
log ||G_ij||_F=(4A-2m+2)w+o(w),
log(sqrt(D)/H^2)=(2A-m-1)w+o(w).                      (13)
```

The product construction exceeds its sufficient obstruction threshold
by `exp((2A-m+3)w+o(w))`. This distinguishes the full-private and
shared-contact calculations; their contact sets cannot be interchanged
inside either height estimate.

This calculation does not rule out specially chosen combinations of
these forms with unusually large common divisors or cancellations.
It identifies the next required condition precisely: construct an
indefinite or rank-one element of `Gamma` below (10), using additional
incidence or residue information. Neither the subset counts nor the
obvious product construction establishes that condition.

## 4. Verification

`python3 docs/check_primal_dual_gram_contact_audit.py` checks the exact
linear map, its integral reconstruction and parity, the determinant
identity, full-depth contact equivalence, product-form identities, and
the cluster counts for `3<=m<=10`. It does not assume or certify a
uniform first-lift obstruction.
