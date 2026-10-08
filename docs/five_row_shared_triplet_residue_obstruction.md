# A shared-triplet factorization for nonnodal five-row restrictions

This note records one exact compatibility identity in the unresolved
five-row irreducible case. It keeps the actual primitive pair residues,
the oriented Gaussian cut blocks, and the common complement reflection.
The identity exposes a precise obstruction to turning several generic
boundary congruences into one coefficient modulus. It does **not**
prove the positive individual coefficient-height bound sought in
[five_row_boundary_node_height.md](five_row_boundary_node_height.md),
Section 4.

Use its actual source notation

```text
P_i=K_i product_(U containing i) H_U,
n_U=Norm(H_U),       log|K_i|<=sigma,
0<|t_ij|<=T,         D_ij=Im(P_i bar(P_j))
                        =Norm(gcd_G(P_i,P_j)) t_ij.        (1)
```

The `H_U` and their conjugates have disjoint odd split-prime support.
For one inside pair `S`, write the outside labels as `a,b,c`. At that
boundary line let

```text
F_(S,Q)=alpha(z_a-z_c)+beta(z_b-z_c),
z_i=Y_i/P_i,
alpha beta (alpha+beta) != 0.                         (2)
```

The condition in (2) is the nonnodal conclusion of the boundary-node
note for a sufficiently short irreducible numerical zero.

## 1. One exact Gaussian determinant identity

Multiplying (2) by `P_a P_b P_c` gives its restriction numerator

```text
V_S=alpha D_ac P_b+beta D_bc P_a.                     (3)
```

The three binary vectors `(P_i,Y_i)` have one real determinant
syzygy. With the signed convention in (1), it is

```text
D_ab P_c+D_bc P_a-D_ac P_b=0,
V_S=alpha D_ab P_c+(alpha+beta)D_bc P_a.              (4)
```

Both expressions in (4) keep two nonzero terms under (2). The
syzygy changes a basis of outside differences; it does not make
`V_S` a scalar multiple of one pair bracket. Each pair bracket in
(4) is exactly its common Gaussian norm times the small signed
integer `t_ij`, rather than an arbitrary independent residue.

The actual one-pass numerical-zero congruence is
`H_S | K_a K_b K_c V_S`. The reflected point yields a second
congruence at the complementary core block, with the same labelled
`alpha,beta` and pair residues negated. The two congruences concern
the same invariant section but different normalized Gaussian
numerators; neither says an ordinary integer coefficient in (2) is
divisible by its core modulus.

## 2. Exact core cancellation in the one-pass numerator

For a pair put

```text
G_ij=product_(U containing i,j) n_U,
r_ij=Norm(gcd_G(P_i,P_j))/G_ij in Z_(>0).             (5)
```

The core product divides the Gaussian gcd. Any extra gcd
valuation comes from the corrections, giving
`1<=r_ij<=Norm(K_i K_j)<=exp(4sigma)`.
Thus `D_ij=G_ij r_ij t_ij` exactly, with no discarded prime power.

Classify a cut `U` by its three membership bits on `(a,b,c)`.
The following Gaussian core factors are literal:

```text
C_abc = product_(abc(U) in {110,101,011}) H_U
        product_(abc(U)=111) H_U^2 bar(H_U),
A_abc = product_(abc(U)=010) H_U
        product_(abc(U)=101) bar(H_U),
B_abc = product_(abc(U)=100) H_U
        product_(abc(U)=011) bar(H_U).                 (6)
```

In particular,

```text
G_ac product_(U containing b)H_U=C_abc A_abc,
G_bc product_(U containing a)H_U=C_abc B_abc,
V_S=C_abc [alpha r_ac t_ac K_b A_abc
             +beta r_bc t_bc K_a B_abc].              (7)
```

Each cut type `010,101,100,011` has four subsets, since the two
inside labels are free. Therefore `A_abc` and `B_abc` each contain
exactly eight distinct oriented cut blocks. If their norm logarithms
lie in `[(1-eta)w,(1+eta)w]`, then

```text
log|A_abc|, log|B_abc| in [4(1-eta)w,4(1+eta)w],
log|V_S/C_abc|
 <=log(|alpha|+|beta|)+4(1+eta)w+5sigma+log T. (8)
```

The upper estimate uses `|r_ij K_k|<=exp(5sigma)` and preserves
`|t_ij|<=T`. It remains valid when cancellations make the bracket
much smaller.

There is a sharper compatibility hidden in that cancellation. Put
`U=D_ac P_b/C_abc=r_ac t_ac K_b A_abc` and
`V=D_bc P_a/C_abc=r_bc t_bc K_a B_abc`, both Gaussian integers.
For a cut containing `r` of the three outside labels, the exponent
of `n_U` in `Norm(C_abc)` is `binomial(r,2)`, exactly the number of
outside pairs contained in that cut. Hence

```text
Norm(C_abc)=G_ab G_ac G_bc,
delta_abc=Im(U bar(V))
         =D_ac D_bc D_ba/Norm(C_abc)
         =-r_ab r_ac r_bc t_ab t_ac t_bc,
0<|delta_abc|<=exp(12sigma) T^3.                    (9)
```

The nonzero bound follows from the actual nonzero pair residues.
It also proves `V_S=C_abc(alpha U+beta V)` cannot vanish when
`alpha,beta` are nonzero real coefficients: `U,V` are independent
over `R`. The possible obstruction is a *very short nonzero*
combination, not an exact cancellation. Its squared norm is the
positive integral binary quadratic form

```text
Norm(alpha U+beta V)
 =Norm(U)alpha^2+2Re(U bar(V))alpha beta+Norm(V)beta^2,
discriminant=-4 delta_abc^2.                         (10)
```

The coefficients of this form can be of order `exp(8w)`, while its
discriminant is at most `4exp(24sigma)T^6`. Completing the square
gives, for an integer shear `t`,
`Norm(V+tU)=Norm(U)(t+Re(U bar(V))/Norm(U))^2
+delta_abc^2/Norm(U)` exactly. This is the same short-vector
geometry as the fixed-column shear in
[contact_lattice_shortest_lift.md](contact_lattice_shortest_lift.md),
Section 3. The small discriminant supplies no deterministic lower
bound on the least integer shear representative; the existing
contact-lattice theorem likewise isolates, rather than excludes,
an exceptionally short CRT representative. This identity does not
claim that `U,V` are themselves that contact frame.

This residual core does not repeat at another inside pair. In
`A_abc`, the oriented singleton block `H_{ {b} }` and the opposite
pair block `bar(H_{ {a,c} })` occur; in `B_abc` the corresponding
private signature is `H_{ {a} }` and `bar(H_{ {b,c} })`. Either
signature identifies the entire outside triple `{a,b,c}`, hence its
inside complement `S`. All twenty residual monomials from the ten
inside pairs are distinct. A literal stacking of different cut
congruences therefore cannot act on the same residual monomial;
any common scalar must come from a further, genuinely arithmetic
relation between their two-term cancellations.

The inside cut `H_S` has membership `000`, so it divides none of
`C_abc,A_abc,B_abc`; the one-pass congruence genuinely constrains
the bracket in (7), up to the small correction factor. The
complementary cut `H_(S^c)` has membership `111`, so it already
divides `C_abc` as `H_(S^c)^2 bar(H_(S^c))`. Multiplying a one-pass
`H_S` congruence by this *trivial* original numerator divisibility
does not reproduce the reflected arithmetic modulus.

Even if one could transfer both nontrivial tests to the same bracket
in (7) at a cost `O(sigma+log T)`, the paired cut product would have
Gaussian modulus about `exp(w)`, while the bracket has the allowed
upper modulus `exp(4w+O(sigma+log T+log C_Q))`. A direct
archimedean size-versus-divisibility argument would still miss by
`3w`. Combining three or five boundary restrictions needs an
additional identity tying their different brackets or their unit
residue classes; (4) and (7) alone do not provide one.

## 3. What the complement reflection does to the same numerator

Use the exact complement reflection of
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md):
`P_i'=Z bar(P_i)/(Norm(P_i)d_i)`, where `Z` is one common Gaussian
integer and `d_i` is a positive ordinary content. Put
`rho_i=1/(Norm(P_i)d_i)`. Then the reflected outside determinant is

```text
D_ij'=-Norm(Z) rho_i rho_j D_ij,
V_S'=alpha D_ac' P_b'+beta D_bc' P_a'
    =-Norm(Z) Z rho_a rho_b rho_c bar(V_S).            (11)
```

This is a shared-chart identity, not an independent assignment of
the reflected residues. After primitive normalization,
`t_ij'=-t_ij` exactly. The large factor and the three rational
denominators in (11) remove the original trivial
`H_(S^c)^2 bar(H_(S^c))` factor and expose the reflected
restriction test. They also carry cut-scale height, so (11) cannot
turn the two generic congruences into divisibility of one small
ordinary coefficient. In the nodal single-difference case,
the different scalar argument in the boundary-node note does yield
that coefficient modulus; condition (2) is exactly the surviving
case where it does not.

## Verification and limit

Run [check_five_row_shared_triplet_residue_obstruction.py](check_five_row_shared_triplet_residue_obstruction.py).
It checks (4), the full 32-cut exponent table and eight-block
residuals in (6)--(7), exact integer primitive residues, (9)--(10),
and (11)
on literal conjugate-primitive Gaussian rows. The fixture does not
have controlled small `T` or endpoint arc angles; those are retained
as hypotheses in the algebraic estimates, not inferred by the
finite checker. No cross-boundary coefficient-height gap is claimed.
