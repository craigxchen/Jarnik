# Smaller-cut congruences at every depth

The first-smaller-cut calculation extends to a finite hierarchy. At
depth `k`, a relation whose preceding restrictions vanish specializes
to a squarefree translation-invariant polynomial of degree `2k`.
Its maximal coefficient minors have the same leading complementary
modulus `exp(2w)` as at depth one. The correction loss is independent
of depth; the residue loss grows linearly with depth.

This is a conditional restriction theorem at each level. A proper
restriction image is not a zero image, so the theorem does not put
short relations successively into the next kernels. That implication,
needed to turn the finite hierarchy into a contradiction, is unproved.

The elementary polynomial certificates are in
[check_all_cut_invariant_relation_hierarchy.py](check_all_cut_invariant_relation_hierarchy.py).
The arithmetic uses the already audited local gcd bounds and
[complement reflection](complement_reflection_invariant_relations.md).

## 1. The finite restriction hierarchy

Let `m=2q>=6`. Let `V` be the rational simultaneous SL2 invariant
space of degree two in each labelled binary row. For an inside set
`S`, specialize its rows to `(0,1)` and the other rows to `(1,z_j)`;
write the resulting polynomial as `F_(S,Q)`. This choice of the two
coordinate vectors is equivalent to the swapped convention in earlier
notes by one common SL2 transformation.

Define `K_k` to be the subspace on which all restrictions with
`|S|>=q-k+1` vanish. Thus `K_0=V` by diagonal weight balance, and
`K_1` is the balanced kernel. For `Q in K_k`, take `|S|=q-k`, and put

```text
n=q+k,       d=2k.
```

Then `F_(S,Q)` is homogeneous of degree `d`. Indeed an invariant has
total second-coordinate degree `m`, and the collapsed rows account
for `2|S|` of it. The coefficient of `z_j^2` is the restriction at
`S union {j}`, so it is zero. Finally the common unipotent map
`(x,y)->(x,y+tx)` fixes `(0,1)` and translates every outside `z_j`.
Therefore `F_(S,Q)` belongs to the space

```text
W_(n,d)={squarefree homogeneous f of degree d: sum_j partial_j f=0}.
```

The map consisting of all these restrictions on `K_k` has kernel
exactly `K_(k+1)`.

For completeness, let `D` delete a variable from squarefree monomials,
and let `U` add an absent variable. Use the inner product making those
monomials orthonormal. They are adjoint, and on degree `d`,

```text
D U-U D=(n-2d) Id.
```

If `d>n/2`, this identity implies
`||D f||^2=||U f||^2+(2d-n)||f||^2`, so `D` is injective.
If `d<=n/2`, `U` is injective on degree `d-1`, so its adjoint `D`
onto degree `d-1` is surjective. Consequently

```text
dim W_(n,d)=binomial(n,d)-binomial(n,d-1)  if d<=n/2,
W_(n,d)=0                              if d>n/2.       (1)
```

In particular a nonzero level requires `3k<=q`. For larger `k`,
`K_k=K_(k+1)`. Propagating to the empty inside set, whose restriction
determines the whole multihomogeneous polynomial, gives the exact
termination statement

```text
K_(floor(m/6)+1)=0.                                    (2)
```

This proves termination of the algebraic kernels, not membership of
small numerical relations in those kernels.

The later [Specht filtration theorem](kernel_specht_filtration.md)
identifies the exact permutation constituents and dimensions of every
`K_k`. Independently, the
[higher-kernel content theorem](higher_kernel_evaluation_heights.md)
computes their conductor orders and primitive evaluation heights by
monomial support. Neither theorem supplies the missing membership of
short relations in the next kernel.

The [full-restriction theorem](full_invariant_restriction_image_and_index.md)
shows that each active positive-depth restriction is integrally onto its
matching lattice. At distinct directions the full numerical-zero space
still has full rational restriction image; its integral image has cyclic
quotient of order equal to an exact evaluation-gcd ratio. These facts
separate full-space rank from the unresolved short-relation coverage.

The index is sharp. Write `m=6ell+2r`, with `r in {0,1,2}`, and
take a product of `2ell` disjoint triangle invariants and `r`
disjoint squared brackets. A nonzero restriction can collapse at
most one row in each component, hence at most `2ell+r=q-ell`
rows. Collapsing exactly one from each component gives a product
of outside differences, so this polynomial is in `K_ell` and is
not in `K_(ell+1)`. It is nonzero at all distinct actual directions;
this is an algebraic witness, not a numerical zero relation.

## 2. An integral matching basis

For `d<=n/2`, index subsets `B={b_1<...<b_d}` of `{1,...,n}` by
the condition `b_j>=2j`. For each such subset choose distinct
`a_j<b_j` from its complement; greedy selection is possible since
there are at least `j` complement elements below `b_j`. Set

```text
w_B(z)=product_(j=1..d)(z_(b_j)-z_(a_j)).                (3)
```

These are squarefree and translation invariant. There are
`t=binomial(n,d)-binomial(n,d-1)` such subsets by the elementary
ballot reflection count. Order squarefree monomials by the binary
weight of their index set. Every monomial of `w_B`, except `z_B`,
has smaller weight: replacing any `b_j` by its smaller `a_j`
strictly decreases the sum of powers of two. Thus the coefficient
matrix on the indicated subsets is triangular with diagonal one.
By (1), these polynomials form a basis.

The triangular matrix is integral unimodular, so any integer
polynomial in `W_(n,d)` has integer coordinates in this basis.
Its off-diagonal entries have absolute value at most one. Successive
triangular elimination gives the convenient fixed bound

```text
F_(S,Q)=sum_B c_(B,Q) w_B,       c_(B,Q) in Z,
sum_B |c_(B,Q)| <= L_(m,k) C_Q,
L_(m,k)=2^(t-1),                                        (4)
```

where `C_Q` is the coefficient sum of the original invariant. The
restriction itself does not increase this sum. The constant in (4)
is deliberately loose and depends only on the fixed row count and
depth. At depth one the earlier basis has the better constant one.

Every basis element is a single matching on `2d=4k` distinct outside
rows. Thus its numerical value is nonzero when the actual directions
are distinct. This avoids assuming that a sum of matching products
is a unit or even nonzero.

## 3. Actual integer moduli for every maximal minor

Use the actual Gaussian data of the preceding notes:

```text
P_i=K_i product_(T containing i)H_T,
gcd_G(P_i,bar(P_i))=1,
log|K_i|<=sigma,       0<|t_ij|<=T,
log N(H_T)>=(1-eta)w.
```

All core blocks and their conjugates have disjoint odd split-prime
support. Let `Q in K_k` be an integer numerical zero and fix an
inside set of size `q-k`. Put `H=H_S`, and on the outside rows put

```text
z_j=Y_j/P_j,
v_B=(product_outside P_j) w_B(z),
kappa=product_outside K_j.
```

The `v_B` are Gaussian integers. In fact each equals, up to sign,
the product of `2k` distinct-row brackets `D_ab`, multiplied by
the `n-4k` unused outside rows `P_j`. It is nonzero.

The common determinant-one change `(X,Y)->(X+iY,Y)` means the
numerical zero can be evaluated on `(P_i,Y_i)`. Reducing its inside
first coordinates modulo `H`, and removing the inside `Y_i` and
outside core units exactly as in the first-smaller-cut proof, gives

```text
H | kappa sum_B c_(B,Q) v_B.                            (5)
```

For any `t` numerical zero relations in `K_k`, let `delta_S` be
the determinant of their coefficient matrix in (4). An integer
adjugate applied to (5) gives `H | kappa delta_S v_B` for each `B`.
Fix one basis element `B_0`. Since the determinant is an ordinary
integer, the exact norm-divisibility rule gives

```text
M_S | delta_S,
M_S=N(H)/N(gcd_G(H,kappa v_(B_0))).                     (6)
```

This same modulus divides every maximal minor of any rectangular
coefficient matrix of numerical zero relations at this cut.

The logarithmic losses in (6) are at most `2n sigma` from `kappa`,
`2(n-4k)sigma` from the unused outside rows, and
`4k sigma+2k log T` from the `2k` outside brackets. Core factors in
the unused outside rows are units at `H`. For brackets use the exact
bound `gcd(N(H),D_ab)<=exp(2sigma)T`, retaining prime powers.
Since `n=q+k`, their sum simplifies to

```text
log M_S >= log N(H_S)-2m sigma-2k log T.                 (7)
```

Apply complement reflection to the same polynomials and labelled
restriction matrix. The reflected rows have correction bound
`(2m-1)sigma`, preserve the pair residues up to sign, and retain the
complementary block with norm loss at most `2m sigma`. The second
modulus divides `N(H_(S^c))`, so it is coprime to the first. Thus

```text
M_S M_S' | delta_S,
log(M_S M_S') >= log N(H_S)+log N(H_(S^c))
                -2m(2m+1)sigma-4k log T.               (8)
```

The correction-height coefficient is independent of the depth.
Reflection does not change membership in the algebraically defined
`K_k`, or any coefficient of any invariant relation.

By (4) and Hadamard,

```text
|delta_S|<=L_(m,k)^t product_(j=1..t) C_(Q_j).           (9)
```

It follows that relations in `K_k` with coefficient heights at most
`h` have proper restriction images at every cut of size `q-k` if

```text
t(h+log L_(m,k))+2m(2m+1)sigma+4k log T < 2(1-eta)w.   (10)
```

No individual coordinate or bracket has been inverted modulo a
prime power. The only removed powers are charged in (7)--(8).

## 4. What would still be required for uniformity

Equation (2) makes an induction attractive, but (10) does not provide
its inductive step. It bounds the rank of the restriction image by
`t-1`; it does not make that image zero. Even numerical-zero spaces
obtained by multiplying smaller-support zero relations by nonzero
triangle invariants can have proper images at many cuts.

For `m=6`, the special exterior-product injectivity says that proper
restriction images of a two-dimensional relation space at every
first-smaller cut would force global proportionality, as proved in
the earlier six-row note. No analogue closing the hierarchy for arbitrary `m`
is proved here. The result extends the available actual-integer
constraints while leaving the uniform point-count objective open.

The later [eight-row cofactor congruence](eight_row_segre_cofactor_congruence.md)
uses the same conductor moduli to force a cubic identity in the raw
maximal minors of four short restriction rows. This is a nonlinear
necessary condition beyond proper rank; it allows boundary points and
does not close the induction.

The [Chow basepoint criterion](first_smaller_cut_chow_basepoint.md)
extends this nonlinear constraint to general first-smaller cuts. For
fixed `m=400` it applies to a large space of relations supplied by
successive minima. Common boundary points and zero images are still
allowed, so the kernel-membership gap remains.

Independent audit: `symmetry_descent` checked the kernel termination,
integral basis, actual specialization, and both prime-power loss
calculations. The root agent checked the exact polynomial certificates.
