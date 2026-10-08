# A rational lift forces the full branch translation kernel

Let `C/Q` be the genus-five curve `y^2=f(x)`, with twelve distinct
branch roots, and suppose their Galois group is S_12. Let K be the
branch splitting field. If a normalized maximal unramified two-cover
`X -> C` has a rational point and

```text
disc(f) is not -1 in Q*/Q*2,                            (1)
```

then the affine Galois action on every branch fiber has full
translation kernel `J[2]`. All branch points on X form one orbit
of degree `12*1024=12288`. By the injectivity theorem in
[the osculating hyperplane note](six_point_descent_osculating_hyperplane.md),
their order-twenty hyperplanes have that same field degree.

Thus the degree-12 and degree-132 osculating cases cannot occur on
these **pointed** normalized twists. This is a restriction on the
applicability of the low-degree height criterion, not an endpoint
height bound. It makes no claim for arbitrary twists without a
rational point, or when the discriminant squareclass is -1.

## 1. The rational point identifies the based cover

Let `x in X(Q)` lie above `P in C(Q)`. The multiplication-by-two
pullback based at P is

```text
X_P = {(R,Q) in C x J : 2Q=[R-P]}.
```

It has the rational point `(P,0)`. Both covers are geometrically the
connected maximal elementary abelian unramified two-cover of C.
Any geometric C-isomorphism between them can be translated by a deck
transformation to send x to `(P,0)`. There is exactly one such pointed
isomorphism: the deck group acts simply transitively on the fiber.
Its Galois conjugates have the same pointed property, so it descends
to Q. Consequently

```text
(X,x) is Q-isomorphic to (X_P,(P,0)).                    (2)
```

This uses the actual rational point, not a hypothetical rational
identification with a previously fixed based cover. The cover and
base-point constructions are recorded in
[the base-point audit](six_point_two_cover_basepoint_audit.md).

## 2. Zero branch kernel would rationalize all four-torsion

Suppose the translation kernel of the fiber above one branch point
`W_j=(alpha_j,0)` were zero. Over `Q(alpha_j)`, the linear action is
S_11 and becomes trivial over K. Zero translation kernel therefore
means that every point of this fiber is K-rational. By (2), choose

```text
Q_j in J(K),           2Q_j=[W_j-P].
```

Because K is normal over Q, the roots form one Galois orbit, and P
is rational, conjugation gives such a Q_j for every branch label.
For any pair j,k,

```text
Q_j-Q_k in J[4](K),     2(Q_j-Q_k)=[W_j-W_k].           (3)
```

The branch differences generate `J[2]`, and `J[2]` is already
K-rational. Thus multiplication by two maps `J[4](K)` onto all of
`J[2]`, with kernel all of `J[2]`. Its cardinality is consequently
`|J[2]|^2=|J[4]|`, proving

```text
J[4] subset J(K).                                       (4)
```

The canonical principal polarization of J and the nondegenerate,
Galois-equivariant Weil pairing then force `mu_4 subset K`; see
[Milne, Abelian Varieties, I Section 13](https://www.jmilne.org/math/CourseNotes/AV.pdf#page=63).
In particular, `i in K`.

An S_12 extension has exactly one quadratic subfield, its sign field
`Q(sqrt(disc(f)))`. Hence `i in K` would imply that disc(f) has
squareclass -1, contradicting (1).

The S_11 module `J[2]` is simple, so its branch-fiber translation
kernel is either zero or all of J[2]. The zero case is excluded,
leaving the full kernel. Translations are transitive on the 1024
points in a branch fiber; root transitivity then gives the asserted
absolute orbit size 12288. The simplicity argument and the complete
unpointed orbit alternatives appear in Section 4 of the osculating
note.

## 3. The existing first S_12 witness satisfies (1)

For the first-ruling witness with central weights

```text
(1,-20,250,-1000,1445,-676),
```

the full S_12 certificate is in
[the constant-cover audit, Section 1](six_point_constant_two_cover_archimedean_obstruction.md).
Use the primitive integral polynomial reconstructed by
`check_six_point_branch_two_torsion.py`. At the good prime 59,

```text
disc(f) = 22 mod 59,
22^29 = 1 mod 59,          (-22)^29 = -1 mod 59.          (5)
```

Therefore `-disc(f)` is not a rational square, proving (1). There
is also a group-theoretic check: the already certified factorization
type `(1,11)` at 59 is even, whereas -1 is a nonsquare modulo 59.
The exact modular discriminant calculation (5) has been added to
the existing checker.

For this witness, every normalized twist carrying an actual rational
lift has branch and osculating-hyperplane degree 12288. This does not
assert the full-S_12 or discriminant hypotheses for arbitrary fair
weights, and it gives no upper bound on an isolated endpoint point.

## 4. The complete contact budget for homogeneous auxiliary forms

The full branch degree also rules out a power saving from the entire
following class of single-orbit contact arguments, not merely the
order-twenty hyperplane. Let F be a homogeneous form of degree d in
the rational beta coordinates, with coefficients in a number field
E of degree e. Assume `F|_X` is not identically zero. After clearing
coefficients, form the rational norm section

```text
S = product_(sigma:E->Qbar) F^sigma |_X,
S in H^0(X,O_X(ed)),                r=ord_B(S).
```

The geometric integrality of X makes S nonzero. Importantly, r is
the **sum** of the contacts of all conjugate forms at B, not just
the contact of the initially chosen F. Thus it already incorporates
every small conjugate value that an algebraic norm could exploit.
Since S is rational, it has the same order r at all D conjugates of
B. A nonzero section of `O_X(ed)` has a zero divisor of degree
`1024ed`, giving the exact budget

```text
D r <= 1024 e d.                                       (6)
```

For fixed auxiliary degree with polynomially controlled coefficients,
the endpoint bound `|zeta|<=U^C T^-10` and local contact give

```text
|S(z)| <= C(C_end) U^C T^(ed-10r),
ed-10r >= ed(1-10240/D).                                (7)
```

When `D=12288`, the last exponent is at least `ed/6>0`.
Consequently increasing the auxiliary degree or choosing a smaller
coefficient field cannot make this contact-only norm comparison
force an upper bound on T. If S(z)=0, its value cannot be used as a
nonzero integral norm in the first place. The calculation needs no
assumption that E contains the residue field of B.

This obstruction is restricted to homogeneous auxiliary sections and
the lower bound from a nonzero integral norm, using only the stated
endpoint window and vanishing at this branch orbit. It does not
exclude extra finite-place divisibility, a stronger actual local
approximation, rational functions with additional pole information,
or constructions comparing multiple target points.
