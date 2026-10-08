# Which symmetries can amplify the same restriction modulus?

Complement reflection supplies two coprime moduli for one labelled
coefficient determinant. Reanchoring does not supply a third large
independent block for that determinant. The statement can be made exact,
with arbitrary prime powers and without bounding the common Gaussian
rotation used in reanchoring.

The new support estimate below applies to common rotations, reflections,
and rational row scalings of the same labelled directions. General
nonisometric projective maps are discussed separately; neither their
exclusion nor endpoint uniformity is asserted here.

## 1. Signed valuations remove all ordinary row content

Let the fixed labelled rows have independent core factorizations

```text
P_i=K_i product_(S containing i) H_S,     1<=i<=m,
log|K_i|<=sigma.
```

The blocks have mutually disjoint odd split prime support, including
conjugate support. At a rational prime `p`, choose an orientation `pi`
used by the old core, if any. Write its old block as `pi^e` with label
`S`; if there is no old core factor, put `e=0`. Set

```text
n_i=v_pi(P_i)-v_barpi(P_i)=e 1_(i in S)+delta_i,
E_p=sum_i (v_pi(K_i)+v_barpi(K_i)).
```

Then `sum_i |delta_i|<=E_p`, and

```text
sum_p E_p log p <= sum_i log N(K_i) <= 2m sigma.       (1)
```

Suppose the same labelled directions are transformed by

```text
P_i'=r_i lambda P_i       or       P_i'=r_i lambda bar(P_i),
lambda in Q(i)^*,       r_i in Q^*.
```

The resulting rows can be made integral and primitive; these choices
do not affect signed valuations. At each prime the exact transformation
is

```text
n_i'=c_p+epsilon n_i,       epsilon in {+1,-1},          (2)
```

because an ordinary rational row multiplier has equal valuations at
`pi` and `barpi`. No estimate on `lambda` or on the row contents is
needed in (2).

Assume the new rows also have an independent core factorization with
corrections `K_i'`, `log|K_i'|<=sigma'`. A new core power over `p`
has exponent `f`, orientation sign `epsilon'`, and dividing-row set
`T`. Thus

```text
n_i'=epsilon' f 1_(i in T)+delta_i',
sum_i |delta_i'|<=E_p',
sum_p E_p' log p<=2m sigma'.                            (3)
```

For every nonempty proper `T`, equations (2)--(3) imply

```text
f <= e 1_(S=T or S=T^c) + E_p+E_p'.                    (4)
```

Here constant old patterns and `e=0` count as incompatible with `T`.
To prove (4) for an incompatible pattern, some two rows have the same
old bit and different new bits. Subtracting their equations cancels
both `c_p` and the old core contribution, and bounds `f` by the four
correction absolute values. For a compatible pattern, subtract rows
on opposite sides: the old core difference has absolute value `e`,
giving `f<=e+E_p+E_p'`. This proves the claim even when the transformed
rows are not conjugate-primitive before their ordinary contents are
removed.

Summing yields the useful fixed-label bound

```text
log N(H_T') <= log N(H_T)+log N(H_(T^c))
                 +2m(sigma+sigma').                   (5)
```

Thus a common isometry cannot move a third block's full weight into
the restriction at `T` while keeping correction height negligible.
This is a statement about actual valuations, not merely Boolean
patterns or a generic transformation.

For several such transformations indexed by `a`, let `M_a` be any
ordinary modulus dividing the corresponding `N(H_T^(a))`. Taking
the maximum prime valuation in (4), rather than summing old block
weights repeatedly, gives

```text
log lcm_a M_a <= log N(H_T)+log N(H_(T^c))
                    +2m sigma+2m sum_a sigma_a.        (6)
```

In a full profile the leading term is `2w`. If the total correction
budget on the right is `o(w)`, any number of these applications still
has leading modulus at most `2w`, not a growing multiple of `w`.
For many transformations their correction costs must actually be
summed; (6) does not discard them.

## 2. What reanchoring does on the same rows

For primitive old anchor numerators, reanchoring at an active row `a`
while retaining the same labelled `m` directions gives

```text
P_i' = P_i bar(P_a)/N(gcd_G(P_i,P_a)).                  (7)
```

This is exactly a common rotation with positive rational row scales.
The pivot row becomes `P_a'=1`. At a core prime, before corrections,
the signed load is `e(1_(i in S)-1_(a in S))`.
Consequently its support is `S` if `a` is outside `S`, and its
conjugate support is `S^c` if `a` is inside. Complementary old
blocks merge into one new dividing-row set excluding `a`.

Equation (5) is the correction-sensitive version of this statement.
It includes the exact common divisor in (7), and so incurs no hidden
denominator or individual row rotation.

The modulus theorem in
[complement_reflection_invariant_relations.md](complement_reflection_invariant_relations.md)
already obtains the full two-block leading weight. Reanchoring is
an alternative access to this support, not an independent third
source of prime powers.

## 3. The larger Boolean anchor group does not act on a fixed relation

For an underlying configuration labelled `0,...,m`, represent an
unoriented cut by binary coordinates with `b_0=0`. A permutation
`sigma` of all `m+1` labels followed by choosing its row zero as
anchor acts by

```text
b_j' = b_(sigma(j)) + b_(sigma(0))  in F_2,
1<=j<=m.                                             (8)
```

These are the permutation action on
`F_2^(m+1)/<all-ones>`, a faithful `S_(m+1)` action for `m>=2`.
For the transposition swapping zero and `a`, the explicit formula is

```text
b_a'=b_a,       b_j'=b_j+b_a  for j!=a.                (9)
```

This can change the normalized Hamming weight, but it preserves
the unordered sizes of the two sides in the full `m+1` row set.
For example, with six active rows the orbit of a size-two cut has
21 members, of normalized sizes two or five; it is not all 64
normalized patterns.

There is a tempting formal amplification: adjoining the translation
`b -> b+all-ones` to this group generates the affine group
`F_2^m semidirect S_(m+1)`. Indeed (9) sends the all-ones vector to
the `a`th coordinate vector, which generates all translations.
This is an exact statement about the formal Boolean operations.
It is not an action preserving an arbitrary fixed numerical
relation on the `m` active directions.

The obstruction is the active row set. To reanchor the complete
`m+1`-point configuration and still write exactly `m` nonanchor
rows, one removes the former active pivot and inserts the former
anchor. An invariant zero in the original `m` directions need not
remain zero after this replacement. Common covariance applies to
the same directions, not to inserting a different point.

An exact example uses `P_i=2i+i_G` for `i=1,2,3,4`, where `i_G`
is the Gaussian imaginary unit, and external anchor `P_0=1`. All
four rows have primitive coordinates and odd norm. The invariant

```text
Q=(4 Delta_12 Delta_34-Delta_13 Delta_24)^2
```

has degree two in every row and vanishes on these four rows.
Replacing the first row by `P_0` gives `Q=16`, not zero.
Applying the common reanchoring rotation afterwards multiplies
this by a nonzero covariance scalar and cannot restore the zero.

The complement reflection on the active rows also moves the external
anchor: its image is the common multiplier rather than the newly
chosen real anchor. Treating this external point as if it were
unchanged is precisely what makes the formal affine action appear
to preserve more data than it does.

## 4. Relabeling coefficients moves the determinant with its modulus

One may of course relabel both the rows and the polynomial. If old
row `i` is renamed `sigma(i)`, define the new polynomial by the
opposite variable substitution. The numerical zero and coefficient
sum are preserved. The restriction at the renamed set `sigma(S)`
then has the same coefficient matrix as the old restriction at `S`,
up to the integral basis change induced by relabeling outside rows.
That change is unimodular, so maximal determinants agree up to sign.

But its block is also the renamed old block `H_S`. If instead the
target label `S` is held fixed, its determinant is the old one at
`sigma^(-1)(S)`, not the original determinant at `S`.
Thus this operation either preserves both determinant and prime
support, or moves both of them together.

Multiplying such different determinant inequalities does not by
itself amplify the threshold. For `r` nonzero determinants formed
from `t` relations of logarithmic coefficient height at most `h`,
their product has upper height `r t h` while the corresponding
two-block lower bounds sum to approximately `2r w`. The factor
`r` cancels. A special relation space with additional permutation
identities could change this conclusion, but no such identity is
implied by a numerical zero or by invariant covariance.

## 5. Nonisometric rational projective transformations

Every common rational `GL_2` transformation preserves the same
invariant zeros. Unlike a common rotation or reflection, it need
not preserve the primitive pair residues. It must therefore be
audited using the new residues and new corrections, not inserted
freely into (6).

There is already an exact useful low-height restriction. After
reanchoring at the image of the old real anchor, write an integral
matrix `A` of entry height `H` as

```text
q=a^2+c^2,    r=ab+cd,    s=ad-bc,
u=q+s-ir,    v=q-s+ir,
V_i=q Re(P_i)+(r+is) Im(P_i).
```

The audited Gaussian divisibility in
[projective_radius_inflation.md](projective_radius_inflation.md)
is

```text
gcd_G(V_i,V_j) | q u v t_ij,
|gcd_G(V_i,V_j)| <= 40 H^6 T        if uv!=0.           (10)
```

Primitive output representatives divide these raw rows. Hence any
new independent block incident to at least two active rows satisfies

```text
log N(H_S') <= 12 log H+2log T+2log 40.                (11)
```

For first-smaller cuts `|S|=m/2-1>=2`, this excludes a full block
of weight `w` in this reanchored chart when `log H,log T=o(w)`.
If `uv=0`, the map is a rotation or reflection and (5)--(6) apply.

High-height nonisometric maps, or a subsequent large rotation of
their outputs changing the arithmetic chart, are not excluded by
(11). No theorem here asserts that such maps cannot produce extra
useful moduli. The precise remaining opportunity is to construct a
new chart with new primitive residues and corrections still small,
and demonstrate extra support for the same determinant after all
of those costs are retained. The Boolean anchor operations alone
do not provide it.

The prime-power support inequality, formal Boolean action, and
explicit relation-replacement failure are checked in
[check_reanchoring_relation_moduli.py](check_reanchoring_relation_moduli.py).

Independent audit: `fresh_algebraic` checked the signed-valuation and
lcm bounds, anchor group, exact replacement example, and unimodular
relabeling statement; no gap was found.
