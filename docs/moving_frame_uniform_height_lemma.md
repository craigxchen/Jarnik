# A uniform height lemma for rational maps of bounded presentation degree

This note proves the height estimate needed for rational frames whose functions may move with the index. The presentation degree is bounded, and the arithmetic height of the coefficients is explicitly charged to the error. No fixed-map height-comparison constant is silently reused for a moving map.

The lemma is a height statement on one fixed elliptic curve. It does not produce a frame from a lattice arc or bound the complexity of any frame obtained by other means.

## 1. Precise statement

Fix an elliptic curve `E/Q`, a rational origin `O`, and a nonsingular projective Weierstrass model. Let `q` denote the canonical height associated to `[O]`, normalized by

```text
h_x(Q) = 2q(Q)+O_E(1).
```

All heights below are absolute logarithmic Weil heights, with the usual normalized local weights over any number field of definition.

Fix a positive integer `D`. Let `f_1,...,f_4` be rational functions on `E` with presentations

```text
f_i = F_i/G_i,
```

where `F_i,G_i` are homogeneous rational polynomials of the same degree at most `D` in the fixed projective coordinates, and `G_i` does not vanish identically on `E`. Affine numerator/denominator presentations of bounded total degree are included by homogenizing numerator and denominator to a common degree. Zero numerators are allowed.

Define the presentation coefficient height by

```text
H = max_i h([coefficients of F_i : coefficients of G_i]).
```

The numerator and denominator coefficients are taken together in each projective vector. Thus their relative scale is included; independent rescaling of the numerator and denominator is not ignored. Other ordinary coefficient-height conventions are equivalent up to constants depending on `D`.

Let

```text
Delta = sum_B max(0,-ord_B(f_1),...,-ord_B(f_4)) [B],
d = deg Delta,
Phi = [1:f_1:f_2:f_3:f_4] : E -> P^4.
```

Here the divisor may be written over an algebraic closure, with multiplicities. The rational map extends to a morphism. Its pullback of `O(1)` is `O_E(Delta)`, because the five sections `1,f_1,...,f_4` have no common zero after clearing the least common pole divisor. In particular `d <= 12D`.

**Uniform height lemma.** There is a constant `C_(E,D)`, depending only on the fixed curve/model and `D`, such that for every such presentation and every `Q in E(Qbar)`,

```text
|h(Phi(Q)) - d q(Q)|
    <= C_(E,D) (sqrt(q(Q)(H+1)) + H+1).                 (1)
```

At poles, `Phi(Q)` means its morphism value. In particular,

```text
h(Phi(Q)) >= d q(Q)
              - C_(E,D) (sqrt(q(Q)(H+1)) + H+1).       (2)
```

The determinant-one identity is not needed for this lemma. It is used separately when contacts force a lower bound on `d`.

## 2. Uniform heights of the pole points

Each pole of `f_i` is among the zeros on `E` of `G_i`. This is a proper intersection because `G_i` does not vanish identically on `E`. The intersection degree is at most `3D`.

Bounded-degree elimination, or the arithmetic Bezout estimate for a fixed projective curve intersected with a hypersurface, gives

```text
h_E(B) <= C_(E,D)(H+1)                                (3)
```

for every point in these four intersections. This applies to every conjugate and therefore to every point in the support of `Delta`. One can see the coefficient dependence directly: eliminate against the fixed cubic, use a nonzero bounded-degree resultant in a coordinate chart, and apply the elementary height bound for roots of a polynomial. All resultant coefficients have height `O_(E,D)(H+1)`. A finite collection of charts/eliminants handles vanishing leading coefficients and points at infinity.

Write the geometric divisor, repeating points according to multiplicity, as

```text
Delta = [B_1]+...+[B_d].
```

Since the projective Weierstrass embedding has divisor class `3[O]`, (3) also implies

```text
q(B_j) = O_(E,D)(H+1).                                 (4)
```

The degrees of all pole points, and the number of their conjugates together, are bounded in terms of `D`. Consequently they all lie in one number field of degree bounded in terms of `D`; for example a splitting field obtained from the at most `12D` geometric denominator zeros has degree bounded by `(12D)!`. No height constant is allowed to depend on an uncontrolled extension degree.

## 3. Centering the divisor in the fixed linear system `|d[O]|`

Assume first that `d > 0`. Let

```text
S = B_1+...+B_d
```

be the elliptic-curve group sum. Choose `T` with `[d]T=S`. Such a point exists over an algebraic extension of degree at most `d^2` of the field already chosen. Positivity and Cauchy–Schwarz for the canonical height give

```text
q(S) <= d sum_j q(B_j),
q(T) = q(S)/d^2 = O_(E,D)(H+1).                        (5)
```

Enlarge the bounded-degree field to contain `T`. Translation formulas on the fixed curve show that each translated function

```text
f_i'(X) = f_i(X+T)
```

has a presentation of bounded degree and coefficient height `O_(E,D)(H+1)`. Indeed, the addition law has fixed degree, and the ordinary coordinate height of `T` is `O_(E,D)(H+1)` by (5). The identity translation and exceptional affine charts are handled using finitely many fixed formulas.

The minimal common pole divisor of these translated functions is

```text
Delta' = sum_j [B_j-T].
```

Its group sum is `S-[d]T=O`. The standard identification `Pic^0(E)=E` therefore gives

```text
Delta' ~ d[O].                                        (6)
```

All support points of `Delta'` still have bounded degree and height `O_(E,D)(H+1)`.

## 4. A height-controlled principalizing function

For each integer `d` in the fixed range `1 <= d <= 12D`, fix a rational basis of `L(d[O])`. By (6) there is a nonzero `g in L(d[O])` with

```text
div(g) = Delta' - d[O].                               (7)
```

The coefficients of `g` in this fixed basis can be chosen with projective height `O_(E,D)(H+1)`. Here is a direct justification, including multiplicities.

Require the section of `O_E(d[O])` represented by `g` to vanish to the specified multiplicity at each support point of `Delta'`. These are homogeneous linear conditions on its basis coefficients. They are jet-evaluation conditions of order at most `d`, hence of bounded order. Use a fixed finite set of nonsingular coordinate charts and local parameters; at `O`, use a fixed parameter and the section trivialization of `O_E(d[O])`. The matrix entries are rational expressions of bounded degree in the coordinates of the support points. Their heights are therefore `O_(E,D)(H+1)`.

The kernel is one-dimensional: it is `H^0(E,O_E(d[O]-Delta'))`, which is the space of sections of a trivial line bundle. A nonzero kernel vector obtained from a nonzero maximal minor consequently has coefficient height `O_(E,D)(H+1)`. This argument also handles a contact at `O`; there the condition is vanishing of the section, not incorrectly imposed vanishing of a rational function already permitted to have a pole.

Now put

```text
s_0 = g,
s_i = g f_i'  (1 <= i <= 4).
```

These are sections of `O_E(d[O])` and have no common zero. This follows equally from (7) and from minimality of the common pole divisor: outside `Delta'`, the original section `1` is nonzero; at every point of `Delta'`, an entry achieving its maximal pole order remains nonzero after clearing that pole.

Expressed in the fixed basis of `L(d[O])`, all five sections have coefficient height `O_(E,D)(H+1)`. To make this quantitative, multiply the bounded-degree presentations for `g` and `f_i'`, clear denominators in the identity with an unknown basis expansion, reduce modulo the fixed Weierstrass equation, and solve the resulting bounded-size linear system by a nonzero minor. Products, elimination, and these minors preserve the asserted linear dependence on `H+1`.

Thus the translated morphism

```text
Psi(X) = Phi(X+T) = [s_0(X):...:s_4(X)]                (8)
```

is represented by a basepoint-free tuple of controlled height in a fixed complete linear system.

## 5. Uniform comparison for a basepoint-free tuple

For `d >= 3`, the fixed basis gives an embedding

```text
iota_d : E -> P^(d-1).
```

The sections in (8) are linear combinations of that basis, with coefficient matrix `M` of height `O_(E,D)(H+1)`. Consequently `Psi=M iota_d` on the fixed embedded curve. Although `M` can approach the basepoint locus as its coefficients vary, the arithmetic height of this approach is controlled by its coefficients.

More explicitly, let `I_d` be a fixed homogeneous ideal defining the embedded curve. Basepoint freeness says that `I_d` together with the linear forms `(MX)_j` has no projective zero. A bounded-degree homogeneous Nullstellensatz gives an integer `K`, bounded solely in terms of the fixed embedded curve, and identities

```text
X_i^K = sum_j A_ij(X) (MX)_j   modulo I_d,
deg A_ij = K-1.                                       (9)
```

The coefficient heights in (9) can be bounded by `O_(E,D)(H+1)`. In fact, for the bounded degree furnished by the Nullstellensatz, the coefficient identities are a finite linear system whose entries are fixed coefficients of `I_d` or coefficients of `M`. At the particular basepoint-free specialization, choose a nonzero rank minor and solve. Its determinant and the requisite ratios of minors have heights linear in the height of `M`, with constants depending only on the bounded matrix size. This is why degeneration does not introduce an uncharged comparison constant.

At each place, (9) bounds the maximum coordinate size of `X` by the maximum of `(MX)_j` times a local coefficient bound. The direct linear equations give the reverse inequality. Summing with normalized local weights and using the product formula yields

```text
|h(M iota_d(X)) - h(iota_d(X))|
    <= C_(E,D)(H+1).                                 (10)
```

All fields involved have bounded degree, although the absolute-height argument itself uses normalized weights and is unchanged by a further field extension.

The fixed line bundle `O_E(d[O])` is symmetric. Its canonical height is `d q`, so the height comparison for this fixed complete linear system is

```text
h(iota_d(X)) = d q(X)+O_(E,d)(1).                     (11)
```

Only finitely many `d` occur. Combining (10) and (11),

```text
h(Psi(X)) = d q(X)+O_(E,D)(H+1).                     (12)
```

The smaller degree cases require no elimination:

- If `d=0`, all four functions are constants and `0 <= h(Phi) <= O_(E,D)(H+1)`, proving (1).
- A nonconstant morphism of this kind cannot have `d=1`: a degree-one line bundle on an elliptic curve has a one-dimensional section space, so no basepoint-free nonconstant tuple exists. Equivalently, a constant map has pullback degree zero.
- If `d=2`, the fixed complete series maps to `P^1`. Basepoint freeness forces `M` to have rank two; otherwise all sections would be proportional to a section of a positive-degree bundle, which has a zero. A nonzero two-by-two minor gives an inverse of coefficient height `O_(E,D)(H+1)`. This proves (10) directly, and (11) remains valid.

## 6. Expansion of the canonical height and the moving-map consequence

Take `X=Q-T` in (12). We obtain

```text
h(Phi(Q)) = d q(Q-T)+O_(E,D)(H+1).
```

Using the canonical pairing,

```text
q(Q-T) = q(Q)+q(T)-2<Q,T>,
|<Q,T>| <= sqrt(q(Q)q(T)).
```

Equation (5), together with the bounded range of `d`, proves (1) and hence (2).

In particular, fix a nontorsion rational point `P`, put

```text
h_N = q(NP) = N^2 q(P),
```

and allow a new rational frame map `Phi_N` at each index. If their presentation degrees are bounded by the same `D` and their coefficient heights satisfy `H_N=o(h_N)`, then

```text
h(Phi_N(NP)) >= (deg Delta_N) h_N - o(h_N).            (13)
```

Thus a separate argument forcing `deg Delta_N >= 25` would contradict an upper bound `h(Phi_N(NP)) <= 4h_N+o(h_N)`. The present lemma supplies precisely the uniform moving-map height comparison needed for that implication. It does not itself force the contacts or provide a bounded-degree, subpower-coefficient rational presentation of arbitrary pointwise frames.
