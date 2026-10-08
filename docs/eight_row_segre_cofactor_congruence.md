# A nonlinear cofactor congruence at an eight-row boundary

Four exact numerical invariant relations on eight actual Gaussian rows
give a cubic congruence in the maximal minors of their first-smaller-cut
restriction matrix. The congruence keeps full prime powers and needs no
rank hypothesis modulo a conductor prime. At a sufficiently small
coefficient height it becomes an exact Segre-cubic identity.

This strengthens the local necessary conditions beyond proper restriction
rank. It does not prove the endpoint bound: the resulting rational point
can lie on the boundary or singular locus, and no contradiction from the
simultaneous conditions on overlapping cuts is established here.

The restriction arithmetic is from
[the all-cut hierarchy](all_cut_invariant_relation_hierarchy.md).
The Segre cubic itself is already recorded in
[the six-row matching calculation](segre_gradient_arithmetic.md).
The additional statement here applies that cubic to the **cofactor vector
of numerical relations**, with an exact common-content loss.

## 1. The five-dimensional restriction image

For eight rows, let `K_1` be the degree-two-each invariant space vanishing
on every balanced cut. Fix a three-row inside set. Its restriction lies
in `W_(5,2)`, the squarefree translation-invariant quadrics in the five
outside variables `z_1,...,z_5`. An integral basis is

```text
A=(z_2-z_1)(z_4-z_3),
B=(z_2-z_1)(z_5-z_4),
C=(z_4-z_1)(z_3-z_2),
D=(z_3-z_2)(z_5-z_4),
E=(z_5-z_2)(z_4-z_3).
```

Its coefficient matrix on the five monomials
`z_1z_2,z_1z_3,z_1z_4,z_2z_3,z_2z_4` has determinant of absolute
value one. Since `dim W_(5,2)=5`, these matching polynomials form an
integral basis. Their values satisfy

```text
Phi(A,B,C,D,E)=B C E-A D(A+B+C+D+E)=0.                 (1)
```

One may check (1) by expansion, or append an infinity row to the five
outside binary rows and use the existing six-row matching identity.

For actual outside Gaussian rows `P_j=X_j+iY_j`, set `z_j=Y_j/P_j`
and multiply all five coordinates by `product_j P_j`. The resulting
Gaussian vector `v` still satisfies (1). Its coordinates are products
of two real brackets and one unused `P_j`, for example

```text
v_A = Delta_12 Delta_34 P_5,
Delta_ij=P_i Y_j-Y_i P_j=X_i Y_j-Y_i X_j.
```

All five coordinates are nonzero when the actual outside directions are
distinct. They need not have Gaussian gcd one.

## 2. Exact modular cofactor lemma

Let `R` be an integral `4 by 5` matrix, and let `c` be its signed
four-minor vector, so `R c=0`. Let `v` be a nonzero Gaussian integer
vector satisfying (1). Suppose, for Gaussian integers `H!=0,kappa!=0`,

```text
H divides every coordinate of kappa R v.
```

Put `G_v=gcd_G(v_1,...,v_5)`. Then the exact divisibility is

```text
H/gcd_G(H,kappa G_v) divides Phi(c).                   (2)
```

If `H` is conjugate-primitive, this gives the ordinary integer divisor

```text
M = Norm(H)/Norm(gcd_G(H,kappa G_v)),
M divides Phi(c).                                    (3)
```

**Proof.** Divide `v` by its full Gaussian gcd. The new vector `u` is
primitive and satisfies `Phi(u)=0`. Put
`H_0=H/gcd_G(H,kappa G_v)`; then `H_0 | R u`. Cofactor-adjugate
identities give, for all `i,j`,

```text
H_0 divides c_j u_i-c_i u_j.
```

Work at any prime power dividing `H_0`, and choose a coordinate `u_j`
which is a unit there. Modulo that prime power, `c` is a scalar multiple
of `u`. Homogeneity gives `Phi(c)=0` modulo the same full prime power.
This proves (2). A conjugate-primitive Gaussian divisor of an ordinary
integer contributes its full norm to that integer, proving (3).

No four-minor was inverted. If `R` drops rank modulo a prime, some or
all cofactors can acquire content, but the proof and (2) remain valid.
The loss involves **one** copy of the vector gcd `G_v`; using an arbitrary
un-normalized coordinate and then cubing would unnecessarily lose three.

The raw cofactor vector matters. If `q=gcd(c_1,...,c_5)` and `c!=0`,
the guaranteed divisor for the primitive vector is only

```text
M/gcd(M,q^3) divides Phi(c/q).                         (4)
```

For example let `R=5[I_4 | (1,2,3,4)^t]`, take the primitive Segre
point `v=(1,1,3,1,3)`, and let `H=2+i`, `kappa=1`. The hypotheses
hold, while the primitive cofactor vector is `(-1,-2,-3,-4,1)` and
its cubic value is `42`, a 5-unit. The raw vector has common factor
`5^4`, and its cubic has the required divisibility. This is an explicit
rank-drop example, not an exception to (2).

## 3. Complementary conductor blocks and the height threshold

Let the rows of `R` be the five coefficients of four integral numerical
zero relations in `K_1`, restricted at the selected inside set `S`.
The actual Gaussian congruence of the all-cut hierarchy is

```text
H_S | kappa R v,       kappa=product_(j outside S) K_j.
```

Therefore (3) applies. The common coordinate gcd divides every selected
matching coordinate. Hence its loss is no larger than the loss already
audited for one matching coordinate in that hierarchy. In its notation,
for `m=8`, correction bound `log|K_i|<=sigma`, and primitive bracket
residue bound `T`,

```text
log M_S >= log Norm(H_S)-16 sigma-2 log T.              (5)
```

Apply the existing complement reflection to the same four numerical
relations. Their restriction coefficient matrix is unchanged, and the
reflected five-coordinate vector again lies on (1). Thus a second modulus
`M'_S`, supported on `Norm(H_(S^c))`, divides the **same** integer
`Phi(c)`. The two moduli are coprime. The already audited correction
losses give

```text
M_S M'_S | Phi(c),
log(M_S M'_S) >= log Norm(H_S)+log Norm(H_(S^c))
                -272 sigma-4 log T.                  (6)
```

This multiplication does not assume that any cofactor or matching value
is a unit at every conductor prime.

If each row of `R` has coefficient sum at most `L`, Hadamard gives
`|c_i|<=L^4`. The cubic has six signed monomials, so

```text
|Phi(c)| <= 6 L^12.                                  (7)
```

For full core blocks with log norms at least `(1-eta)w`, it follows that

```text
12 log L+log 6 < 2(1-eta)w-272 sigma-4 log T
          implies Phi(c)=0.                          (8)
```

If `rank_Q R=4`, then `c!=0`, and its rational projective kernel point
lies exactly on the Segre cubic. If the rank is smaller, all cofactors
vanish and this particular test is vacuous. In the limiting full-profile
regime the coefficient threshold is `log L<(1/6-o(1))w`.

The missing next step is now explicit. One must either produce sufficiently
short numerical relations whose restriction has rank four and whose
cofactor point is off the Segre cubic, or prove that the recovered cubic
points at many overlapping cuts cannot be simultaneously compatible with
the full actual Gaussian data. Neither assertion is established here.
An exact point of the Segre cubic, especially a singular or boundary point,
does not itself recover five distinct rational outside directions. The
existence of four relations at the required height is also not supplied
by this lemma. These limitations prevent claiming any new radius exponent.

## Verification

`python3 -B docs/check_eight_row_segre_cofactor_congruence.py` checks the
integral matching basis, the universal polynomial identity (1), exact
prime-power and multi-prime instances of (2)--(3), and the primitive
cofactor counterexample in (4). The modular fixtures use actual matching
values; they are not asserted to be endpoint configurations or restrictions
of global numerical zero relations.
