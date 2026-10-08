# Saturation of normalized local generators at a simple coalition

Let `R` be a DVR with uniformizer `pi` and fraction field of
characteristic zero. Choose any nonzero nonunit `s in R`; in
particular `s=pi^e` is allowed for every `e>=1`.
Put `m=nd`, `n,d>=2`, and partition the labels
into a coalition `T` and its complement `O`. In an affine binary
chart suppose

```text
z_i=s a_i       (i in T),
z_j=b_j         (j in O),
```

where the `a_i,b_j` belong to `R`, every difference `a_i-a_k`
within `T` is a unit, every difference `b_j-b_l` within `O` is a
unit, and every `b_j` is a unit. Higher-order terms are allowed in
the `a_i,b_j`. These assumptions describe one simple coalition
whose two residual configurations have distinct directions.

For `|I|=2n`, set

```text
L_I=det[V|diag(z)V]_I,
g_I=max(0,|I intersection T|-n),
Lbar_I=s^(-g_I)L_I.
```

The coefficients of `Lbar_I` belong to `R`. Let `K_R` be the
saturated lattice of multilinear invariant polynomials in the
generic coherent kernel, and let `N_R` be the `R`-span of

```text
Lbar_I times products of d-2 ordinary n-row determinants
on a partition of the complementary labels.
```

Then `N_R` has uniformly bounded index in `K_R`. More precisely,
put `q=n(d-2)`, let `lambda=((d-2)^n)`, and let

```text
C_(n,d)=q!/f^lambda,
```

the hook product of that rectangle, with `C_(n,2)=1`. Then

```text
C_(n,d) K_R subset N_R subset K_R.                    (1)
```

In particular, if `C_(n,d)` is a unit in `R`, the normalized local
generators span the saturated kernel exactly. At `(n,d)=(3,4)`,

```text
C_(3,4)=6!/5=144.                                    (2)
```

Thus exact saturation holds in every residue characteristic except
possibly 2 and 3. At those two characteristics the quotient is
killed by 144, so its length is at most `341 v_R(144)`. There is
no exceptional algebraic relation among otherwise distinct
residual directions in this statement. The bound is independent
of `v_R(s)`, not just uniform for collisions of depth one.

## 1. The normalized ideal is a linked determinantal ideal

In `S=R[v_(i,a)]`, order the coalition rows first and define

```text
A_s = [ V_T        s diag(a) V_T ]
      [ V_O          diag(b) V_O ],

B_s = [ V_T          diag(a) V_T ]
      [ s V_O        diag(b) V_O ].
```

For a `2n`-row set `I` containing `k` coalition rows, row and column
scaling over the fraction field gives

```text
det(B_s)_I = s^(n-k) det(A_s)_I.                      (3)
```

Consequently the ideal `J=(Lbar_I)` is exactly

```text
J=I_(2n)(A_s)+I_(2n)(B_s).                            (4)
```

Indeed, for `k<=n` its generator is the `A_s` minor and the
corresponding `B_s` minor is its multiple by `s^(n-k)`; for
`k>=n` use the `B_s` minor instead.

This is a two-term linked determinantal locus in the sense of
[Murray–Osserman, Definitions 1.1–1.2 and Theorem 1.3](https://arxiv.org/pdf/1412.3818).
Explicitly, the two source bundles have rank `2n`; their linking
maps are `diag(s I_n,I_n)` and `diag(I_n,s I_n)`. Their composites
are `s` times the identity, and their kernels and images agree
as required on the special fiber. The endpoint maps are

```text
g_1=[V_T | diag(a)V_T],       g_2=[V_O | diag(b)V_O].
```

The two induced maps to the direct sum of the endpoint targets
are precisely `B_s` and `A_s`. The rank bound is `2n-1`, so the
expected codimension is

```text
c=(2n-(2n-1))(m-(2n-1))=m-2n+1.                     (5)
```

The cited theorem says that every component has codimension at
most `c`, and that expected codimension over a Cohen–Macaulay
base gives a Cohen–Macaulay locus. The nonempty proper coalition
case fits its endpoint-rank hypotheses. Empty or full coalitions
reduce directly to the ordinary determinantal ideal after a
common column scaling.

## 2. Special-fiber dimension for every allowed residue tuple

Work over an algebraic closure `k` of the residue field and write

```text
A_0 = [ V_T     0             ],
      [ V_O     diag(b) V_O   ]

B_0 = [ V_T     diag(a) V_T   ].
      [  0      diag(b) V_O   ]
```

Let `X_0` be the locus where both matrices have rank at most
`2n-1`, and put `D=mn-m+2n-1=mn-c`. Three incidence charts cover
`X_0` set-theoretically.

* If `ker A_0` has a vector `(u,w)` with `u!=0`, use the projective
  pair `[u:w]` and the equations

  ```text
  V_T u=0,       V_O u+diag(b)V_O w=0.                (6)
  ```

  They automatically put `(u,0)` in `ker B_0`.

* If `ker B_0` has a vector `(u',w')` with `w'!=0`, use `[u':w']`
  and the equations

  ```text
  V_T u'+diag(a)V_T w'=0,       V_O w'=0.             (7)
  ```

  They automatically put `(0,w')` in `ker A_0`.

* Otherwise the two nonzero kernels are supported in their
  opposite coordinate blocks. Choose independent projective
  parameters `[u']` and `[w]` (the vectors need not be linearly
  independent in `k^n`) and impose

  ```text
  V_T u'=0,       V_O w=0.                           (8)
  ```

For (6), all coalition row normals are nonzero because `u!=0`.
If `u,w` are linearly independent, all outside row normals are
also nonzero. Its incidence has dimension
`(2n-1)+m(n-1)=D`. The dependent-pair base has dimension at most
`n`, with all row equations present except possibly at one
outside label. Distinctness of the `b_j` is exactly what ensures
at most one missing equation. At such a missing-equation stratum
the base dimension falls to `n-1`, cancelling the added free row
dimension. Thus the dependent incidences have dimension at most
`m(n-1)+n<D`.

The same argument applies to (7), using distinctness of the
`a_i` and the fact that its outside row normals `w'` are nonzero.
Finally, (8) has base dimension `2n-2` and all `m` independent
row equations, so its incidence dimension is `D-1`. Therefore

```text
dim X_0 <= D                                        (9)
```

for every residual configuration satisfying the unit assumptions,
in every residue characteristic. No generic-residue inference is
used here.

## 3. Flatness and saturation of the polynomial ideal

On the generic fiber, all slopes `s a_i,b_j` are distinct. The
[determinantal generation theorem](distinct_configuration_determinantal_kernel_generation.md)
gives codimension `c` for the common generic ideal of `A_s,B_s`.
Its closure therefore has dimension `D+1` in `Spec S`.

By (9), a component supported entirely in the special fiber would
have dimension at most `D`, hence codimension at least `c+1`.
The linked determinantal codimension bound rules this out. Thus
every component dominates the DVR and has codimension `c`.
The linked determinantal Cohen–Macaulay theorem now applies.

There are no embedded associated primes, and `pi` belongs to no
minimal prime. Hence `pi` is a nonzerodivisor on `S/J`; equivalently,
`S/J` is flat over `R`. In particular,

```text
J = J[pi^(-1)] intersection S.                       (10)
```

For `Q in K_R`, the all-distinct theorem on the generic fiber
puts `Q` in `J[pi^(-1)]`. Its coefficients already belong to `R`,
so (10) gives an integral presentation

```text
Q=sum_I Lbar_I F_I,       F_I in S.                  (11)
```

Row multigrading makes `F_I` independent of the rows in `I` and
multilinear on `I^c` without introducing any denominator. This
is saturation of a polynomial ideal; it still remains to extract
the invariant coefficient polynomials.

## 4. A fixed denominator for the invariant projection

The coefficient polynomials in (11) lie in `q=n(d-2)` tensor
positions. Over the fraction field, the `SL_n` invariant part is
the `GL_n` representation `det^(d-2)`. Schur–Weyl duality identifies
the projection to this part with the symmetric-group central
idempotent

```text
P_lambda = (f^lambda/q!) sum_(sigma in S_q)
              chi_lambda(sigma^(-1)) sigma,
lambda=((d-2)^n).                                    (12)
```

The group acts by permuting the complementary tensor positions.
Its irreducible characters are integer-valued. Therefore
`C_(n,d) P_lambda` has integral matrix entries, with
`C_(n,d)=q!/f^lambda`. The central-idempotent formula is recorded,
for example, in
[Ram, equation (1.2)](https://soimeme.org/~arunram/Publications/1991InvMatv106p461.pdf).

To justify applying different complementary-position projections
in one equation, let `P_m` be the invariant projection on the full
multilinear degree-`m` module. Multiplication by `Lbar_I` is an
`SL_n`-equivariant map from the complementary tensor module to
that full module; it carries the `GL_n` twist `det^2`, which is
trivial on `SL_n`. Naturality of invariant projection gives

```text
P_m(Lbar_I F_I)=Lbar_I P_lambda(F_I).
```

The symmetric-group operator in (12) acts only on the complement
of this particular `I`; it need not be one common permutation
operator on all `m` labels. Apply `P_m` to (11), use `P_m Q=Q`,
and multiply by `C_(n,d)`. This gives

```text
C_(n,d) Q = sum_I Lbar_I · (C_(n,d) P_lambda F_I).    (13)
```

The new coefficients are integral invariant polynomials.
They are integral linear combinations of standard determinant
products. One can see the absence of another lattice denominator
directly: standard tableaux give a characteristic-zero basis,
and their determinant products have distinct diagonal leading
monomials, each with coefficient one. With these leading
monomials ordered, the corresponding coefficient minor is
triangular with diagonal entries one. Thus their span is a
saturated submodule of the full monomial lattice.

Equation (13) proves (1). At twelve rows `lambda=(2,2,2)` and
`f^lambda=5`, giving (2). The five six-row gradient invariants
used in the existing generator matrix differ from a determinant
product basis by an integral unimodular change, so the same
bound applies to that displayed 4620-column family.

## 5. Scope

The theorem concerns a single coalition whose normalized residual
brackets are all units. It supplies an exact local spanning theorem
away from a fixed set of residue characteristics, and a fixed
annihilator for the remaining local quotient. It does not apply
unchanged when residual directions themselves collide: that case
requires another degeneration or an estimate in the residual
brackets. It also does not by itself compare global Archimedean
covolumes or bound the cost of joining presentations at different
primes.

The subsequent [residual-discriminant estimate](residual_discriminant_local_index_bound.md)
removes the unit-residual restriction at a cost bounded by a fixed
power of that discriminant. The resulting
[global index theorem](global_primitive_local_generator_index.md)
controls the index of the joined primitive generators under the
full profile. Neither result is an Archimedean shortest-vector bound.

The existing optional audit
`python3 docs/check_generic_terminal_coherent_syzygies.py --extended-core`
checks column-normalized local generators at every coalition size
`0,...,12`; all thirteen generator ranks are 341 over `F_101`.
These are distinct from the row-normalized coherent-map ranks.
The linked determinantal argument above gives the missing
all-residue theorem under the stated unit hypotheses.

The small [linked-ideal checker](check_linked_determinantal_coalition_saturation.py)
checks 945 minor-scaling and normalized-limit identities and the
integral determinant/gradient lattice changes. An independent
Gram-matrix computation finds denominator 72 for the six-row
invariant projection, which in particular verifies the safe
bound 144 used above.
