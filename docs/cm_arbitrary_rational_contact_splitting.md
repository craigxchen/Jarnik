# CM splitting for arbitrary finitely many rational contacts

The finite-shift construction in
[rational_cm_inert_mass.md](rational_cm_inert_mass.md) does not require
all contact points to be integer multiples of one rational point. An
abelian-group Chinese remainder argument extends it to arbitrary fixed
rational contacts on a fixed curve with the same Gaussian CM action.
This supplies contact norms and Gaussian allocations, not their common
oriented lift to the prescribed rational configuration.

## Statement

Let `E/Q` be a fixed elliptic curve with the full `Z[i]` action defined
over `Q(i)`, with complex conjugation sending `[i]` to `[-i]`. For
example, `y^2=x^3-cx`, with `c!=0`, has `[i](x,y)=(-x,iy)`.
Assume `E(Q)` contains a nontorsion point `P`. Fix distinct rational
points `B_1,...,B_r` and any fixed finite set of exceptional primes.
Use one fixed integral Weierstrass model, and write `b(Q)` for the
positive square-root denominator of `x(Q)`.

There are rational points `Q_j` of heights `H_j=hat h(Q_j)` tending
to infinity such that, simultaneously for all contacts,

```text
log b(Q_j-B_e)=H_j+o(H_j),
sum_(p=3 mod 4) v_p(b(Q_j-B_e)) log p=o(H_j),
max_p v_p(b(Q_j-B_e)) log p=o(H_j).                    (1)
```

The fixed exceptional-prime parts are `exp(o(H_j))`. Outside those
primes and the finitely many primes where the sections `B_e` intersect,
the contact denominators have pairwise disjoint prime support. After
discarding a total `o(H_j)` logarithmic mass at each contact, each
retained denominator can be split into any fixed finite list of
positive prescribed proportions by whole Gaussian denominator cells.
The resulting Gaussian factors and all their conjugates have disjoint
support, across contacts as well as within one contact.

The sequence need not lie in one fixed cyclic subgroup or one fixed
translate of it. This is the additional freedom supplied by the
abelian-group CRT below.

## CRT in the rational point group

Fix a tolerance `epsilon>0`. Choose pairwise coprime squarefree products
`m_e` of odd split primes such that

```text
rho(m_e)=product_(p|m_e)(1-2/p+2/p^2)<epsilon.
```

They can avoid any prescribed finite set. The existence follows from
the split-prime reciprocal divergence used in the earlier CM note.
Choose integers `c_e` with `c_e=1 mod m_e` and `c_e=0 mod m_f` for
`f!=e`. Put

```text
M=product_e m_e,        Q_0=sum_e [c_e]B_e,
S_e=sum_f [(c_f-1_(f=e))/m_e]B_f.
```

All the coefficients defining `S_e` are integers, and exactly

```text
Q_0-B_e=[m_e]S_e.
```

For the stage parameter `n`, set

```text
Q(n)=[Mn]P+Q_0,
R_e(n)=[(M/m_e)n]P+S_e,
Q(n)-B_e=[m_e]R_e(n).                                 (2)
```

No basis of the Mordell--Weil group, division of a rational point, or
assumption that the contacts lie in one cyclic subgroup is needed.

## All omitted local terms can be bounded at each fixed stage

For each `e`, let `J_e` index the Gaussian primes over `m_e`, and put
`alpha_S=product_(pi in S)pi` for `S subset J_e`.
The finitely many relevant maps are

```text
z |-> [alpha_S]([(M/m_e)]z+S_e).
```

Enlarge a finite set `Sigma` of rational primes to contain the exceptional
primes, all bad model and CM-action primes, all primes dividing `2M`,
and all primes where different contact sections meet. Omit every
Gaussian prime above each prime in `Sigma`, so that the omitted set is
invariant under conjugation. Include both complex embeddings in the
archimedean bounds. Choose an integer `n_0` for which none of these
finitely many maps sends `[n_0]P` to the identity. This excludes only
finitely many integers: the coefficient multiplying `n_0P` is a
nonzero isogeny and `P` is nontorsion.

Compact recurrence of the diagonal cyclic orbit in
`E(R) times product_(p in Sigma) E(Q_p)` gives infinitely many positive
`n` returning near `[n_0]P`. Choose the neighborhoods sufficiently small
that every displayed map stays away from the identity at all relevant
places, including after the Gaussian CM maps at complex or extended
finite places. Their local denominator heights are then bounded.
This replaces any need for a logarithmic bound on approximation to a
fixed noncyclic translate. All bounds here are allowed to depend on
the fixed stage.

## The unchanged Boolean cell calculation

For a fixed `e`, form the denominator ideals of `[alpha_S]R_e(n)` over
`Q(i)`, omitting `Sigma`, and their Boolean differences `C_T` exactly as
in the earlier note. Its local argument applies to every rational
point `R_e(n)`: the annihilator at a denominator prime divides `(m_e)`,
full denominator depth is constant on the upper Boolean interval, and
the resulting cells are integral and pairwise coprime. Frobenius at an
inert prime forces its cell label to be conjugation-stable.

The bounded omitted local terms and canonical-height functoriality give

```text
log Norm C_T
 =2 hat h(R_e(n)) product_(pi in T)(Norm(pi)-1)+O_stage(1).
```

Summing the stable cells and bounding a single cell therefore give

```text
inert_mass(b(Q(n)-B_e))
 <= rho(m_e) hat h(Q(n)-B_e)+O_stage(1),

max_p v_p(b(Q(n)-B_e)) log p
 <= 2 L(m_e) hat h(Q(n)-B_e)+O_stage(1),
L(m_e)=product_(p|m_e)(1-1/p)^2 <= rho(m_e).             (3)
```

The factors two are the same as in the earlier norm calculation: a
rational denominator has logarithmic Gaussian ideal norm twice its
rational logarithm. The maps with `S=J_e` include `Q(n)-B_e`, so the
same local bounds also give

```text
log b(Q(n)-B_e)=hat h(Q(n)-B_e)+O_stage(1).
```

As `n` tends to infinity within the fixed-stage subsequence,

```text
H(n)=hat h(Q(n))=M^2 n^2 hat h(P)+O_stage(n),
hat h(Q(n)-B_e)=H(n)+O_stage(n).                       (4)
```

Use tolerances tending to zero and choose one sufficiently large
return time at each stage so that all its errors divided by `H(n)`
tend to zero. This proves (1). The stable-cell removal and cumulative
allocation of whole nonstable conjugate pairs then give the stated
proportional Gaussian factors exactly as in Section 6 of the earlier
note. No stage-uniform error estimate has been assumed.

## Consequence and remaining scope

If a fixed CM elliptic curve of this type has a birational moduli map
whose boundary pullbacks are distinct reduced rational points, the argument
supplies balanced contact denominators with negligible inert mass and
no dominant prime power for all those contacts simultaneously. Fixed
integral-model comparison transfers them to the boundary contact
contents with only fixed exceptional-prime changes.

It does not prove that such a CM six-point moduli curve exists. The
verified disjoint six-curve currently in the project has not been
identified with a Gaussian CM curve. Even if a CM example is found,
the theorem does not recover its rational moduli point by multiplying
the chosen Gaussian blocks. Compatible orientations, private factors,
small original Gaussian residues and the short determinant-one frame
remain additional requirements.

The [independent proof audit](cm_arbitrary_rational_contact_splitting_audit.md)
checks the compact recurrence, full prime-power depths, height
normalization, and diagonalization.
