# The collective phase field of rational CM cells

This note concerns fixed auxiliary products in
[the rational CM denominator construction](rational_cm_inert_mass.md).
Its new statement is a function-field theorem: the phase functions of
whole nonstable Boolean cells, grouped into distinct complementary cut
pairs at distinct rational contacts, collectively generate `Q(i)(E)`.
It does not establish small residues for the associated rational moduli
points or uniform estimates as the auxiliary products vary.

In particular, the twenty blocks assigned to the ten boundary contacts
of the five-row model omit ten peripheral singleton/quadruple blocks.
Those omitted blocks have substantial profile height. Their phases
cannot be absorbed into a correction of height `o(N^2)`.

## 1. Primitive torsion divisors

Let `E:y^2=x^3-2x`, `P=(2,2)`, and `K=Q(i)`. For a fixed squarefree
product `m` of odd split rational primes, let `J` index its Gaussian
prime factors, including both conjugate orientations. Set
`alpha_S=product_(pi in S) pi`.

The divisor

```text
Z_T=sum_(S subset T) (-1)^(|T|-|S|) [alpha_S]^*[O]       (1)
```

is effective: it is exactly the sum of the torsion points whose
annihilator ideal in `Z[i]` is `(alpha_T)`. This follows by applying
Boolean inversion to each point of `E[m]`. In particular

```text
deg Z_T=product_(pi in T)(Norm(pi)-1),
Z_empty=[O],
Z_T and Z_U have disjoint supports if T!=U.             (2)
```

These divisors are defined over `K`, because the CM endomorphisms are
defined over `K` and Galois preserves their annihilators. For nonempty
`T`, their supports have odd order and are stable under negation, with
no fixed point of negation. Their group sum is therefore `O`.
Consequently

```text
Delta_T=Z_T-Z_(bar T)
```

has degree zero and group sum `O`, so it is principal. If `T` is not
conjugation stable, choose the function `f_T` with this divisor and
`f_T(O)=1`. The origin is outside its support. The normalization makes
the function defined over `K` and gives

```text
bar(f_T)=1/f_T,       f_(bar T)=1/f_T.                  (3)
```

Here conjugation acts on coefficients and on `i`; on real points the
function has complex modulus one. Stable labels give the trivial
divisor and may be represented by the constant function one.

## 2. Local audit of the cell phase

Let `R=tP` and let `C_T(R)` be the denominator cell after removing all
primes over `2m`, exactly as in the denominator note. Choose a Gaussian
generator `c_T(R)` and put
`gamma_T(R)=c_T(R)/bar(c_T(R))`.

For every prime `v` away from `2m`, the curve has good reduction and
each `[alpha_S]` is etale. Pullback of the identity section therefore
identifies the intersection of `R` with its kernel divisor with the
full formal denominator depth of `[alpha_S]R`. Applying (1) gives
the intersection with `Z_T`. The horizontal divisor of `f_T` then
identifies the difference of the `T` and `bar T` depths with
`v(f_T(R))`, up to a vertical constant independent of `R`.

That constant is zero. All nonzero `m`-torsion sections are disjoint
from `O` at residue characteristics prime to `m`; compare with the
section `O`, where `f_T(O)=1` and both horizontal intersections
vanish. Thus, away from `2m`,

```text
v(gamma_T(tP))=v(f_T(tP)).                             (4)
```

The same reasoning may be made after a finite extension splitting the
torsion points, and descends because both sides are defined over `K`.
Distinct torsion sections of order prime to the residue characteristic
remain distinct, so no collision is being silently discarded here.

At the finitely many omitted primes, the discrepancy has valuation
`O_m(log(2+|t|))`. For an explicit bound on the function side, pullbacks
of the identity under the finitely many fixed maps `[alpha_S]` have
horizontal torsion divisors and fixed vertical corrections on a fixed
regular model. Their intersection depths at `tP` are bounded by the
denominator depths of `t[alpha_S]P`, plus constants. The latter are
`O_m,v(log(2+|t|))` by the formal-group logarithm after passage to a
fixed deep open subgroup. The shallower depths and vertical components
contribute only constants. Taking the finite alternating sums in (1)
preserves this bound. The truncated cell has zero valuation at these
primes, by definition.

It follows that the exact phase law has the form

```text
gamma_T(tP)=s_T(t) f_T(tP),
s_T(t) in K*,      s_T(t) bar(s_T(t))=1,
h(s_T(t))=O_m(log(2+|t|)).                             (5)
```

The factor `s_T(t)` is supported on the fixed primes over `2m`, apart
from Gaussian units. The norm-one condition makes its archimedean
modulus one and excludes any valuation at the prime over two. If a
fixed additional exceptional set is trimmed, its fixed-prime depths
obey the same logarithmic bound and can be included in `s_T`.

This supplies the needed local phase audit. It is stronger than an
identity of radicals, and it does not assert that the omitted-prime
phase factor belongs to a finite set.

## 3. Several contacts and whole-cell allocations

Fix at least two distinct integer contact labels `k_e`. Choose pairwise
coprime auxiliary products `m_e>1`, each as above, and choose `N_0` with
`N_0=k_e mod m_e`. Put

```text
M=product_e m_e,
N=N_0+Mj,
a_e=(N_0-k_e)/m_e,       b_e=M/m_e,
R_e=a_e P+[b_e]z,       z=jP.                          (6)
```

For each contact choose one cell from every nonstable conjugate pair
`{T,bar T}`, and divide the chosen cells into two groups. Assign the
groups to complementary nontrivial Boolean cut labels `b` and
`bar b`. Different contacts must receive different unordered cut pairs.
All cells are used. The groups may, for example, be those giving
balanced complementary weights in the denominator note.

For fixed auxiliaries there are only finitely many such orientation and
grouping choices. Any infinite sequence of allocations therefore has
a subsequence with fixed choices. Define the row-to-anchor phase
functions by multiplying

```text
f_T(a_eP+[b_e]z)
```

with the exponents specified by those row memberships. These are fixed
functions `F_1,...,F_(r-1)` on `E` over `K`, where `r` is the number
of row labels. Formula (5) expresses the actual cell-product phases
at `z=jP` as these functions times factors of height `O_m(log j)`.

For a cut `b`, its vector of exponents relative to row `r` is
`c_b=(b_1-b_r,...,b_(r-1)-b_r)`. Distinct nontrivial complementary
pairs give distinct vectors up to sign. Thus every nonzero valuation
vector in the collective map is primitive and has the form `±c_b`.
It occurs only in the one contact assigned to that complementary pair.

## 4. The collective map is birational

Choose `Q_e in E(Qbar)` satisfying `[b_e]Q_e=-a_eP`.
All zeros and poles belonging to contact `e` lie in

```text
Q_e+E[M].                                             (7)
```

Different contacts have nontorsion differences between these centers:

```text
[b_e b_f](Q_e-Q_f)
 =M(k_e-k_f)/(m_e m_f) P != O modulo torsion.           (8)
```

Let `pi:E->C` be the map to the smooth projective curve with function
field `K(F_1,...,F_(r-1))`, and write `d=deg(pi)`. The curve has genus
zero or one by Riemann--Hurwitz. At any point supporting a nonzero
valuation vector, primitivity forces ramification index one. Every
other point in its fiber has the same valuation vector and hence
belongs to that same contact coset in (7).

The image cannot have genus zero. If it did, all fibers would be
linearly equivalent divisors of degree `d` on `E`. A fiber in contact
`e` has group sum `d Q_e` plus torsion; a fiber in another contact
`f` has group sum `d Q_f` plus torsion. Linear equivalence would make
these sums equal, contradicting (8).

Thus `C` has genus one and `pi` is an isogeny followed by a translation.
Its kernel translations preserve all the phase functions and their
valuation vectors. They cannot permute distinct contact cosets, by
(8). It remains to show that a translation preserving each contact's
support belongs to `E[b_e]`.

Under the CM decomposition

```text
E[m_e] = product_(q|m_e)(F_q x F_q),
```

the full support of nonstable cells is the complement of

```text
B_e=product_(q|m_e) ({(0,0)} union (F_q* x F_q*)).       (9)
```

Indeed the annihilator support is conjugation stable exactly when, for
every `q`, both coordinates vanish or both are nonzero. Each factor
in (9) has `q^2-2q+2` points, which is `2 mod q`. A nonzero additive
translation of `F_q^2` has order `q`, so any subset invariant under
it has cardinality divisible by `q`. Each factor therefore has trivial
translation stabilizer. Projection to the factors gives the same
conclusion for `B_e`, and hence for its complement.

The contact support upstairs is the inverse image of that complement
under `[b_e]`, translated by `Q_e`. A translation preserving this
nonempty support must have `[b_e]h in E[m_e]` and must stabilize the
complement in (9). Thus `[b_e]h=O`. Any common kernel is consequently
contained in every `E[b_e]`. Since the `m_e` are pairwise coprime,
`gcd_e b_e=1`; that common kernel is trivial. We have proved

```text
K(F_1,...,F_(r-1))=K(E).                               (10)
```

This argument uses the entire nonstable support of every contact.
It need not remain true if arbitrary cells are omitted or if individual
prime factors inside a cell are reoriented independently.

## 5. Exact scope of a height-compression consequence

For fixed auxiliary products and fixed allocations, (10) expresses
`x(z)` as a fixed rational function of the collective phases. Suppose
these particular phase functions at `z=jP` were all within
`exp(-(2hhat(P)+epsilon)j^2)` of a tuple in `K` of height `o(j^2)`,
for some fixed `epsilon>0`. The collective image equations first
force that tuple onto the image, by the product formula. The rational
inverse then supplies a target for `x(jP)` of height `o(j^2)`, with
the same exponential accuracy up to `exp(o(j^2))`. The Gaussian
numerator estimate in
[the phase audit, Section 4](cm_denominator_phase_lift_audit.md)
contradicts `h(x(jP))=2hhat(P)j^2+O(1)`. Finitely many inverse
exceptional targets are covered by the fixed-target elliptic-logarithm
bound. Thus this explicitly stated simultaneous approximation is
impossible.

This implication applies to the whole-cell product phases just defined.
It applies to actual full row phases only if every additional
nontrivial cut block is also included in the fixed phase description,
or its remaining phase contribution is a target of height `o(j^2)`.
The ten five-row boundary contacts alone do not meet that condition:
their ten peripheral singleton/quadruple blocks each have substantial
height. One could instead describe those cut pairs by additional
distinct rational contacts, but that is a further specified construction,
not a consequence about arbitrary lifts of the ten boundary contacts.

All function degrees, coefficient heights, inverse constants, and
exceptional targets in this note depend on the fixed auxiliary products
and the chosen contacts. The denominator theorem uses a diagonal
sequence with growing products. No uniform control of these phase
constants along that diagonal has been established here.

The [finite support checker](check_rational_cm_cell_phase_field.py)
verifies the stable-support counts and translation stabilizers for
the split primes five and thirteen, their product, and the Boolean
valuation-vector incidence. Product stabilizers are checked factor
by factor; the total product cardinality need not be nonzero modulo
every constituent prime. The divisor and function-field arguments
above are separate from those finite checks.
