# Handoff to Claude: uniform lattice-point counts at exponent 1/2

Updated 2026-10-07 for the user's Claude review request, through item
533 of the research status. The latest checked results and exploratory
leads are distinguished below.

Please independently review the arguments below and then pursue the main
goal. Treat the recorded proofs as mathematical claims to audit, not as
axioms. A fresh route is welcome; there is no requirement to continue the
latest approach if it does not address the remaining quantifiers.

## Quick start for the reviewer

1. Open `/Users/cxc/Github/Jarnik`, or clone the repository and check out
   `agent/uniform-gcd-formalization`. Preserve any existing working-tree edits.
2. Read the exact goal and the unrestricted formulations below. Start
   with `inert_prime_cofactor_bound.md`,
   `exact_nested_profile_residual_extraction.md`, and
   `integer_cotangent_lcm_height_target.md`. The older
   `endpoint_full_profile_quantifiers.md` supplies the alternative
   prime-disjoint extraction; the template theory is optional unless
   you choose to build on it.
3. Independently audit the reduction you intend to use. Report a gap
   before relying on the claim. Then seek an argument that controls
   arbitrary circles and improves dependence on the radius.
4. If pursuing the template branch, audit the packet-projection proof,
   especially preserved coefficients and the kernel of the complete
   composition. A uniform template count still needs a transfer to
   arbitrary configurations.

The user welcomes a fresh approach and symmetry arguments, including
working backward from a sufficient contradiction. They do not require
following the most recent route, a sharp constant, or a formalization
before a correct mathematical proof. Return an independent review plus
the strongest new result you can actually justify. If the full theorem
remains open, state that plainly and identify the exact missing lemma.

### Suggested first review deliverable

Before extending a route, write a short audit identifying which claims
you checked, any gap or correction, and which conclusions depend on an
unproved input. In particular, audit the simultaneous tuple selection,
exact Gaussian gcd identities, and aggregate residual bound in item 532.
Then either prove a new arithmetic lower bound with all uniformity
parameters explicit, or explain precisely why your proposed method does
not yet supply one. No classification of arbitrary configurations by
the available templates has been proved.

For orientation, the latest reduction has
`H_sum <= k(k-1)W/[4(M-1)]`, where `W=log(R^2)` after primitive
normalization, `k` is a fixed selected tuple size, and `w=W/2^(k-1)`.
Here `M,R` are measured after retaining one common Gaussian-unit class
and dividing that class by its Gaussian gcd. The class retains at least
one quarter of the points, and its primitive radius does not exceed
the original radius. This differs from the normalization in the older
prime-disjoint extraction described below.
Its blocks permit nested prime overlap. A lower bound
`H_sum >= f(w)` with `f(w)/log w -> infinity` would improve the current
growth order; a fixed positive linear bound in `w` would prove
uniformity. These are sufficient targets, not established results.
The exact hypotheses and thresholds are in the linked extraction note.

## Repository and working state

- Local checkout: `/Users/cxc/Github/Jarnik`.
- GitHub remote: <https://github.com/craigxchen/Jarnik>.
- Research branch: [`agent/uniform-gcd-formalization`](https://github.com/craigxchen/Jarnik/tree/agent/uniform-gcd-formalization).
- The 2026-10-07 research snapshot includes the notes, checkers, certificates,
  Lean additions, and this handoff through item 533. Use this branch rather
  than the repository's default branch. The snapshot is integrated with
  the research already on the remote; those additional notes are a parallel
  line of investigation, not a proof of the uniform theorem.
- For provenance, commit `c2b8337` saves the local research before integration;
  the merge incorporates the prior remote history through `a2664da`.
- Three preexisting, untracked GitHub Actions templates from June 2026 are
  excluded: `create-release.yml`, `lean_action_ci.yml`, and `update.yml`.
  They add release, documentation-publication, or dependency-update actions
  unrelated to the research snapshot. The tracked Lean workflow is retained.
- `GaussianChain/` contains Lean, `src/numerics/` contains Rust numerical
  tools, and `docs/` contains the proof notes and exact Python checkers.

The main index is [endpoint_continuation_status.md](endpoint_continuation_status.md).
It is long: start with its opening status, items 514--533, and the final
continuation entry. Follow targeted links rather than reading every old
route in order. The older [quantitative handoff](codex_handoff_quantitative_uniformity.md)
concerns a previous Subspace Theorem route, not the current research state.

### Unfinished follow-ups at the snapshot

These are exploratory directions, not additional proved results:

- Extend item 533 to polynomial rows with variable imaginary parts. Clearing
  denominators gives exact pair factorizations, but a residual-degree gap
  still needs a global error bound. An extra section at the first center
  of a repeated-root collision can change inertia labels farther along
  the resolution, so a naive cost per bad root is insufficient.
- Extend the Kummer calculation to nested collision trees. Local blowup
  bookkeeping is available, but final branch-node counts and centers with
  two incoming branch exceptional curves still need a global estimate.
- Transfer the geometric obstruction to arithmetic profiles. The audit
  has not established an arithmetic analogue of the required Chern-number
  inequality: archimedean contributions and moving-family constants remain
  uncontrolled. This does not give a lower bound for the residual height.

On 2026-10-07, `lake build` completed successfully both before and after
remote integration (with linter warnings; the merged build had 8,514 jobs).
The merged `GaussianChain/ChordRatio.lean` also passed a separate
`lake env lean GaussianChain/ChordRatio.lean` check after its sign
normalization was repaired; both theorem APIs are preserved.
Both latest checkers for items 532 and 533 passed, and all 459 new Python
files and nine JSON files passed syntax parsing. These checks do not
certify every prose argument or the unfinished directions above.

## The exact goal and what is actually known

For every fixed `C>0`, prove the existence of a finite `M(C)` such that
every arc of length at most `C sqrt(R)` on the **origin-centered** circle
`x^2+y^2=R^2` contains at most `M(C)` lattice points, independently of
`R` and arc position. The radius need not be an integer. Equivalently,
the lattice points are Gaussian integers of one common modulus `R`.

It suffices to prove this for one positive constant, such as `C=1/2`:
subdivide longer arcs into a bounded number of pieces. The user wants
an improvement in the growth rate, ultimately a constant; improving
only the leading coefficient of a growing bound does not complete the
task. No sharp constant is needed.

**The uniform theorem is unproved.** The strongest current general prose
bound is

```text
M <= (1+epsilon) log R / loglog R
```

for fixed `C,epsilon>0` and sufficiently large `R`, uniformly in arc
position. Its proof is [inert_prime_cofactor_bound.md](inert_prime_cofactor_bound.md).
The Lean-verified asymptotic theorem retains leading constant `8+epsilon`.
Recent results are restricted theorems, reductions, or finite checks;
none improves this general dependence on `R`.

## Read these first

| File | What it supplies |
| --- | --- |
| [inert_prime_cofactor_bound.md](inert_prime_cofactor_bound.md) | Strongest current general growth bound; a natural starting point for independent review. |
| [endpoint_full_profile_quantifiers.md](endpoint_full_profile_quantifiers.md) | The actual extraction from a hypothetical unbounded circle cluster, with all scale losses. |
| [exact_nested_profile_residual_extraction.md](exact_nested_profile_residual_extraction.md) | New exact extraction: no corrections, nested prime overlap, and aggregate residual height `O(k^2 W/M)`; a superlogarithmic lower bound in this class would improve general growth. |
| [boolean_section_kummer_uniform_bound.md](boolean_section_kummer_uniform_bound.md) | New geometric theorem: at most 35 nonanchor rows in the equal-degree, constant-imaginary Boolean polynomial model, including repeated roots. No integer-profile transfer is proved. |
| [residual_height_gap_growth_quantifiers.md](residual_height_gap_growth_quantifiers.md) | Exactly which new height inequality would improve growth or prove uniformity. |
| [integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md) | An equivalent integer divisibility/lcm target, without a template-classification assumption. |
| [ordered_pair_triple_gap_denominator.md](ordered_pair_triple_gap_denominator.md) | A second sufficient target on five ordered actual points: control the forced pair-versus-triple conductor gap below exponent `1/15`. |
| [thin_edge_divisor_matching_criterion.md](thin_edge_divisor_matching_criterion.md) | Latest direct-divisor restriction: a four-witness bridge needs four disjoint source edges; incidence alone does not supply thin outputs. |
| [resonant_fibonacci_cyclotomic_reduction.md](resonant_fibonacci_cyclotomic_reduction.md) | Factorization, converse, and the history of restricted counts; Section 10 links the latest theorems. |
| [cyclotomic_affine_rank_reduction.md](cyclotomic_affine_rank_reduction.md) | Exact allocation map and the improved whole-fiber estimate `2rho+1`. |
| [cyclotomic_uniform_rank.md](cyclotomic_uniform_rank.md) | Degree-independent count for all odd layer indices, arbitrary prime powers and multiplicities. |
| [mixed_affine_cyclotomic_rank_scope.md](mixed_affine_cyclotomic_rank_scope.md) | Why the one-class rank theorem does not extend class by class; exact mixed cancellation, reciprocity, and the bounded-integer-coefficient gap. |
| [common_slope_cyclotomic_uniform_count.md](common_slope_cyclotomic_uniform_count.md) | New uniform count for arbitrary intercepts sharing a slope: integer coefficient bounds and three-prime amplification force classwise vanishing. |
| [frequency_separated_affine_uniform_count.md](frequency_separated_affine_uniform_count.md) | Extension to arbitrarily many slope values with distinct 2-adic valuations; also supplies the weighted colored-kernel estimate used by the next theorem. |
| [divisibility_chain_affine_uniform_count.md](divisibility_chain_affine_uniform_count.md) | Integral packet projection handles divisibility chains and two arbitrary slopes per 2-adic bucket, even with overlapping frequencies. General incomparable families remain open. |
| [three_class_affine_uniform_count.md](three_class_affine_uniform_count.md) | Further uniform extension: at most three original affine classes per bucket, with arbitrary incomparable slopes; exact two-term phase separation replaces prime isolation. |
| [multiprime_packet_cost_obstruction.md](multiprime_packet_cost_obstruction.md) | Conditional one-shot factor-three integral split, an unavoidable two-prime cost increase, and unbounded cost for complete base-packet refinement. |
| [cyclotomic_five_row_subendpoint_family.md](cyclotomic_five_row_subendpoint_family.md) | Actual five-point family with normalized arc length tending to zero; it refutes the suggested universal four-row model bound. |
| [gaussian_moat_composite_block_transfer.md](gaussian_moat_composite_block_transfer.md) | The audited transfer from the user's `openai/math` source, including its limitations. |
| [four_row_one_coordinate_divisor_reduction.md](four_row_one_coordinate_divisor_reduction.md) | Arithmetic tests on the genuine four-row profile; Sections 5--7 prove the sign restriction, dependence among the three gcd tests, and an exact residual-content identity. |
| [positive_hankel_full_profile_budget.md](positive_hankel_full_profile_budget.md) | A positive determinant sensing repeated cotangent residues; its guaranteed core-divisor budget has strict slack for every proper minor in every fixed dimension. |

## A sufficient target for the full circle problem

The following reduction applies to arbitrary hypothetical unbounded
clusters, not merely to Pell families. Subdivide to `C_0<sqrt(2)` and
remove the entire cluster's Gaussian gcd. Write

```text
M = remaining cluster cardinality,
W = log(R^2), lambda_0 = log(sqrt(2)/C_0) > 0,
k = m+1, b = 2^m, w = W/b.
```

For any fixed number `m` of nonanchor rows, an unbounded cluster supplies
the following genuine Gaussian integers, for each nonempty `T subset [m]`:

```text
P_i = K_i product_(T contains i) H_T = X_i+iY_i,
gcd_G(P_i,bar(P_i)) = 1,
n_T = Norm(H_T),                  log n_T = w+o(w),
G_ij = product_(T contains i,j) n_T,
Im(bar(P_i)P_j) = G_ij t_ij,      t_ij in Z minus {0},
Y_i != 0.
```

The blocks are conjugate-primitive and have pairwise disjoint
rational-prime supports. The `K_i` can overlap core primes: do not
silently add coprimality assumptions. The `P_i` are primitive anchor
numerators with `z_i/z_0=P_i/bar(P_i)`; they **need not have equal norms**.

Define logarithmic residual height

```text
H = max_i log max(1,|K_i|,|Y_i|),
    also maximized with max_(i<j) log|t_ij|.
```

For `M>=16k^2`, the audited bounds are

```text
W >= 4 lambda_0 (M-1),
epsilon_profile <= 8 k b / sqrt(M),
H <= 4 k W / sqrt(M) = (4 k b / sqrt(M)) w.
```

Thus a fixed positive bound `H>=c(m)w`, for one fixed `m>=4`, uniformly
on a fixed neighborhood of the balanced profile at sufficiently large
`w`, would prove the original uniform theorem. At `m=4`, there are
fifteen blocks, eight incident to each row and four shared by a pair.
Small-height examples already exist at lower row counts, so a universal
gap there is false.

A weaker lower bound can still be valuable. At fixed dimension,
`H>=w/log w` would give `M=O((loglog R)^2)`. To beat the current growth
rate by this comparison, a regular lower bound must exceed
`sqrt(w log w)` by an unbounded factor. These are conditional targets,
not established inequalities. Read the quantifier note before claiming
an implication for `R`.

There is now an alternative exact extraction with prime overlap along
two nested chains at each prime. It has `K_i=1` and
`H_sum<=k(k-1)W/[4(M-1)]`, with exact multirow Gaussian gcd identities.
Its pattern error remains `O(2^k/sqrt(M))`. For this class, a hypothetical
fixed-dimensional lower bound `H_sum>=f(w)` would imply
`M<=max(Q,1+k(k-1)W/[4f(W/2^m)])`. Thus `f(w)/log w -> infinity`
would improve growth, and `f(w)=c w` would give uniformity. The
[new note](exact_nested_profile_residual_extraction.md) proves the
extraction and conditional implications. No growing arithmetic lower
bound has been supplied, and disjoint-support height theorems do not
automatically cover the nested class.

## Other formulations for a fresh attack

**An equivalent integer problem.** Take distinct integers `X_i>=L>0`
such that every `(X_i X_j+L^2)/(X_i-X_j)` is integral. Adjoin an anchor,
and let `Q_0i=X_i`, `Q_ij=(X_i X_j+L^2)/(X_i-X_j)`. For every edge put

```text
g_e=gcd(Q_e,L), a_e=Q_e/g_e, b_e=L/g_e,
n_e=(a_e^2+b_e^2)/(2 if a_e,b_e are both odd, else 1),
N=lcm_(all edges) n_e, A=min_i X_i.
```

The associated primitive circle has squared radius `N` and exact
normalized arc length `2 N^(1/4) arctan(L/A)`. The full uniform theorem
is equivalent to the existence of one fixed `k` and `c>0` such that
every such system of `k` finite coordinates satisfies `N>=c(A/L)^4`.
Both `L` and the coordinates may grow. Keep all edges: the anchor-only
lcm can miss opposite Gaussian orientations. The known five-point
subendpoint family rules out `k=4` (as well as smaller choices);
an affirmative version would need at least five finite coordinates.
The linked note proves the equivalence and describes the endpoint-swap
symmetry. This target remains unproved.

**An ordered five-point sufficient criterion.** For a primitive ordered
five-tuple of squared radius `N_5`, let `N_outer` be the squared radius
after primitively normalizing its first and last points, and `N_inner`
the squared radius after independently normalizing its middle three.
Set `K=N_5/gcd(N_5,N_outer N_inner)`. A uniform estimate

```text
log K <= delta log N_5+B,       delta<1/15,
```

for every such tuple on a `(1/2)sqrt(R_5)` arc would prove the full
uniform count. The elementary estimate only reaches exponent `1/4`.
See the [exact denominator note](ordered_pair_triple_gap_denominator.md)
and [ordered sampling proof](ordered_conductor_index_uniformity_criterion.md).
No estimate below `1/15` has been proved.

## Latest work and unfinished leads at handoff

**A new geometric route bounds the polynomial model independently of
degree.** The [Kummer-cover theorem](boolean_section_kummer_uniform_bound.md)
proves `m<=35` when every nonempty Boolean Gaussian block has the same
positive polynomial degree, distinct blocks and their conjugates are
coprime, and every row has nonzero constant imaginary part. Repeated
roots within blocks are allowed. Pair differences determine a section
arrangement on a ruled surface. A smooth double Kummer cover would have
positive canonical square but a negative Bogomolov--Miyaoka--Yau deficit
if `m>=36`, giving a contradiction. The proof and local covering
construction have been independently audited. This is a new obstruction
to arbitrarily large polynomial constructions; it does not interpolate
finite integer configurations into such families. Extending the argument
to varying imaginary parts and nested overlap is a concrete next
algebraic question. A separate arithmetic height argument would still
be needed for the original goal.

**Checked direct-divisor restriction.** The
[new edge-divisor note](thin_edge_divisor_matching_criterion.md) proves
that four corrected divisors of row/pair words cannot form the
common-support Hadamard quartet if two source edges share a vertex,
unless one changing group's norm divides the product of those two
source correction norms. This includes full prime-power depth and
correction/core overlap. With subpower corrections and three growing
groups, the sources must be four vertex-disjoint edges: at least eight
original points are needed. Eight vertices permit the incidence
pattern, but division does not preserve small imaginary parts. The
proof has been independently audited; it does not establish the
missing bridge to the four-witness height gap.

The same note proves an exact restriction on thin division. If
`V=DU` in `Z[i]` and `D` is nonreal, then

```text
|U|^2 <= |Im V| |U| + |V| |Im U|.
```

If both imaginary parts are subpower and `|V|<=exp(alpha w+o(w))`,
then `|U|<=exp(alpha w/2+o(w))`. The scale is sharp even for
conjugate-primitive factors of coprime norms:
`(a+2+i)(a-i)=(a+1)^2-2i` for positive even `a`.
An additive construction still needs compatible outputs, rather than
an individual divisor estimate.

**Checked Hankel refinement.** Section 4 of the
[Hankel note](positive_hankel_full_profile_budget.md) now computes the
leading residues at singleton-core primes. The norm equations alone
do not force an extra singleton divisor; exact fixtures have nonzero
second-minor residues. Those fixtures are not endpoint configurations.
Simultaneous constraints across primes together with small imaginary
parts remain a possible input. The extended checker was rerun for this
handoff and passes 44 combinatorial budgets, four positive minors,
all asserted core divisibilities and 12 leading-residue identities.

**The nested extraction is now proved and independently audited.**
The [new note](exact_nested_profile_residual_extraction.md) selects one
tuple simultaneously satisfying balanced pattern weights and a small
sum of all pair-distance excesses. Keeping every threshold layer gives
exact row products, exact multirow gcds, and the aggregate bound stated
above. This strengthens the earlier many-pair `1/M` estimate in
[the near-real divisor reduction](kurlberg_wigman_shrinking_arc_audit.md)
to one balanced tuple with all its pair residuals controlled.
The remaining task is an arithmetic height lower bound for this exact
nested class. Its overlap-aware Hankel divisor is valid and can be
stronger at individual primes, but supplies no universal exponent gain.

**Most recent exploratory checks supplied no further height bound.**
The quotient-row lattice index remains controlled by the pair residuals
in the nested class, but a short lattice basis can require coefficients
of exponential height in `w`; small residuals do not remove that cost.
Likewise, fixing the residual integers does not bound the coefficients
of the central weighted-moment curve arising from the cotangent
eigenvectors. A finiteness theorem for a fixed curve is therefore not
yet a uniform height estimate for these varying curves. See
[small_imaginary_row_lattice_transference.md](small_imaginary_row_lattice_transference.md),
[residue_content_descent.md](residue_content_descent.md), and
[integer_cotangent_eigenvector_heights.md](integer_cotangent_eigenvector_heights.md).
These observations are recorded to prevent repeating an unsupported
inference; they are not additional progress on the general growth rate.

**Exploratory ordered-Ptolemy idea.** In the exact positive identity
`T_2 X_2=T_1 X_1+T_3 X_3`, a matching conductor product dominating
the sum of the other two forces a large residual matching product.
Summing these logarithmic dominance costs over ordered quadruples
could supplement inert-prime divisibility. The obstacle is a lack
of guaranteed dominance: a formal three-matching cut fixture with
products `5,13,17` satisfies the necessary pairwise endpoint distances
but `17<5+13`. It is not an endpoint circle. This idea has no general
growth consequence at present; actual phase compatibility would need
to exclude or charge near-ties. See
[ordered_residue_growth.md](ordered_residue_growth.md),
[ordered_inert_sieve_route.md](ordered_inert_sieve_route.md), and
[composite_inert_crt_collision_budget.md](composite_inert_crt_collision_budget.md)
before combining costs already counted by the inert-prime argument.

## Restricted branch: exact resonant Fibonacci factorization

Let `F_j` be Fibonacci numbers and, for positive odd `d`, set

```text
G_d = F_((d+1)/2) + i F_((d-1)/2),   Norm(G_d)=F_d.
```

The exact oriented identity is

```text
gcd_G(G_d,G_e) is associated to G_gcd(d,e).
```

This matters for unequal affine rates: apparently different factors
can share growing Gaussian divisors in the same orientation. Raw
factor degrees do not give the primitive radius.

For one proportional class `d_s(n)=m_s d_0(n)`, with fixed odd `m_s`,
write `phi=(1+sqrt(5))/2`, `alpha=sqrt(phi)`, `beta=i/alpha`. There are
integral Gaussian layers

```text
C_1 = G_d0,
C_e = Phi_e(alpha^d0,beta^d0), e>1 odd,
G_(m d0) = product_(e|m) C_e.
```

Distinct fixed layers and opposite orientations have bounded gcds,
including valuation depth. Nonproportional affine classes have bounded
mutual gcds. After removing the complete growing common product, the
true primitive radius satisfies

```text
log R_n = (log phi)/2 * sum_c d_c(n) sum_e totient(e) W_(c,e) + O(1),
```

where `W_(c,e)` is the range of row exponents at that layer. Parameters
and the implied constant are fixed while `n` grows.

For one class, even base-layer parity and equal-norm Gaussian prefactors
give monic reciprocal integer polynomials

```text
A_i(z) = (1-z)^a_i1 (1+z)^(W_1-a_i1)
         product_(e>1 odd) Phi_e(z)^a_ie Phi_e(-z)^(W_e-a_ie),
0 <= a_ie <= W_e,
W_1 and every a_i1 even,
L = sum_e totient(e) W_e,
z = i^d0 phi^(-d0).
```

All have constant coefficient one, even degree `L`, unit-circle roots,
and the same product `A_i(z)A_i(-z)`. The endpoint condition forces

```text
h_ij = ord_0(A_i-A_j) >= ceil(L/4).
```

Every `h_ij` is odd. The coefficients at positions **strictly below**
`ceil(L/4)` agree; this is not an assertion of that many equal
nonconstant coefficients. Their common odd jets can be nonzero: a
degree-40 four-row fixture has common coefficient `-1` at position one,
so the arc center drifts much earlier than pairwise separation.

The converse is proved: arbitrary finite layer orientations with the
base parity lift to actual fixed Gaussian templates, using divisor-poset
inversion, common exponent shifts and fixed equal-norm prefactors.
There is no extra incidence obstruction. However, the resulting endpoint
constant can depend on the template; arbitrarily many polynomial rows
with uncontrolled constants would not alone disprove the original
fixed-`C` theorem.

### The useful linearization and proved counts

For a pair difference `b_e=a_ie-a_je`, let `c_e(k)` be the Ramanujan sum:

```text
S_k = sum_e b_e c_e(k) = sum_(d|k) d B_d,
B_d = sum_(e:d|e) mu(e/d) b_e.
```

Contact is exactly the first odd `d` with `B_d!=0`. In particular,
low contact constraints fix `B_d` below the minimum pair contact `h`.
Projection of the row family to coordinates `e>=h` is injective by
descending inversion.

Proved in the new note, with arbitrary multiplicities:

- All active nonbase indices odd prime powers: **at most four rows**,
  sharply. Low moments force each low prime-power chain to vary as one
  coordinate; every higher coordinate costs at least `2h`, against a
  total budget `L<=4h`.
- At most `r` distinct prime divisors per active index: the newer rank
  argument gives at most `2 floor(4/c_r)+1` rows, where
  `c_r=product_(j<=r)(1-1/p_j)` for the first `r` odd primes. This is
  **15 for r=2**, **17 for r=3**. The older projection counts were
  128 and 256 respectively.
- **Unrestricted odd layer indices and arbitrary multiplicities:
  at most 2251 rows.** The new prime-power compression proves this
  for the entire one-class model, independent of polynomial degree.
  It supersedes the intermediate degree bounds `exp(O(sqrt(log L)))`
  and `O(sqrt(log L))`. For actual growing fixed templates in that
  class, the alternative small-arc estimate gives
  `min(2251,1126 max(1,ceil(C/sqrt(2))))`. These counts are uniform
  over the template parameters, but the asymptotic entry threshold
  may depend on the fixed template. They do not cover arbitrary circles.
- Every profile through degree 40 was exhaustively checked: 20,886
  profiles, 1,506,841 orientations, maximum fiber four. This is finite
  evidence, not a theorem in unbounded degree.

The one-class fiber now has a uniform upper bound. A universal
bound of four is **disproved** by the five-row construction with
`L=32812`, `h=8211`, and `L<4h`. Its explicit Gaussian rows lie on
arcs of length `o(sqrt(R))`, so `M(C)>=5` for every fixed `C>0`.
It supplies only five points, not unbounded cardinality. The proved
uniform model bound was extended in items 523--528 below, including
some overlapping slope families. General incomparable slopes and the
transfer to arbitrary circle configurations remain open.

## Other results and barriers worth retaining

- [quadratic_unit_template_bound.md](quadratic_unit_template_bound.md)
  proves at most **eight** rows for its general fixed-shift quadratic-unit
  templates, including both unit norms under its stated nonreal-phase
  hypothesis. This is stronger than the rough 4096 count in the later
  specialized Pell note; do not mistake that rough count for the best
  fixed-shift result. There is no uniform entry threshold for varying
  templates and no classification of arbitrary circles by these templates.
- [primitive_eight_point_translated_pell_family.md](primitive_eight_point_translated_pell_family.md)
  supplies eight actual primitive points with global gcd one and a fixed
  endpoint constant, attaining that restricted bound. It is not an
  unbounded-count family. The circle center remains the origin.
- [four_signed_witness_matching_norm_gap.md](four_signed_witness_matching_norm_gap.md)
  gives a positive height gap for four compatible signed products on one
  common support. The missing bridge is producing those witnesses from
  genuine rows with different supports. Its absolute small-imaginary
  hypothesis cannot simply be replaced by containment in a translated
  thin rectangle. The
  [parity extension](four_signed_witness_parity_gap.md) still admits
  unbounded abstract Walsh families.
- [full_profile_monomial_height_classification.md](full_profile_monomial_height_classification.md)
  shows that cheap multiplicative words at original row height yield
  only rows and pair quotients; other integer exponent patterns cost
  at least a factor `3/2` in logarithmic norm. Additive constructions and
  controlled Gaussian division remain separate possibilities.
- [endpoint_triangle_frame_holonomy.md](endpoint_triangle_frame_holonomy.md)
  gives actual small-determinant triangle frames. Their fixed-word loop
  traces depend only on the residuals; a split-prime fixture disproves
  automatic preservation of a second isotropic line.
- Section 5 of [the one-coordinate note](four_row_one_coordinate_divisor_reduction.md)
  proves `gcd(A_i,X_i t_ij+Y_i x_ij)=h_i`, where
  `h_i=Norm(K_i)n_{ {i} }`, `e_i=n_(I\{i})`, `A_i=h_i e_i`.
  It recovers the co-singleton norm and permits at most one primitive
  nearest-square sign when `0<|t_ij|<n_{ {i} }`. It does not force a
  height gap or prove that the surviving sign realizes all rows.
- [four_row_real_polynomial_profile_certificate.md](four_row_real_polynomial_profile_certificate.md)
  certifies a full fifteen-block real algebraic polynomial profile.
  Its coefficients are not proved rational. Thus real positivity alone
  cannot exclude this profile, but it is not an integer counterexample.
  The [rational quartic construction](four_row_rational_quartic_construction.md)
  uses only eight supports and gives a fixed five-point family, a
  different statement.
- [integer_cotangent_lcm_height_target.md](integer_cotangent_lcm_height_target.md)
  gives an alternative exact integer formulation. Retain **all edges**
  and their reduced denominator norms in its lcm formula. A naive lcm
  of anchor norms can miss opposite Gaussian orientations.

## The user's external source

The user explicitly requested ideas from <https://github.com/openai/math>.
The inspected revision is `adc7f1241b42e322a6451854ab7e4b4c146bf78a`.
The main transfer came from
[Bounded-Step Walks on Gaussian Primes](https://github.com/openai/math/tree/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Bounded-Step-Walks-on-Gaussian-Primes-September-26-2026).

Its random-conjugation thin-rectangle lemma was extended to disjoint
conjugate-primitive **composite** blocks. This gives exponentially sparse
exceptional sign sets, but actual circle orientations are not independent
random signs. The common signed factors of two actual chords meet a
critical quarter-weight threshold, without the positive slack required
by the argument. On a short circle arc, a nonzero displacement determines
its starting point, preventing the source's cheap shared-displacement
entropy hypothesis. The literal zero-avoidance sieve is constant on a
fixed norm fiber. The transfer note proves these limitations explicitly.

The pinned
[Jacobsthal manuscript](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/A-quadratic-bound-for-Jacobsthals-function-September-25-2026/paper.pdf)
was also inspected. Its independent product-box hypothesis has not been
obtained here; its curve estimate does not improve our circle bound.
Do not import either manuscript's whole conclusion as an established
input to this problem. New source claims need an exact hypothesis and
dependency audit. Local temporary source copies may exist under
`/tmp/openai-math-reference/`, but the permanent notes and pinned URLs
should be sufficient.

## Verification and suggested next work

From the local repository root, the latest dependency-free exact checks are:

```bash
python3 docs/check_resonant_fibonacci_cyclotomic_reduction.py
python3 docs/check_resonant_cyclotomic_fibres.py
python3 docs/check_four_row_integer_five_star_lift.py
python3 docs/check_cyclotomic_weighted_nullity.py --full
python3 docs/check_cyclotomic_prime_power_compression.py
python3 docs/check_cyclotomic_five_row_subendpoint.py --k 19 --d 1 3
python3 docs/check_cyclotomic_odd_subset_search.py
python3 docs/check_mixed_affine_prime_amplification.py
python3 docs/check_three_class_affine_promotion.py
python3 docs/check_multiprime_packet_cost.py
python3 docs/check_positive_hankel_full_profile_budget.py
python3 docs/check_thin_edge_divisor_matching_criterion.py
python3 docs/check_exact_nested_profile_residual_extraction.py
python3 docs/check_boolean_section_kummer_uniform_bound.py
```

Recorded runs of the first eight pass. The preceding handoff reran the
third, fifth, and eighth after their latest edits; the other five are
earlier recorded checks. The three newly appended checkers accompany
items 528--530; their latest run results are recorded in the status
journal. Do not confuse a finite check with the general proof.
For this handoff the extended Hankel checker and the new direct-divisor
checker were rerun. The latter covers all 28,035 four-edge sets in
`K_5,...,K_8`, 1,000 sharp factorization fixtures and 144 sign/depth
cases. This verifies finite examples and incidence, not the missing
analytic thinness of a four-witness bridge.
The nested-extraction checker was also rerun after its final extension:
5,680 local allocation profiles, 73,947 multirow Gaussian gcd checks,
5,681 original-circle gcd checks, 29,708 exact pair-numerator identities,
five reanchoring comparisons with 50 nonzero pair residues, 5,680
first-moment identities, and four overlap-aware Hankel divisibilities.
The tuple-selection and conditional growth proofs were audited separately;
these finite checks do not supply their missing arithmetic lower bound.
The Kummer-cover checker was rerun after the repeated-root and genus
extensions. It passes 385 squarefree formula cases, 2,018 inertia
subset checks, 780 local repeated-root cases, 260 mixed profiles,
and 84 genus substitutions. These verify exact arithmetic and local
label bookkeeping; the covering and general-type arguments are prose
proofs, independently audited, rather than consequences of the tests.
The first checks Gaussian gcds, layer identities, normalization,
ramified factors and converse lifts. The second is the complete degree-40
search, including 70,528 within-fiber distance checks and high-coordinate
injectivity. The third includes 144 actual-coordinate gcd recoveries and
144 sign exclusions with correction overlaps, plus the new gcd-collapse
and valuation-depth checks (720 arbitrary-candidate checks and 144
depth-safe roots). It now also passes 72 residual-content identities,
including 44 nontrivial pairs, four common-content fixtures and four
failures of virtual norm preservation. The full nullity search passes every checked odd
`h<=31`, including 1,201,789 supports at `h=31`, and also verifies
the larger rank-four counterexamples. The compression checker passes
34,188 exact rational-rank comparisons, plus 4,096 arbitrary support/row-set
pairs and 10,496 compression comparisons for selected rows. The five-point and odd-subset checkers verify actual primitive circles
and a complete search within the eight odd subset masks of four fixed
rates. The mixed-affine checker verifies the Ramanujan congruences, quadratic
conjugation, threshold arithmetic, and weighted colored root groups.
It now also passes 639 packet checks and 36 exact polynomial identities,
including successive projections and a counterexample to choosing
prefix slices separately at each coordinate. These are exact finite checks,
not substitutes for the general proofs.

The new results are prose proofs independently reviewed by multiple
agents; **they are not Lean formalizations**. Items 532 and 533 did not
add Lean proofs. The full merged `lake build` passed on 2026-10-07,
as recorded above; that success checks only the Lean
modules imported by that build, not these prose results. Historical checker
counts in the status note are not all fresh reruns for this snapshot.

Suggested order:

1. Audit the exact nested extraction's quantifiers or one of the alternative
   unrestricted formulations above. Audit the cyclotomic count proofs
   if using that branch. Report any gap before building on a claimed theorem.
2. Prioritize an arithmetic residual-height lower bound for the exact
   nested class. Even a fixed-dimensional superlogarithmic lower bound
   would improve growth through item 532; retain its prime overlaps,
   fixed-width profile neighborhood, and exact gcd hypotheses.
   Alternatively, in the older prime-disjoint formulation, one concrete
   route is to exclude every admissible primitive first-row Gaussian divisor
   with imaginary part below `exp(epsilon w)`, using the exact gcd/sign
   test and the reconstruction congruences, or force a correction
   or determinant residual to that height. Another is to construct the
   matching-norm witness pattern, or two sufficiently small independent
   relations, from actual Boolean rows with a proved coefficient budget.
   The direct-edge divisor criterion now rules out a common-support
   quartet from overlapping source edges; an eight-point matching
   escapes that incidence obstruction but still needs a thinness proof.
   A different sufficient criterion for arbitrary circles is welcome.
3. The model count now covers, within each 2-adic slope bucket, either
   a divisibility chain, at most two arbitrary slopes, or at most three
   original affine classes. Read both packet theorems before revisiting
   overlapping frequencies. The triple `6,10,14` is now covered when
   each slope has one class; three slope groups with unbounded numbers
   of intercept classes remain open. Seek an operation that preserves each original class's
   low coefficients and does not increase weighted integer width cost.
   Bounding the kernel of the full composition avoids a loss at each
   step. Arbitrary-selected-row compression alone does not preserve
   omitted coefficients. Rational kernel dimension can grow, so retain
   the coefficient bounds on actual integer differences. A conditional
   multiprime split costs at most three times the original budget once;
   iterating that loss is not justified. Complete base refinement has
   an exact unbounded cost ratio, so a stopping or compensation argument
   must use more than its formal factorization.
   Alternatively seek an unbounded counterfamily
   with the endpoint constant controlled. The present five-point family
   cannot be extended to six merely by selecting more of the eight odd
   masks at the checked rates. Keep the arbitrary-circle transfer explicit.
4. Record each new assertion with its exact assumptions, uniformity
   parameters, proof, and a counterexample search where useful. Keep
   conjectures, finite evidence, restricted results, and the full goal
   explicitly separate.

Do not assume uniform exceptional sets from fixed-`S` Subspace Theorems,
discard Gaussian gcd or denominator costs, replace a drifting common
phase by zero, or treat a template-dependent threshold as uniform in
the radius. Global Gaussian gcd one does not imply that each point is
conjugate-primitive, and a small correction overlap does not authorize
discarding an entire core prime at every valuation depth. The earlier
[quantitative audit](codex_quantitative_audit_result.md)
documents the first of these pitfalls. The user encourages autonomous
work and fresh approaches; the standard of completion is a proof of
the stated radius-independent count, not another constant improvement.

## Latest arithmetic and rank results

**The three private-norm gcd tests are largely dependent.** Fix a row
`i`, put `A=A_i`, and for any candidate integers `X,Y` define

```text
L_j=X t_ij+Y x_ij,
d_j=gcd(A,L_j),       d=gcd(A,L_j for all j!=i).
```

The triangle identities underlying equation (11) of the one-coordinate
note give `A | t_ik x_ij-t_ij x_ik`. Hence
`A | t_ik L_j-t_ij L_k`, and therefore `d_j | t_ij d_k`
for every `k`. Taking the gcd over `k` proves

```text
d_j/d divides |t_ij|.
```

Thus the outputs agree up to residual-sized factors, exactly agreeing
in valuation at primes not dividing `t_ij`. This argument applies to
arbitrary candidate coordinates and permits correction/core overlap.
It prevents charging the three gcd tests as independent large-modulus
conditions. It does not invalidate the sign-exclusion lemma or prove a
height gap.

**Actual residuals satisfy an exact imaginary-content identity.**
For conjugate-primitive rows with `Im(bar(P_i)P_j)=G_ij t_ij`,

```text
gcd(|Y_i|,|t_ij|)=gcd(|Y_i|,|Y_j|)=gcd(|Y_j|,|t_ij|).
```

It follows from coordinate primitivity and `G_ij|Norm(P_i),Norm(P_j)`.
If `g=gcd_i |Y_i|`, then `g` divides every residual and is coprime to
each `G_ij`. Dividing every imaginary coordinate and residual by `g`
preserves determinant congruences, but generally destroys the Gaussian
block norm factorization. This is a local filter, not a height gain or
a valid normalization of the full profile. Section 7 gives an exact
counterexample to the latter interpretation.

**Rank gives a count, but the proposed constant three is false.**
For a finite set `E` of active odd layer indices and odd `h`, define

```text
M_(d,e)=mu(e/d) if d divides e, and 0 otherwise,
             d odd, 1<=d<h, e in E,
D(E)=2*1_(1 in E)+sum_(e in E, e>1) totient(e).
```

Write `rho=dim_Q ker(M)`. The exact allocation bridge is now proved:
endpoint phase alignment forces all equal-norm prefactor ratios to
equal `(bar(F)/F)^((A_i1-A_01)/2)`, `F=1+2i`. Thus prime allocations
are exact affine functions of the layer orientations, including all
overlaps and the true common gcd. Actual growing templates on a
`C<=sqrt(2)` arc have at most `rho+1` rows. Explicitly assume positive
effective degree: stationary templates with arbitrary fixed prefactors
do not satisfy the phase-limit argument.

There is also a wholly algebraic bound. The canceled reciprocal degree
gives weighted pair distance at least `2h`; the affine box embedding
then has pairwise nonpositive inner products. A Gram-matrix argument
gives at most `2rho+1` polynomials in any fiber, or `rho+1` if `L<4h`.

The [rank-three conjecture is false](cyclotomic_nullity_four_counterexample.md):
take `h=4199` and the divisor union of `(4199,4275,4335,4845)`.
Its charged cost is `16792<4h=16796`, and its nullity is exactly four.
There is also a squarefree counterexample. The saved exhaustive check
through `h=31` still has maximum three; its finite range is essential.

The positive replacement for **all odd supports** is `rho<=1125`.
Prime-power compression cannot increase charged cost or decrease nullity
and reaches a divisor downset, whose nullity counts high indices.
For a prime `p`, write support slices `E_j` at exponent `j` and process
levels `a` downward. Put `B=intersection_(j<=a) E_j`; replace `E_a`
by `B` and each lower `E_j` by `E_j union E_a`. On kernel vectors,
subtract `f_a-f_(a+1)` from all slices through `a`. The lost kernel
embeds independently in the new kernel, proving nullity monotonicity.
The geometric totient sum pays for all added lower columns, including
the special base charge. This works for arbitrary prime powers.

The resulting downset's root-of-unity union
has size at most `4h`. Joining high orders with gcd at least `2h/15`
gives graph degree at most 224; six independent vertices would have
union size greater than `4h` by inclusion-exclusion. Hence there are
at most 1125 high orders. This gives the whole-fiber bound 2251 without
any index or multiplicity restriction. The proof and its exact rational
rank checker have been independently audited. It still does not classify
arbitrary circle configurations.

The same prime-power compression now works for **any selected set of
odd rows** `R` of the Möbius matrix. It produces a divisor downset `F`
without increasing either ordinary totient cost or the charged base
cost, with `rho_R(E)<=rho_R(F)=|F\R|`. No interval condition on `R`
is needed. It does not preserve values at omitted rows; that limitation
is crucial when trying to transport coupled equations.

**Mixed classes retain reciprocity but require bounded integer coefficients.**
For `d_c(n)=A_c n+B_c`, put `xi=phi^(-n)` and
`gamma_c=i^B_c (-1)^(A_c n/2) phi^(-B_c)` on a parity subsequence.
The coefficient at frequency `t` couples the classes with
`A_c|t` and odd `t/A_c`:

```text
sum_c S_(c,t/A_c) gamma_c^(t/A_c)/(t/A_c)=0.
```

The three classes `(A_c,B_c)=(2,1),(2,3),(2,5)` with base widths
`(2,6,2)` give actual leading cancellation and a two-point subendpoint
family: effective degree 20, first contact 6. Classwise cancellation
is therefore false. The automorphism fixing `i` and taking
`phi` to `-1/phi` sends `gamma_c` to its inverse. Together with total
base parity, this proves the same pair-degree/contact inequality and
the `2rho+1` count for the full coupled rational kernel.

However, take `N=8K` classes with `A_c=2`, `B_c=2c+1` and width one.
The `K` relevant low frequencies impose at most `2K` rational equations,
leaving nullity at least `6K` under the formal endpoint budget. Yet
no nonzero vector in `{-1,0,1}^N` cancels the first frequency:
the leading term of `sum b_c(-phi^(-2))^c` dominates its entire tail.
Thus unbounded rational nullity neither supplies a large actual fiber
nor permits a uniform-rank extension of the one-class proof. The exact
hypotheses and argument are in
[mixed_affine_cyclotomic_rank_scope.md](mixed_affine_cyclotomic_rank_scope.md).

**The bounded-integer gap is now closed for common slopes, and more.**
For common slope A, all coupled kth moments take the form `P_k(q^k)`,
where `q=-phi^(-2)` and the integer Laurent coefficient mass is at most
the normalized width degree D. Conjugation forces `P_k=0` once
`phi^(2k)>D`. For `D>=2^32`, three primes recover every remaining
small class moment via `c_e(dp)=c_e(d) mod p`; their product exceeds D.
Colored ordinary-cost compression then gives affine dimension at most
1125 and at most 2251 actual rows. For smaller D, the coordinate
dimension gives `2D+1`, so a uniform common-slope bound is `2^33+1`.

Distinct slope values with distinct 2-adic valuations have disjoint
odd-harmonic frequency sets. Applying the same argument to groups with
`L/A>=2^32` and charging the remaining coordinates against L gives
at most `2(1125+2^32)+1=8589936843` rows. These deliberately crude
constants are uniform over all parameters of the stated model classes.
The asymptotic entry threshold may depend on the fixed template.
This is now a special case of the packet theorem below. General
incomparable slopes and arbitrary circles remain unresolved; the
general growth bound in R is unchanged.

**Integral packets extend the count to overlapping slope families.**
The initial packet proof is
[divisibility_chain_affine_uniform_count.md](divisibility_chain_affine_uniform_count.md).
It retains the bound `2(2^32+1125)+1=8589936843` while allowing each
fixed-`v_2(A)` bucket to be a divisibility chain of arbitrary length,
or to contain at most two arbitrary distinct slopes. Different buckets
may satisfy different alternatives. Layer indices and multiplicities
are unrestricted. All data are fixed as the template parameter grows;
the entry threshold can still depend on that data.

The mechanism and its important audit points are:

1. Put `L=sum_(c,e) A_c totient(e) W_ce` and let `tau` be the minimum
   pair contact frequency. Original endpoint geometry gives `L<=4tau`.
   Use the original class labels throughout, even when slopes coincide.
2. For a current slope `A<=L/2^32`, find an odd prime `p` such that
   `v_p(A)<v_p(B)` for every other slope `B` in its 2-adic bucket.
   At p-free odd harmonics only the A-group contributes. Quadratic
   conjugation separates sufficiently large harmonics; four Bertrand
   primes, discarding p and retaining three, recover all smaller ones.
   This proves classwise vanishing at the p-free low harmonics.
3. Within one class write `e=p^j s`, `p` not dividing `s`, and
   `u_s=b_(0,s)-b_(1,s)`. Then `sum_s u_s c_s(k)=0` at every low odd
   harmonic, including those divisible by p. Equalize the **entire**
   level-0 and level-1 slices to the cheaper of the two weighted slices.
   Coordinatewise choices can destroy cancellations and are invalid.
4. The exact identities `Psi_s(z)Psi_ps(z)=Psi_s(z^p)` and
   `Psi_(p^j s)(z)=Psi_(p^(j-1)s)(z^p)` for `j>=2`, with
   `Psi_1=1-z`, promote `(A,B_c)` to `(pA,pB_c)` by integral coordinate
   copies. Weighted width cost does not increase. Each original class's
   low physical phase coefficients are preserved. The coupled phase
   coefficient is `S_(c,k) gamma_c^k/k`; after clearing the common
   physical frequency the coefficient is `A_c S_(c,k) gamma_c^k`.
   Omitting these factors can produce false cancellation examples.
5. In the stated slope families, eligible promotions terminate with
   every active slope greater than `L/2^32`, leaving fewer than `2^32`
   active coordinates. Intermediate rows need not retain base parity
   or an actual Gaussian-circle realization; the final count uses the
   geometry of the **original** rows.
6. If `T` is the complete projection and `V` the original difference
   span, `ker(T|V)` lies in the direct sum of the original classwise
   truncated Möbius kernels. The weighted colored-root argument bounds
   that entire direct sum by 1125. Apply rank-nullity once, then the
   original reciprocal/obtuse count. Do not add a 1125 loss per promotion.

Section 6 states the more general conditional promotion criterion.
Every promoted group must satisfy both prime isolation and
`A<=L/2^32`; an abstract sequence that requires raising a group above
that threshold is not a proof. The triple `6,10,14` has no eligible
isolation prime when all three groups remain active. The two-term
argument below bypasses this obstacle when there is only one original
class at each slope.

Even resolving all affine slope sets would still leave the transfer
from arbitrary circle clusters to these templates, or the alternative
full-profile height inequality, to be proved.

**At most three original classes per bucket now suffice.** At physical
frequency `t`, the phase of a contributing class is
`(-1)^(tn/2) i^(B_c t/A_c) phi^(-B_c t/A_c)`. Two nonproportional
classes have distinct odd exponents, so their phase ratio is
`+/- phi^(nonzero even integer)`, which is irrational. An equation
with at most two integer coefficients therefore separates classwise,
with no coefficient-size threshold.

Choose a numerically smallest current slope A and a larger B in a
bucket with at most three original colors. Some odd p has `v_p(A)<v_p(B)`.
At p-free A-harmonics the B-color is absent, leaving at most two
contributors. The ordinary cost-preserving packet is therefore valid.
Every promoted slope divides the original least common multiple, and
each affected color's slope increases, so finite repetition leaves
one slope value. These steps require no `L/A>=2^32` condition. Apply
the old arithmetic cutoff only afterward, retaining formal original
colors and charging the composed kernel once. This proves the same
`8589936843` count for these new buckets, and they may be mixed with
the old bucket types. Intermediate rows still need not preserve parity
or actual circle realization. No arbitrary-circle classification follows.

**Multiprime splitting has a precise cost limitation.** Assuming an
individual class's low moments vanish at harmonics coprime to a product
P of odd primes, an integral triangular packet split preserves all low
coefficients with cost at most three times its original width budget.
This is one simultaneous linear map on the entire row family. The
arithmetic hypothesis has not been proved uniformly for general coupled
classes. Cost preservation is false even locally: the parity-valid
exponent difference `(2,0,0,-2)` on indices `(1,5,7,35)` costs 50,
while every indicated two-prime packet gauge costs at least 66. The
example is not endpoint-compatible on its own. More broadly, the exact
base-packet expansion of `b_P=1, b_1=-mu(P)` costs `sigma(P)-1`,
against original cost `phi(P)+1`; their ratio is unbounded. This rules
out unrestricted complete refinement with a universal cost factor,
not every possible stopped projection scheme.

**A full-profile positive determinant remains below the needed budget.**
Clear the actual Y and t denominators with subpower L, put
`X_i=L x_i/y_i`, `C_i=X_i^2+L^2`, and form
`H_r[a,b]=sum_i X_i^(a+b)/C_i`. The positive integer
`D_r=(product_i C_i)det(H_r)` is a sum over r-subsets S of squared
Vandermondes times complementary C-factors. A block T forces exponent
`f_T(S)=|T|-|T intersect S|+|T intersect S|(|T intersect S|-1)`.
Its guaranteed divisor has total exponent `sum_T min_S f_T(S)`,
whereas every positive summand has size exponent
`sum_T f_T(S)=(m-r)2^(m-1)+r(r-1)2^(m-2)`. The former is strictly
smaller for every `r<m`, by the singleton blocks, and equal for `r=m`,
the usual squared Vandermonde. This holds for every fixed m and actual
nonconstant Y coordinates, allowing correction/core overlaps. Additional
divisibility or a different use of the joint norm equations is needed;
positivity and these prescribed valuations alone give no height gap.
