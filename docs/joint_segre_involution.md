# A rational joint involution of six points

The quadratic deck transformation on the Igusa quartic induces a genuine
rational involution of six labelled points on the projective line. It is
generically different from every permutation of the points. Its formulas,
exceptional divisors, and arithmetic effects can all be made explicit.

For an input with the full sixty-four-block endpoint profile, the output
matching height is between `42w-o(w)` and `48w+o(w)`, compared with the
input height `18w+o(w)`. Consequently every endpoint-scale circle realization
of the output moduli has `log R_out>=28w-o(w)`. For the actual pullback
fixing three input directions, sixteen old cuts give a stronger result:
the arc-preserving ordering forces `R_out>=R_in^(49/48-o(1))`, even after
dropping the original anchor. Every labeling of this fixed-three pullback
has normalized output arc length at least `exp(w/3-o(w))`. Thus this
specified candidate cannot provide an endpoint descent. A subsequent
arbitrary fractional-linear change, other joint maps, and the uniform
circle-point bound remain unresolved.

The main exact certificate is
[check_joint_segre_involution.py](check_joint_segre_involution.py).
The actual local common-content bound uses
[check_joint_segre_local_content.py](check_joint_segre_local_content.py).

## 1. The deck operation and its Segre formula

Use the matching coordinates and cubic from
[segre_gradient_arithmetic.md](segre_gradient_arithmetic.md):

```text
A=12·34·56, B=12·36·45, C=14·23·56,
D=16·23·45, E=16·25·34,
Phi=BCE-AD(A+B+C+D+E).
```

Here `ij` denotes `Delta_ij=det(v_i,v_j)`. Write
`(p,u,v,q,w_0)=grad Phi`. Its dual quartic is

```text
J=(pq-uv-u w_0-v w_0)^2
       -4uv w_0(u+v+w_0-p-q).
```

As a quadratic in `w_0`, its leading and constant coefficients are
`(u-v)^2` and `(pq-uv)^2`. Its rational deck involution is therefore

```text
(p,u,v,q,w_0) -> (p,u,v,q,(pq-uv)^2/((u-v)^2 w_0)).        (1)
```

Introduce five linear forms in the original Segre coordinates:

```text
F=C-B,
H=B+C,
U=B+C+2D,
V=2A+B+C,
G=2A+B+C+2D+2E.                                         (2)
```

Pulling (1) back through the inverse quartic gradient gives the simple
homogeneous quartic map

```text
T(A,B,C,D,E)=(-A UFG, -B UVG, C UVG, -D VFG, -E F^3).     (3)
```

On the Segre cubic, exact polynomial identities give

```text
grad Phi(T(M))
 =-HUVG (pF^2,uF^2,vF^2,qF^2,w_0G^2),                  (4)
(u-v)G w_0=(pq-uv)F.
```

Thus (4) is exactly the deck operation (1), followed by the inverse
polar map, whenever the displayed denominators are nonzero. This is
not merely the standard forward/inverse polar composition.

The five forms transform, on `Phi=0`, as

```text
F(T(M))=HUVG,    G(T(M))=FHUV,    H(T(M))=FUVG,
U(T(M))=FHVG,    V(T(M))=FHUG.                            (5)
```

It follows directly that

```text
T(T(M))=(FGHUV)^3 M.                                    (6)
```

The certificate checks (3)--(5) by exact polynomial reduction modulo
the cubic. Equation (6) follows by substitution of the factored forms.

## 2. Rational coordinates, collisions, and a nonpermutation example

Normalize the first three projective points to `infinity,0,1`, and write
the remaining three as `a,b,c`. The matching coordinates are

```text
((1-a)(b-c), (1-c)(a-b), c-b, b-a, b(a-1)).
```

The forms (2) become

```text
F=c(a+1-b)-a,
G=c(a+b-1)-a,
H=c(b+1-a)+a-2b,
U=c(b+1-a)-a,
V=c(a+b-1)-a(2b-1).
```

The involution fixes the first three normalized points and sends

```text
a'=aH/V,       b'=bF/V,       c'=cH/F.                    (7)
```

The collision identities are

```text
a'-1=-(a-1)F/V,           b'-1=-(b-1)U/V,
c'-1=(c-1)U/F,           a'-b'=(b-a)G/V,
c'-b'=(c-b)UG/(FV),      c'-a'=(c-a)HG/(FV).              (8)
```

In these normalized coordinates the form transformations are

```text
(F',G',H',U',V')
 =(HUG/(FV), HU/V, UG/V, HG/V, HUG/V^2).                 (9)
```

Thus the open set where all original points are distinct and
`FGHUV!=0` is preserved, and (7) is an involution on it. These statements
follow from the displayed identities, not just from numerical testing.

An exact example is

```text
(infinity,0,1,2,3,5)
 -> (infinity,0,1,6/5,-3/5,-15).                         (10)
```

The two unordered six-point sets are not projectively equivalent.
The certificate checks all `6·5·4=120` ways to send three original
points to `infinity,0,1`; none sends the original set to the target
set. Checking sets covers the remaining permutations. Since there
are finitely many permutations, this example also rules out generic
identification of the rational map with any permutation action.

## 3. Returning to the actual circle coordinates

The formulas use a rational projective normalization, not a new choice
of Euclidean circle. For the actual input binary vectors, let

```text
chi(v_i)=Delta_i2 Delta_31/(Delta_i1 Delta_32).
```

Then `chi(v_1),chi(v_2),chi(v_3)=infinity,0,1`. A normalized output
coordinate `(N:D)` lifts exactly to

```text
v(N:D)=N Delta_32 v_1-D Delta_31 v_2.                    (11)
```

Taking determinants verifies `chi(v(N:D))=N/D`. Use (11) with the
three numerator/denominator pairs in (7), then clear rational contents
and take the actual primitive Gaussian numerators. The resulting six
rational directions admit a Gaussian circle realization by the least
common-denominator construction.

This exact lift does not assert that the output remains in the original
short arc, that its least radius decreases, or that all other points of
a larger cluster can be moved compatibly. Those are extra metric and
arithmetic requirements. The detailed lift and geometric interpretation
are in [segre_igusa_symmetry_scope.md, Section 5](segre_igusa_symmetry_scope.md).

## 4. A preserved unit constrains the output scale

The rational projective unit

```text
X=a(c-1)/(c(b-1))
 =Delta_42 Delta_63 Delta_15/(Delta_14 Delta_62 Delta_53)
 =(G+F)/(G-F)
```

satisfies the exact relation

```text
X'=-X.                                                   (12)
```

This follows either from (7)--(8) or from (5). Its full sixty-four-cut
height norm is twelve, so the actual full-profile input has
`h(X)=12w+o(w)`. If the output again had a full near-uniform profile
with block scale `w'` and negligible primitive residues and corrections,
the same height identity would force `w'=w+o(w)`.

Each of `F,G,H,U,V` is a sum or difference of two actual matching
monomials with disjoint edges:

```text
F=C-B,          H=C+B,
U=(B+D)+(C+D),
V=(A+B)+(A+C),
G=(L-B)+(L-C),       L=A+B+C+D+E.                         (13)
```

Every parenthesized term is, up to sign, one of the fifteen matching
products. Their ratios also have full-profile norm twelve. Hence none
of the five forms can vanish in an increasingly accurate full-profile
input: vanishing would give a matching ratio equal to `1` or `-1`.
This places those actual inputs in the domain used above.

Further conserved units and the scope of the height invariant are
recorded in [joint_involution_unit_arithmetic.md](joint_involution_unit_arithmetic.md).

## 5. The actual output has raw height 192w

Retain the central arithmetic notation of the Segre height note:

```text
|Delta_ij|=c_ij b_ij,
c_ij=product_(S containing i,j)n_S,
1<=b_ij,       log b_ij<=beta=2s+log T,
(1-eta)w<=log n_S<=(1+eta)w.
```

The value `2s` uses the actual central pair-gcd correction bound. If
only the abstract individual row correction bound is available, replace
it by `4s` throughout. Let `m_0,m_1` be the minimum and maximum absolute
values of all fifteen input matching products. Then

```text
exp(48(1-eta)w)<=m_0<=m_1<=exp(48(1+eta)w+3beta).          (14)
```

For every pair of the five forms, a signed sum is twice one matching
product. The ten identities are

```text
H-F=2B,        U-F=2(B+D),       V-F=2(A+B),
G-F=2(A+B+D+E),                 U-H=2D,
V-H=2A,        G+H=2L,          V+U=2(L-E),
G-U=2(A+E),    G-V=2(D+E).                               (15)
```

Consequently at most one form has modulus less than `m_0`. Every
form has modulus at most `2m_1` by (13). For each possible small form,
there is an output coordinate in (3) avoiding it: use `B'` if `F`
is small, `E'` if `G,U`, or `V` is small, and any coordinate if
`H` is small. Thus for the raw integer vector `T(M)`,

```text
m_0^4 <= max|T_i(M)| <=8m_1^4.                            (16)
```

In particular the logarithm of its maximum is `192w+o(w)` with
no assumption against cancellation in an individual form.

## 6. Actual common content: six exceptions to the generic count

At a core prime associated to `S`, the generic order of an invariant
is computed by colliding the rows in `S`. Exact polynomial expansion
computes this as the smallest sum of affine exponents on those rows.
Call the generic common order of (3) `k_S`. The exact table gives

```text
sum_S k_S=144.                                           (17)
```

These orders give actual divisibility by `n_S^k_S`: after an integral
Gaussian change of basis, each relevant term contains at least that
many factors from the inside numerators `P_i`, and the evaluated
invariant is an ordinary integer. Arbitrary core prime powers are
retained.

The generic count is not always the actual upper bound. For example,
take the local degeneration `a=1+epsilon`, with `bc=1` and
`b notin {0,1,-1}`. Only the old pair `34` collides. Then

```text
ord(A,B,C,D,E)=(1,0,0,0,1),
ord(F,H,U,V,G)=(0,0,1,0,1),
ord T(M)=(3,2,2,1,1).
```

The common order is one, whereas the generic `34` cut order is zero.
No additional old pair collision is present. This is a real extra
base-locus condition, so it cannot simply be charged to an old pair's
small primitive residue.

There is also an exact arithmetic instance. At `(a,b,c)=(6,2,3)` and
the prime `5`, the old matching vector is `(5,-8,1,-4,10)` and
`(F,H,U,V,G)=(9,-7,-15,3,15)`. Only the old pair `34` collides
modulo `5`, while the five output valuations are `(3,2,2,1,1)`.
Thus the extra common factor `5` is an actual integer phenomenon,
not just a formal degeneration. This single-prime example does not
assert the simultaneous full sixty-four-block profile.

The corrected upper table follows from two exact sets of constraints:
each individual decomposition (13), and the ten pair identities (15).
If two summands have different valuations, the valuation of their sum
is their minimum; if equal, it is at least that value. The exhaustive
certificate checks all sixty-four cuts and 640 feasible valuation words.
Its largest common order `K_S` satisfies

```text
K_S=k_S+1 for S=34,25,16 and their complements,
K_S=k_S otherwise,
sum_S K_S=150.                                           (18)
```

This finite check covers real-valued orders as well. Ceiling preserves
the constraints: unequal valuations have minimum equal to the integral
order of the matching on the other side; equal valuations have common
ceiling at most that order. At most one form can have order greater
than three. It can be replaced by twenty without decreasing the
possible common order; some output coordinate omits it and has order
at most twelve. Hence the finite choices `0,1,2,3,20` suffice.

Here is the transfer to actual correction valuations. At a rational
prime put

```text
B_p=sum_(all15 edges) v_p(b_ij)+v_p(2).
```

For an old core prime with exponent `e>B_p`, replace each actual form
valuation `x` by `ceil((x-B_p)/e)`. The shifted ceiling satisfies
the same exact table constraints. For instance, an unequal pair of
shifted ceilings comes from unequal actual valuations, whose minimum
is `q e+delta`, with `0<=delta<=B_p`; its shifted ceiling is exactly
`q`. The individual unequal-summand constraints are preserved for the
same reason. Generic lower orders are also preserved because `B_p<e`.
An output monomial has one input matching factor and three linear-form
factors, counted with multiplicity, so

```text
v_p(gcd T(M))<=e K_S+4B_p.                               (19)
```

If `e<=B_p`, the pair identities imply that at most one form has
valuation above `3e+B_p<=4B_p`. Choose an output coordinate omitting
that form. Its valuation is at most `12e+4B_p<=16B_p`.
At a prime absent from the core the same argument gives `4B_p`.
Since `sum_p B_p log p<=15beta+log2`, (17)--(19) prove

```text
144(1-eta)w <= log gcd T(M)
 <=150(1+eta)w+240beta+16log2.                            (20)
```

The generic count `144w` alone would have suggested output height
`48w`; the six explicit exceptions are why that equality is not
claimed. The exact tables are
[joint_segre_cut_valuations.json](joint_segre_cut_valuations.json) and
[joint_segre_local_content_audit.json](joint_segre_local_content_audit.json).

## 7. Output height and the circle-radius implication

Subtracting (20) from (16) gives the actual finite bounds

```text
42w-342eta w-240beta-16log2 <= h([T(M)])
 <=48w+336eta w+12beta+log8.                              (21)
```

In particular the output cannot again satisfy the same full-profile
normal form with negligible corrections and residues. Such an output
would have matching height `18w'+o(w')`; (12) would give `w'=w+o(w)`,
contradicting the lower bound `42w-o(w)` in (21).

There is also a radius consequence requiring no output-profile assumption.
For any six Gaussian lattice points on an arc of length at most
`C_0 sqrt(R_out)`, all three-chord matching products are Gaussian integers
and have one common complex phase up to sign. Their ratios are real and
lie in `Q(i)`, hence are rational. Write their vector as `gamma m`, where
`m` is a primitive integer vector. Integer Bezout gives `gamma in Z[i]`,
so `|gamma|>=1`. Every chord has modulus at most the arc length, proving

```text
h_matching <= (3/2)log R_out+3log C_0.                    (22)
```

Projective matching height is invariant under every real fractional-linear
change of coordinates, so this applies to any circle realization of the
output moduli point. Combining (21) and (22), for fixed `C_0`, gives

```text
log R_out >=28w-o(w).                                    (23)
```

The inherited input radius has `log R=32w+o(w)`. Thus (23) rules out
a radius reduction with exponent strictly below `7/8` on that inherited
scale. It does not rule out every fixed-power reduction, and it does not
assert that an endpoint-scale output realization exists. Extending this
six-point operation to all the points of a large cluster remains a
separate unresolved requirement.

The ordered arc-preserving choice proved in
[joint_involution_real_arc_control.md](joint_involution_real_arc_control.md)
gives a stronger radius floor for its actual angular realization. Its
angular diameter satisfies `Theta_out<=Theta_in<=C exp(-16w)`.
For this finite comparison fix `w=W/64`, with `W=2log R` the original
inherited norm scale, and use the corresponding core-block error `eta`.
An arbitrary comparison weight would require retaining its additional
error in the displayed angular bound.
The same common-phase argument used in (22), now bounding each chord
by `R_out Theta_out`, gives

```text
h_out<=3log R_out+3log Theta_out
     <=3log R_out-48w+3log C.
```

Together with (21), this proves

```text
log R_out>=30w-114eta w-80beta-(16/3)log2-log C
         =30w-o(w).                                     (24)
```

This requires no endpoint assumption on the output radius; it uses its
preserved old angular arc. This preliminary floor is strengthened by
keeping the actual common Gaussian factor of the matching chords.

The eight cuts `12,13,23,124,135,145,16,25` and their eight complements
each force at least half a block of common Gaussian matching content,
up to the negligible correction error. The local proofs allow both
orientation errors, arbitrary prime powers, and unbounded cancellation
depths. Complement symmetry uses the exact raw Gaussian covariants and
the standalone six-point allocation range. Consequently

```text
log|gamma|>=8w-o(w),
h_out+log|gamma|<=3log(R_out Theta_out).
```

For the arc-preserving ordering, this proves

```text
log R_out>=(98/3)w-o(w),
R_out>=R^(49/48-o(1)),
Theta_out sqrt(R_out)>=exp(w/3-o(w)).                    (25)
```

The full proof, including the finite correction bound, is
[joint_involution_standalone_inflation.md](joint_involution_standalone_inflation.md).
The old anchor need not be retained. Keeping it cannot reduce the
required radius, so (25) also strengthens the earlier anchored five-cut
result.

The last conclusion in (25) holds for every labeling of this fixed-three
pullback, even without arc preservation: the unchanged three input
directions give `Theta_out>=exp(-16w+o(w))`, and the first displayed
inequality gives `log(R_out Theta_out)>=(50/3)w-o(w)`. Adding these
inequalities proves the asserted normalized-arc growth.

This settles the arithmetic descent question for this specified operation
and actual inverse normalization. It does not exclude subsequent arbitrary
fractional-linear changes of the output directions, different joint maps,
or constructions acting on an entire cluster. The uniform bound remains
unproved.

## Independent verification

The root agent independently derived the deck ratio, checked the exact
example and nonpermutation issue, and supplied the radius inequality.
The symmetry agent supplied the homogeneous linear forms, the inverse
projective lift, the pair identities, and the exact local-content table.
The fresh-algebraic agent independently derived and checked the chart and
quartic maps, involution, collision identities, and their common-factor
normalizations. The uniformity-audit agent independently verified the
conserved unit, radius inequality, shifted-ceiling transfer, and all finite
constants in (20)--(23). The verified result is a new joint rational symmetry
with quantitative arithmetic restrictions. The standalone content proof
in (25), including all representative local cases, two-orientation
correction bounds, and exact complement covariance, received independent
cross-audits across all four agents. This is not a proof of uniformity.
