# An individual coefficient bound at every five-row boundary node

A five-row irreducible numerical-zero invariant whose anticanonical
curve passes through **any** boundary node satisfies

```text
log C_Q >= 2(1-eta)w-110sigma-2log T.                   (1)
```

Here `C_Q` is its polynomial coefficient sum, `w` is the inherited
full-cut norm-logarithm scale, `sigma` bounds the row corrections, and
`T` bounds the primitive pair imaginary residues. This is a bound for
one relation. It does not require a second independent numerical zero.

Consequently every irreducible five-row numerical zero below (1) has ten
distinct, non-nodal rational boundary intersection points. This excludes
the whole nodal stratum, beyond the special binomial pencil treated in
[five_row_del_pezzo_arithmetic.md](five_row_del_pezzo_arithmetic.md).
It does not exclude the remaining irreducible relations, and it does not
prove the endpoint uniformity statement.

## 1. Actual hypotheses and the non-scalar restriction

Use the five-row actual Gaussian system and complement-reflection
normalization from
[odd_support_relation_rigidity.md](odd_support_relation_rigidity.md):

```text
P_i=X_i+iY_i=K_i product_(U containing i) H_U,
gcd_G(P_i,bar(P_i))=1,
log|K_i|<=sigma,
log N(H_U)>=(1-eta)w,
0<|t_ij|<=T.
```

All core blocks and their conjugates have disjoint odd split-prime
support. The directions are distinct. Let `Q` be a nonzero integer
simultaneous `SL_2` invariant of degree two in each row, with `Q(P)=0`.
No smoothness assumption on its zero curve is made.

For an inside pair `S`, put its binary rows at `(0,1)`, and write the
three outside rows as `(1,z_a),(1,z_b),(1,z_c)`. Its restriction is a
translation-invariant linear polynomial. Suppose it is a nonzero single
difference,

```text
F_(S,Q)=c_S(z_a-z_b),              c_S in Z\{0}.        (2)
```

Direct specialization does not increase the polynomial coefficient sum,
so `|c_S|<=C_Q`. Let `D_ab=Im(P_a bar(P_b))`, and set

```text
v=D_ab P_c,
kappa=K_a K_b K_c.
```

Up to a harmless sign,
`v=(P_a P_b P_c)(z_a-z_b)` at `z_j=Y_j/P_j`. The actual numerical-zero
restriction congruence therefore gives

```text
H_S | kappa c_S v.                                     (3)
```

This is the one-relation congruence from the odd hierarchy, specialized
to one nonzero basis difference. Nothing is inverted modulo a prime
power, and `D_ab P_c` is nonzero at the actual configuration.

## 2. The exact coefficient modulus and its reflected partner

Since `c_S` is an ordinary integer, (3) implies

```text
M_S | c_S,
M_S=N(H_S)/N(gcd_G(H_S,kappa D_ab P_c)).                 (4)
```

The logarithmic norm loss is bounded as follows:

* `kappa` contributes at most `6sigma`;
* the outside row `P_c` contributes at most `2sigma`, since its core
  factors are units at `H_S`;
* the outside bracket contributes at most `2sigma+log T`, by the exact
  primitive-residue identity
  `D_ab=N(gcd_G(P_a,P_b))t_ab` and the outside-core coprimality.

Prime by prime, the valuation of a gcd with a product is at most the sum
of the separate gcd valuations. Thus

```text
log M_S >= log N(H_S)-10sigma-log T.                   (5)
```

Complement reflection preserves the labelled polynomial `Q` and the
same restriction coefficient `c_S`. Its row-correction bound is
`9sigma`, it preserves the pair residues up to sign, and it retains the
complementary core block with norm-logarithm loss at most `10sigma`.
Applying (4)--(5) to it gives a second ordinary modulus `M_S'` with

```text
M_S' | c_S,
log M_S' >= log N(H_(S^c))-100sigma-log T.              (6)
```

The two moduli have disjoint norm-prime support, so they multiply:

```text
M_S M_S' | c_S,
log|c_S| >= log N(H_S)+log N(H_(S^c))
                          -110sigma-2log T.            (7)
```

Equation (1) follows from (7). This argument retains arbitrary prime
powers and charges every reflected content loss. It is stronger than
merely asserting that the restriction image has rank at most one:
there is only one scalar coefficient in (2), and it is nonzero.

## 3. Why every boundary node supplies such a restriction

The five-row invariant quotient is the split degree-five del Pezzo
surface used in the earlier notes. Its ten boundary lines `D_S` are
labelled by pairs. Two meet precisely when their pairs are disjoint;
there are fifteen boundary nodes. A degree-two-each invariant defines
an anticanonical section, whose intersection degree with each boundary
line is one.

Suppose the section passes through `D_S intersect D_T`, with
`T={a,b}` disjoint from `S`, and remaining label `c`. On `D_S`, this
node is the specialization `z_a=z_b`, with `z_c` distinct. A
translation-invariant linear polynomial vanishing there is necessarily
a scalar multiple of `z_a-z_b`. Thus its restriction has form (2),
provided it is not identically zero.

If `Q` is irreducible, that restriction cannot vanish identically.
Indeed, vanishing for every configuration with the two rows of `S`
proportional would make the bracket `Delta_S` divide `Q`. Irreducibility
and the degree-two-each multidegree exclude this. Hence (1) applies at
every boundary node of an irreducible numerical-zero section.

Below (1), the section therefore contains no boundary node and no
boundary component. It intersects each boundary line in one rational
point, because its restriction is a nonzero rational section of degree
one on that rational line. These ten points are distinct: any coincidence
would be an intersection of two boundary lines, hence a boundary node.

Equivalently, writing each restriction in a fixed outside ordering as

```text
F_(S,Q)=alpha_S(z_a-z_c)+beta_S(z_b-z_c),
```

every sufficiently small irreducible numerical zero satisfies

```text
alpha_S beta_S (alpha_S+beta_S) != 0
                         for all ten inside pairs.    (8)
```

This is an actual arithmetic consequence of the full-cut system,
not just a genericity assertion about invariant polynomials.

The two incident restrictions at a node need not have the same slope
coefficient. On the five-dimensional vector space of anticanonical
sections vanishing at that node, the two slope functionals together
have rank two. The exact checker verifies this for all fifteen nodes.
Thus one may apply (7) to each coefficient separately, but may not
multiply their disjoint core moduli into one coefficient without a new
arithmetic relation between those slopes.

## 4. Inherited scale and the precise remaining case

If this five-row irreducible factor came from a degree-two-each relation
on `M` rows, support projection retains

```text
w_5=2^(M-5)w_ambient,
```

while its coefficient logarithm costs only `O_M(1)` beyond that of the
ambient relation. Its original row corrections and pair residues remain
unchanged. The nodal lower bound is therefore

```text
log C_Q >= 2^(M-4)(1-eta)w_ambient-110sigma-2log T.      (9)
```

The amplification is not discarded. Equation (9) excludes nodal
five-row factors of sufficiently short ambient relations, but does not
force the remaining factor to be nodal.

A precise next arithmetic statement for the unresolved five-row case
would be: for some fixed positive `c` and finite constants `A_0,B_0,B_1`,
every irreducible numerical zero satisfying the actual system above and
the nonvanishings (8) obeys

```text
log C_Q + A_0 sigma + B_0 log T >= c w-B_1.             (10)
```

There are further [singular-point and ramification gaps](five_row_singular_point_height.md):
a sufficiently small irreducible relation is geometrically integral,
its actual point is smooth, and all five forgetful maps are unramified
there. These restrictions still leave smooth unramified points on
smooth genus-one curves or on irreducible singular curves.

The full statement (10) is **unproved**. It is an individual-relation statement, so the
known two-relation minors do not imply it. Extending such a result to
arbitrary support sizes, or proving a descent to five rows, would still
be necessary for a global argument.

The example in
[five_row_elliptic_height_compatibility.md](five_row_elliptic_height_compatibility.md)
explains the necessary arithmetic scope of (10): a fixed smooth
irreducible section with ten distinct non-nodal boundary points has
actual rational points realizing all ten aggregate boundary heights,
with disjoint good-prime contact ideals. Its equation has bounded
coefficient height. Thus forgetting the oriented Gaussian blocks and
the small actual primitive residues would make a proposed positive
height bound false.

To get beyond (8), one must control the common lift back to the circle:
all complementary boundary contacts must be realized at the same pair
of conjugate isotropic directions in a single rational chart, while
the actual primitive imaginary residues remain small. A relation
between the two non-scalar restrictions must use this shared chart
and its conjugation; independent divisor heights and independent
boundary-contact congruences do not supply it.
