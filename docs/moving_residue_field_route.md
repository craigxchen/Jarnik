# Moving residues, the small-height field, and the missing place transfer

This investigation treats the primitive residues of a hypothetical
asymptotically uniform fixed-size cluster as moving coefficients of
negligible height. It proves a quantitative linear-independence result
that is useful for this viewpoint. It also identifies an exact obstruction
to applying existing moving-target gcd theorems with varying conductor
support. It does not prove a positive exponent of residue growth or the
uniform endpoint count.

## 1. Setup and small-height half-angle coordinates

Take `k=m+1>=3` rows, with reference row `0`, in one common Gaussian
unit class. At a split rational prime `p`, let the allocations of its
chosen Gaussian prime be `a_0(p),...,a_m(p)`. Write

```text
A_j = product_p pi_p^max(a_j-a_0,0)
                  conj(pi_p)^max(a_0-a_j,0),
z_j/z_0 = A_j/conj(A_j).
```

This definition permits different rows to use opposite conjugate primes.
No squarefree assumption or common orientation of the blocks is made.

Put `W=log(R^2)`. Assume every inherited oriented threshold pattern
on the `k` rows has weight between

```text
(1-eta) W/2^k    and    (1+eta) W/2^k,       0<=eta<1.
```

The reference-to-row cut distance then gives

```text
(1-eta)W/4 <= log|A_j| <= (1+eta)W/4.
```

On an endpoint arc with angular width at most `C exp(-W/4)`, choose
the signs of the `A_j` so their arguments are in the corresponding
short half-angle interval about the positive real axis. Writing
`A_j=X_j+iY_j`,

```text
|Y_j| <= (C/2) exp(eta W/4).                          (1)
```

Thus along a hypothetical sequence with `eta->0` and `W->infinity`,
the integer `Y_j` has height `o(W)`. For large indices it is nonzero:
otherwise the conjugate-coprime `A_j` would be a unit, contradicting
its positive logarithmic size.

The primitive residues have the same negligible-height property from
the extraction estimates. Both types of quantities can therefore be
treated as coefficients of small height, with that assertion proved
rather than assumed.

## 2. A quantitative new linear-relation obstruction

Let

```text
G = gcd_G(A_1,...,A_m),
B_j = A_j/G,
D_j = gcd_G(B_i : i!=j).
```

The gcds are Gaussian. Then

```text
log|G| >= (1-eta) W/2^k,
log|D_j| >= (1-eta) W/2^k,
gcd_G(D_j,B_j)=1.                                   (2)
```

Consequently every nonzero Gaussian-integer coefficient relation

```text
c_0 + sum_(j=1..m) c_j A_j = 0
```

satisfies the explicit bound

```text
max_(0<=j<=m) log max(1,|c_j|)
    >= (1-eta) W/2^k.                              (3)
```

This is a positive height lower bound for *linear relation
coefficients*. It is not yet a lower bound for the primitive residues:
the known reconstruction relations also contain conductor-sized
coefficients, and so do not violate (3).

### Prime-power verification of the gcd bounds

For a fixed split prime write `a_i=a_i(p)`. The two valuations of
`G` are

```text
v_pi(G)     = max(0,min_(i>=1) a_i-a_0),
v_conjpi(G) = max(0,a_0-max_(i>=1) a_i).
```

These count exactly the two threshold patterns in which the reference
row is on one side and all other rows are on the other. Each has
weight at least `(1-eta)W/2^k`. The factor `1/2` converting norm
logarithms to Gaussian moduli gives the first bound in (2).

For the private gcd `D_j`, the corresponding valuations are

```text
v_pi(D_j)
 = max(0,min_(i!=0,j) a_i-max(a_0,a_j)),

v_conjpi(D_j)
 = max(0,min(a_0,a_j)-max_(i!=0,j) a_i).
```

They count exactly the two threshold patterns in which rows `0,j`
are together and every other row is opposite. These give the second
bound in (2). This calculation retains the grouping of thresholds
inside each prime, so it does not infer a prime-power assertion from
a squarefree cut picture.

Finally, `gcd(B_1,...,B_m)=1` directly implies
`gcd(D_j,B_j)=1`.

### Proof of the coefficient bound

If `c_0!=0`, the relation makes `G` divide `c_0`, and (3) follows.
If `c_0=0`, divide the relation by `G`. For any nonzero `c_j`,
reduction modulo `D_j` gives `D_j | c_j B_j`; coprimality gives
`D_j | c_j`, proving (3) again.

Only the reference-isolated and reference-plus-one-isolated patterns
were needed here. Full uniformity is a convenient sufficient condition
for their positive weights.

As a finite arithmetic check, all 26,144 allocation vectors with
`2<=m<=5` and prime norm exponent between zero and four satisfied
the common/private gcd valuation formulas and the asserted
coprimality. The threshold proof above, rather than this finite check,
establishes the arbitrary-exponent statement.

## 3. The small-height field is legitimate, but need not consist of constants

For a precise sequence-field formulation, fix a nonprincipal ultrafilter
and take classes of sequences in `Q(i)` modulo equality on an
ultrafilter-large set. Use absolute logarithmic Weil height `h`.
Inside this ultraproduct, the classes with `h(x_n)=O(W_n)` form a
field `F`. The classes with

```text
lim_U h(x_n)/W_n = 0
```

form a subfield `K_small`. This follows from the usual height
inequalities for addition, multiplication, and inversion. The
ultrafilter is used only to avoid zero divisors from unrelated
vanishing subsequences.

Moreover `K_small` is relatively algebraically closed in `F`.
Indeed, a fixed-degree algebraic equation with coefficients of height
`o(W)` gives a root-height bound `o(W)` after division by its leading
coefficient. This follows place by place from the root bound for a
monic polynomial and then by summing the local height bounds.

Clearing the finitely many denominators in a proposed `K_small`-linear
relation multiplies their heights by only a fixed constant. Equation
(3) therefore proves

```text
1,A_1,...,A_m are linearly independent over K_small.  (4)
```

In particular, small-height linear nondegeneracy is an actual
consequence of the extracted profile. It is not a missing assumption
in this formulation.

Equation (1) gives an additional exact simplification:

```text
A_j-conj(A_j)=2iY_j belongs to K_small.
```

After dividing by the nonzero `Y_j`, put `u_j=X_j/Y_j`. Then
`A_j/Y_j=u_j+i`, and `1,u_1,...,u_m` are also linearly independent
over `K_small`. Their normalized heights are positive. Thus the
remaining arithmetic can be viewed as the common-divisor problem for
the functions `u_j+i` and `u_j-i`, with small residues as coefficients.

This does not make the `u_j` algebraic over `K_small`. In fact their
positive heights prohibit that. Nor does linear independence force
their generated field to have transcendence degree one, or equip it
with a controlled curve, genus, or set of places.

## 4. What the moving-target theorems actually allow

[Vojta's exposition, Theorem 7.6](https://library.slmath.org/books/Book37/files/vojta.pdf)
permits moving hyperplanes whose coefficient heights are negligible
relative to the point heights, with linear nondegeneracy over the
coefficient field. Its set of places `S` is one fixed finite set.
The new independence statement addresses the relevant coefficient-field
issue, but does not provide this fixed set of places.

[Grieve--Wang, Theorem 1.2](https://arxiv.org/pdf/1902.09109)
likewise permits moving polynomials of fixed positive degrees and small
coefficient heights. Its arguments are units outside a single fixed
finite `S`. The alternative involving proper torus subgroups and
small-height translates is therefore not available merely because
the residues are small: conductor prime supports may change with
the sequence.

### Dropping fixed support gives a false theorem

Let `p_n` be increasing primes and put

```text
u_n=p_n,             v_n=2p_n-1,
W_n=log p_n,
f(X,Y)=X-1,         g(X,Y)=Y-1.
```

The coefficients are constant, and

```text
gcd(f(u_n,v_n),g(u_n,v_n)) = p_n-1,
max(h(u_n),h(v_n)) = W_n+O(1).
```

Each pair is an `S_n`-unit pair for its finite set of prime divisors,
but no fixed finite `S` works even on an infinite subsequence because
`u_n` is a fresh prime.

Also `gcd(u_n,v_n)=1`. For every fixed pair of integers `(a,b)!=0`,
reduction of numerator and denominator gives

```text
h(u_n^a v_n^b) = c_(a,b) W_n+O_(a,b)(1),
```

where `c_(a,b)=|a|+|b|` if their nonzero signs agree, and
`c_(a,b)=max(|a|,|b|)` if their signs differ. In particular it is
strictly positive.

Thus these points cannot belong to a finite union of fixed proper
torus subgroups translated by negligible-height points. Every proper
torus subgroup has a nonzero fixed character vanishing on it, and a
small-height translate would make that character have height `o(W)`.
The gcd remains a full positive fraction of the height nevertheless.

This is an explicit counterexample to the corresponding unrestricted
moving-support gcd dichotomy, even with constant polynomial
coefficients. It is not an endpoint-circle counterexample.

## 5. Why replacing the prime support by cut labels does not fix it

A fixed tuple has only finitely many cut labels, but a cut block may
contain many prime factors. Their aggregate divisibility cannot be
treated as one valuation.

For an explicit Gaussian example take `pi=2+i`, `rho=3+2i`, and
the block `H_n=pi^n rho^n`. Normalize by `W_n=n(log 5+log 13)`.
The proposed aggregate valuation

```text
v_block(x)=lim (v_pi(x_n) log 5+v_rho(x_n) log 13)/W_n
```

assigns positive values to both `x_n=pi^n` and `y_n=rho^n`.
But `x_n+y_n` is a unit at both primes, so its value is zero.
This violates the valuation inequality
`v(x+y)>=min(v(x),v(y))`.

Hence one cannot convert the bounded number of cut patterns into a
bounded set of function-field places by aggregating their primes.
Keeping the individual places restores the unbounded-support issue.

There is a separate failure if one only keeps limits of the original
fixed places: for the sequence of fresh primes `p_n` and normalization
`log p_n`, every fixed finite place has limiting valuation zero, while
the archimedean growth has exponent one. The balancing mass has escaped
to moving places. A product-formula model must retain that mass rather
than silently discard it.

## 6. Even a valid one-variable model does not automatically bound R

The counterexample in Section 4 has the exact function-field model

```text
u=T,             v=2T-1
```

over `K_small`. The element `T=[p_n]` is transcendental over this
field by its positive height and relative algebraic closedness.
Both functions have degree one, and

```text
gcd(T-1,2T-2)=T-1
```

also has degree one. This is a perfectly ordinary rational curve.
Their function-field degrees stay fixed while the arithmetic
specialization heights grow like `log p_n`.

Accordingly, a function-field estimate with an error depending on the
number of places or on the curve does not become an `o(W)` arithmetic
error merely by passing to a sequence field. In this example one
function-field degree corresponds to an entire positive unit of
normalized arithmetic height. A proof would need a quantitative
comparison and a strict degree inequality adapted to the reconstruction
system, not just existence of a function-field presentation.

For the actual cluster, a field generated over `K_small` by finitely
many conductor products certainly exists. Nothing here makes it a
controlled curve, and the small residues do not force its projective
invariants into `K_small`: the independently proved cross-ratio
heights are positive fractions of `W`.

## 7. Exact remaining task

The new positive statement is the coefficient lower bound (3), and
its consequence (4). It supplies nondegeneracy for the small-height
coefficient viewpoint, uniformly through arbitrary prime exponents.

The failed step is now specific. Existing moving-target gcd theorems
need a fixed place support; a general version dropping it is false.
The alternative sequence-field passage needs an additional arithmetic
comparison preserving all gcd mass and producing a strict inequality
at the endpoint. Neither replacing primes by cut labels nor simply
calling the coefficients constant supplies that comparison.

No fixed positive lower exponent for the primitive residues has been
proved, and uniformity in the radius remains open in this work.

## 8. Multilinear and binary-support polynomial independence

The linear result extends further than the independent-block model.
It also holds for multilinear polynomials through arbitrary prime
powers. The decisive choice is a support of **minimum cardinality**,
not just an inclusion-minimal support.

More generally, fix positive integers `q_1,...,q_m`, and put
`q_min=min_j q_j`. Consider a nonzero polynomial

```text
P(A)=sum_(J subset [m]) c_J product_(j in J) A_j^(q_j),
                       c_J in Z[i].
```

If `P(A)=0`, then

```text
max_J log max(1,|c_J|)
    >= q_min (1-eta) W/2^k.                         (5)
```

The same conclusion applies after multiplying all monomials by one
common monomial. Thus it covers any support having only two allowed
exponents in each variable, with arbitrary fixed positive steps between
them. Taking every `q_j=1` proves independence of all `2^m`
squarefree monomials over `K_small`.

### Proof retaining the nested prime allocations

Among nonzero coefficients choose `c_J` with minimal weighted degree
`sum_(j in J) q_j`. If `J=[m]`, no other support can occur, and a
single nonzero monomial cannot vanish. Hence assume `J` is proper.

At a split prime, define the two ordering gaps

```text
g_p^+(J)=max(0, min_(i notin J) a_i
                 -max(a_0,max_(j in J) a_j)),

g_p^-(J)=max(0, min(a_0,min_(j in J) a_j)
                 -max_(i notin J) a_i).
```

For `J=empty`, the inside maximum and minimum are omitted, leaving
`a_0` in each expression. These gaps count the two threshold
patterns where the anchor and every row of `J` are together, with
all other rows opposite.

Suppose `g=g_p^+(J)>0`, and put `x_i=v_pi(A_i)` and
`L=max_(j in J) x_j`, with `L=0` for an empty `J`. Then

```text
x_j<=L             for j in J,
x_i>=L+g           for i notin J.
```

For any other occurring support `U`, put

```text
r=sum_(i in U\J) q_i,
s=sum_(j in J\U) q_j.
```

Minimal weighted degree gives `r>=s`; moreover `r>=q_min` because
`U!=J` cannot be a proper subset of `J`. Therefore

```text
v_pi(A_U)-v_pi(A_J)
 >= r(L+g)-sL
 >= q_min g.
```

After division by the valuation of the target monomial, every other
term of the relation is divisible by `pi^(q_min g)`. Hence this
power divides `c_J`. The same argument at the conjugate prime uses
`g_p^-(J)`.

Summing these forced coefficient valuations gives

```text
log|c_J| >= (q_min/2)
             sum_p (g_p^+(J)+g_p^-(J)) log p.
```

The two indicated threshold patterns each have weight at least
`(1-eta)W/2^k`. This proves (5).

This argument explicitly permits nonzero valuations among the inside
rows. Equal-degree monomials exchange equally many weighted factors;
higher-degree monomials add at least as much weight as they remove.
That is why nested thresholds do not invalidate the proof.

Exact finite checks covered 141,480 such weighted support comparisons,
using two and three nonreference rows, allocation values from zero
through four, weights `q_i` from one through three, and both conjugate
valuation directions. All satisfied the displayed valuation gap. The
general proof does not depend on the finite ranges.

### All homogeneous quadratic monomials are also independent

There is a further degree-two consequence beyond binary supports.
Every nonzero homogeneous quadratic relation

```text
sum_i b_ii A_i^2 + sum_(i<j) b_ij A_i A_j = 0
```

has some coefficient satisfying the lower bound in (3).

If a square coefficient `b_jj` is nonzero, use the ordering gap for
`J={j}`. The target monomial `A_j^2` has valuation at least one gap
less than every other degree-two monomial. The same divisibility
argument forces `b_jj` to have logarithmic size at least
`(1-eta)W/2^k`. If every square coefficient is zero, the relation is
multilinear and (5) applies.

Thus neither an arbitrary multilinear relation nor a homogeneous
quadratic relation with negligible-height coefficients can occur in
the extracted sequence.

### Where this particular isolation argument stops

This does not establish arbitrary polynomial nondegeneracy. Already
the inhomogeneous quadratic support

```text
{x,x^2,y,y^2}
```

has no unique least-weight monomial for any of the binary cut weights:
at weights `(1,0)` the two `y` monomials tie, at `(0,1)` the two
`x` monomials tie, and at `(1,1)` the two linear monomials tie.
For homogeneous cubics, the six permutations of the exponent vector
`(2,1,0)` similarly give a tie at every binary cut.

These are exact support obstructions to the unique-minimum argument,
not asserted Gaussian solutions. They identify precisely why the
proved binary-support and quadratic statements cannot simply be
renamed a theorem for all bounded-degree polynomials.

## 9. Stronger nondegeneracy still does not remove the fixed-S requirement

The primary agent supplied a stronger version of the counterexample in
Section 4, which was independently checked. For any prescribed
separate-degree cap `D>=1`, take

```text
t_n=2 n!,
u_n=t_n^(D+1)+1,
v_n=t_n^(D+2)+1.
```

Then

```text
gcd(u_n-1,v_n-1)=t_n^(D+1),
max(h(u_n),h(v_n))=(D+2)log t_n+o(1).
```

The gcd is asymptotically the fraction `(D+1)/(D+2)` of the height.
Also `t_n u_n-v_n=t_n-1`, so `gcd(u_n,v_n)` divides two; both
numbers are odd, and hence that gcd is one.

Every monomial `u^i v^j`, with `0<=i,j<=D`, has leading degree
`i(D+1)+j(D+2)` in `t`. These degrees are distinct: equality of two
such degrees would make `D+1` divide the difference of their `j`
indices, whose absolute value is at most `D`.

Consequently these `(D+1)^2` monomials are linearly independent over
the field of sequences of height `o(log t_n)`. A coefficient of that
height has absolute value `t_n^(o(1))`, and its inverse does too if
it is nonzero. Thus the uniquely largest leading degree in any proposed
finite relation cannot be canceled by the remaining terms.

Finally, every rational prime at most `n` divides `t_n` and divides
neither `u_n` nor `v_n`. Their prime supports eventually avoid every
fixed finite set. This proves support escape without a theorem about
prime factors of polynomial values.

This example satisfies arbitrarily prescribed finite-degree
small-coefficient nondegeneracy while retaining a large gcd and
moving prime support. It therefore shows that the new polynomial
independence statements, although useful positive information, do not
by themselves repair the missing fixed-S hypothesis. As before,
these pairs are not asserted to arise from endpoint circle clusters.

## 10. An explicit all-place failure of a rational-function transfer

The rounded system studied by the primary agent has independent
effective Gaussian blocks, small-height row factors, and identities

```text
A_i=product_(S contains i) H_S,
V_i A_i-conj(V_i A_i)=2i y_i,
h(V_i), log max(1,|y_i|)=o(W).
```

Here `S` ranges over every subset of `[m]`, including the empty set.
All block norm weights are approximately equal. The original points
give these moving linear forms; rounded points themselves need not
remain in the original arc.

There is a precise false transfer to avoid: representing all quantities
by rational functions on a curve and preserving their orders at infinity
does not preserve the effectiveness of the blocks. This already fails on
the rational curve, without a genus issue.

Fix any `m>=2`, and put `d=2^(m-1)`. For the empty set and every
nonsingleton set `S`, choose monic Gaussian linear polynomials `L_S(t)`
whose roots and conjugate roots are all distinct. Choose them also
coprime to every polynomial `t^d+j+i` and `t^d+j-i`, `1<=j<=m`.
Such a choice exists by avoiding finitely many Gaussian roots at each
step. Set

```text
A_j(t)=t^d+j+i,
P_j(t)=product_(S contains j, |S|>=2) L_S(t),
H_S(t)=L_S(t)                         if |S|!=1,
H_{ {j} }(t)=A_j(t)/P_j(t).
```

There are `d-1` factors in `P_j`. Consequently **every** `H_S`
has a simple pole at infinity and leading term `t`. The row-product
identities hold exactly, and

```text
Im A_j(t)=1                          for real t.
```

If one uses only infinity orders, this looks like the required model:
with `G=product_S H_S`, one has `|G(t)|~t^(2d)`, while the angles
of `A_j/conj(A_j)` from the anchor `1` are asymptotic to `2t^(-d)`.
Thus the formal normalized arc length tends to two for arbitrarily
many rows. The rational circle coordinates are

```text
z_0=conj(G),
z_j=conj(G) A_j/conj(A_j).
```

But each singleton block has finite poles at the zeros of `P_j`.
Its rational-function height is `d`, although its pole order at
infinity is only one. The all-place independent effective-block
hypothesis has therefore been lost.

The missing denominator can be computed exactly. Define

```text
T=product_(|S|>=2) L_S^(|S|-1).
```

Then

```text
G=L_empty (product_j A_j)/T,
deg T=(m-2)d+1.
```

The last degree identity follows from
`sum_S |S|=m 2^(m-1)` and the count of nonsingleton subsets.
Our coprimality choices prevent cancellation with the numerator.
All rational circle coordinates have the common denominator
`conj(T)`. Clearing it gives polynomial coordinates, and removing
their common factor `conj(L_empty)` leaves the usual coordinates
of radius

```text
R_min(t) comparable to t^(m d),
Delta(t) sqrt(R_min(t)) comparable to t^((m/2-1)d).
```

The radius assertion is also valid up to constant factors at positive
integer specializations. Indeed, differences of the row numerators
are the fixed nonzero integers `j-l`; differences with their conjugates
are the fixed nonzero Gaussian integers `j-l+2i`; and
`gcd(A_j,conj(A_j))` divides `2i`. Hence all Gaussian gcd corrections
in the exact lcm formula are bounded independently of `t`. Each row
numerator has modulus asymptotic to `t^d`, so the least radius is
asymptotic to `t^(m d)` within fixed positive constant factors.

For `m>=3`, the genuine normalized length grows as a positive power
of `t`. The apparent endpoint family arose entirely from retaining
the archimedean orders while discarding the finite denominator budget.
This is an explicit failed transfer, not an example satisfying the
independent effective-block system.

## 11. The effective polynomial model survives at three circle points

Effectiveness alone does not immediately contradict the endpoint.
There is an elementary full-pattern polynomial example with two
nonreference rows. Take

```text
H_{12}=t+i,
H_{1}=t+1-i,
H_{2}=t+2-i,
H_empty=t+4+2i.
```

These four Gaussian polynomials and their four conjugates have
distinct roots. In particular every block is effective, and all
required pairwise and conjugate gcds are one. Their row products are

```text
A_1=H_1 H_{12}=t^2+t+1+i,
A_2=H_2 H_{12}=t^2+2t+1+2i.
```

With `G=product_S H_S`, the genuine polynomial circle coordinates
`z_0=conj(G)` and `z_j=conj(G) A_j/conj(A_j)` are Gaussian
polynomials. They have common radius `R(t)=|G(t)|~t^4`, and the
three points lie on an arc with angular span asymptotic to `4/t^2`.
Thus `Delta sqrt(R)` tends to four. There are no rational-function
denominators in this construction. Its common factor `conj(H_empty)`
can be removed, but doing so changes the inherited angular scale;
one should not silently renormalize the original endpoint condition.

This example rules out an attempted function-field argument claiming
that the exact effective full-pattern identities are impossible for
every fixed number of rows. It leaves open an exclusion for a
sufficiently large fixed number of rows.

The primary agent's stronger central-truncation reduction produces
row prefactors that are integral and have equal norm, with negligible
height. In a faithful polynomial model where negligible height becomes
degree zero, these are complex constants of equal modulus. After
anchoring they give constant phases in the same identities above.
Thus the equal-norm condition is useful exact arithmetic information,
but by itself does not supply the missing faithful all-place transfer
or change this elementary polynomial example. No function-field
exclusion for sufficiently many rows is proved here.
