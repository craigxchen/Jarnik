# Fixed rational frames cannot carry the elliptic contact profile

The disjoint elliptic contact family cannot be produced by evaluating
one fixed rational determinant-one frame on the elliptic source at the
multiples `NP`, even if a common denominator and the failures of the
contact divisibilities have subpower size.  The obstruction is
quantitative: the twenty-five prescribed contacts force the fixed
projective frame map to have degree at least twenty-five, whereas a
short first frame has logarithmic height at most four times the
elliptic height.

This is a statement about a **fixed** rational map.  The auxiliary
common denominators and error integers may vary with `N` subject to the
subpower bounds below.  The result does not apply when the frame or its
rational functions vary with `N`, and it gives no global uniform
lattice-arc bound.

## 1. A denominator-tolerant statement

Use the fixed elliptic curve, nontorsion point `P`, and contact points
from
[the integer contact profile](disjoint_mixed_elliptic_integer_contact_profile.md).
Put

```text
h_N=hat h_[O](P) N^2.
```

For each of the twenty-five distinct rational contact points `B`, let
`d_B(N)` be its good-prime contact integer on the compact-recurrence
subsequence.  Thus

```text
log d_B(N)=h_N+O(N).                                  (1)
```

The four contact directions and their assigned Gaussian linear forms
are

```text
direction       number of contacts       u_B
infinity                 7                z
0                        7                w
1                        7                z+w
rho=-3/2                 4                -3z+2w.       (2)
```

Let

```text
z=a+ib,       w=c+id,
a,b,c,d in Q(E),       ad-bc=1.                       (3)
```

Suppose that, for infinitely many indices in the recurrence
subsequence, there is a positive common denominator `q_N` for the four
values at `Q_N=NP`.  Write

```text
A_N=q_N a(Q_N),   B_N=q_N b(Q_N),
C_N=q_N c(Q_N),   D_N=q_N d(Q_N).                     (4)
```

Assume these are integers and

```text
log q_N=o(h_N),
log max(|q_N|,|A_N|,|B_N|,|C_N|,|D_N|)
    <=4h_N+o(h_N).                                    (5)
```

Fixed denominators in the coefficients of (2) are harmless and may be
absorbed into the fixed bad-prime set.  To allow subpower failures of
the norm divisibilities, assume for every contact `B` that there is a
positive integer `e_B(N)` such that

```text
log e_B(N)=o(h_N),
d_B(N) divides e_B(N) Norm(q_N u_B(Q_N)).              (6)
```

This formulation is weaker than divisibility of the primitive norm
numerator: it permits the whole common denominator in (4) to
contribute to the norm and permits an additional error of subpower
size.

**Theorem.** Conditions (1)--(6) cannot hold for infinitely many `N`.

## 2. A full contact cannot come from an accidental numerator gcd

We use the following elementary fixed-resultant fact.  Let `f` be a
fixed nonzero rational function on `E` and let `B` be a fixed rational
point.  If

```text
ord_B(f)<=0,                                           (7)
```

then the zero divisor of `f` and `B` are disjoint on the generic fiber.
After choosing one proper integral model, their closures have a fixed
finite arithmetic intersection. Include the fixed vertical coefficient
factors of `f` as well: a constant such as `f=p^A` has no generic zero
divisor but contributes a bounded valuation at `p`. Consequently there is a fixed
positive integer `R_(f,B)` such that, for every rational point `Q`,

```text
sum_p min(k_(B,p)(Q),max(0,v_p(f(Q)))) log p
   <=log R_(f,B).                                      (8)
```

Here `k_(B,p)(Q)` is the intersection depth of `Q` with `B`; the
finitely many model and denominator primes are included in the fixed
constant.  In local coordinates this is just the resultant bound for
a local equation of `B` and a reduced numerator of `f`.  Away from the
resultant primes, write `f=t^r epsilon` with `r=ord_B(f)<=0` and
`epsilon` a unit on the residue disk of `B`; then a point meeting `B`
has `v_p(f(Q))<=0`.  At a resultant prime, the minimum in (8) is bounded
by the fixed intersection multiplicity of the two closures.  This
formulation does not require a new recurrence condition or a later
enlargement of the avoided-prime set.

Fix one contact `B` and set

```text
f_B=Norm(u_B) in Q(E).                                 (9)
```

This function is not identically zero.  Otherwise the two real
components of the nonzero rational linear combination `u_B` would both
vanish identically, contradicting the invertibility in (3).  If
`ord_B(f_B)<=0`, then at every prime,

```text
v_p(Norm(q_N u_B(Q_N)))
 =2v_p(q_N)+v_p(f_B(Q_N)).                            (10)
```

The left side is nonnegative because (4) clears the two components of
`u_B`; fixed coefficient denominators only change a fixed constant.
At a prime outside the resultant, the last term in (10) is nonpositive
when that prime contributes to `d_B(N)`.  At the remaining primes, its
possible overlap with the contact depth is bounded by (8).  More
explicitly, from

```text
k_(B,p)(Q_N)
 <=v_p(e_B(N))+2v_p(q_N)+max(0,v_p(f_B(Q_N)))
```

one may replace the last term by its minimum with the left side.
Taking valuations in (6), multiplying by `log p`, and summing therefore
gives

```text
log d_B(N)
 <=2log q_N+log e_B(N)+log R_(f_B,B)
 =o(h_N),                                             (11)
```

contrary to (1).  Hence every prescribed contact forces

```text
ord_B Norm(u_B)>0.                                     (12)
```

This is the point at which fixedness matters.  The possible overlap at
all exceptional primes is bounded by one fixed resultant only because
both `f_B` and `B` are fixed.  The common denominator and error terms in
(5)--(6) are too small to carry a contact of logarithmic size
`h_N+O(N)`.

## 3. Real coefficients turn norm zeros into two simultaneous zeros

Write the real components of `u_B` as `r_B,s_B in Q(E)`.  At the
rational point `B`, a local parameter and all leading coefficients may
be taken over `Q`.  If neither component is identically zero, then

```text
ord_B(r_B^2+s_B^2)
   =2 min(ord_B(r_B),ord_B(s_B)).                      (13)
```

When the two orders agree, the leading coefficient is a sum of two
rational squares and is nonzero; when they differ, the term with
smaller order is unique.  The same formula has the evident
interpretation if one component is identically zero.  Both components
cannot vanish identically by (3).  Therefore (12) implies

```text
ord_B(r_B)>0,       ord_B(s_B)>0.                      (14)
```

Thus the assigned real two-vector `u_B` vanishes at every one of its
contact points.  This rules out poles as well as nonzero regular
values; denominator clearing cannot turn either into a fixed-function
zero.

## 4. Determinant one forces a pole at every contact

For each of the four rows in (2), choose a fixed complementary real
rational combination

```text
u_B=alpha z+beta w,       v_B=gamma z+delta w,
alpha delta-beta gamma !=0.                            (15)
```

For example, use `w,z,w,z` as complements to
`z,w,z+w,-3z+2w`, respectively.  The determinants of the two real
component vectors satisfy the rational-function identity

```text
det(u_B,v_B)
  =(alpha delta-beta gamma)(ad-bc)
  =alpha delta-beta gamma !=0.                         (16)
```

By (14), both components of `u_B` vanish at `B`.  If both components of
`v_B` were regular there, the left side of (16) would also vanish.
Consequently `v_B` has a pole at `B`, and hence at least one of
`a,b,c,d` has a pole at `B`.

Let the least common pole divisor of the four entries be

```text
Delta=sum_X max(0,-ord_X(a),-ord_X(b),-ord_X(c),-ord_X(d)) [X].  (17)
```

All twenty-five contact points are distinct, and the preceding
argument puts each of them in `Delta` with multiplicity at least one.
Therefore

```text
deg Delta>=25.                                         (18)
```

The five sections `1,a,b,c,d` of `O_E(Delta)` have no common zero.  Off
`Delta`, the section `1` is nonzero; at a point of `Delta`, at least one
entry attaining the maximal pole order in (17) is nonzero after the
common local pole is cleared.  Hence they define a morphism

```text
Phi:E -> P^4,       Q |-> [1:a(Q):b(Q):c(Q):d(Q)]
```

with

```text
Phi^*O(1)=O_E(Delta).                                  (19)
```

Height comparison on the elliptic curve now gives

```text
h(Phi(NP))=(deg Delta) h_N+O(N)>=25h_N+O(N).           (20)
```

On the other hand, (4) supplies integral homogeneous coordinates

```text
[q_N:A_N:B_N:C_N:D_N]
```

for the same projective point.  Dividing by their common gcd can only
decrease the maximum, so (5) gives

```text
h(Phi(NP))<=4h_N+o(h_N),                               (21)
```

contradicting (20).  Notice that the cleared matrix has determinant
`q_N^2`; it need not itself lie in `SL_2(Z)`.  The proof only uses the
determinant-one identity for the underlying rational map and the
integral coordinates in (4), so allowing a subpower common denominator
is legitimate.

Already the seven infinity contacts alone force both nonzero
components among `a,b` to vanish at seven distinct points.  Any such
component then has map degree at least seven and height at least
`7h_N+O(N)`, contradicting the bound `4h_N+o(h_N)`.  The determinant
argument above uses all four contact directions and records the
stronger intrinsic statement that the full projective frame map has
degree at least twenty-five.

## 5. Scope

The result strengthens the fixed-section observation in
[the Boolean-fiber note](disjoint_elliptic_boolean_fiber_identities.md):
the obstruction applies to arbitrary fixed rational entries, not only
to low-degree coordinate or Gaussian linear sections.  It is also more
explicit than the finite-cover finiteness result in
[the first-frame barrier](disjoint_elliptic_first_frame_barrier.md),
because it gives a direct height contradiction for a fixed rational
parameterization.

There is no contradiction for a pointwise choice of a new frame at
each `N`.  Such choices need not be evaluations of any fixed rational
map, their zero and pole divisors can move with `N`, and the fixed
resultant argument in Section 2 no longer applies.  Thus this result
does not settle the residual Gram system, exclude higher-dimensional
arithmetic families, or prove a uniform `sqrt(R)` lattice-arc bound.
