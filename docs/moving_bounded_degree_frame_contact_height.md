# Slowly moving bounded-degree frames still fail the contact height test

The fixed-frame obstruction remains valid when the rational frame map
is allowed to vary with the multiple `NP`, provided its representation
degree is uniformly bounded and its coefficient height is negligible
compared with `N^2`.  Two uniform estimates are needed.  A moving
resultant has logarithmic size controlled by the coefficient height,
and the height comparison for the moving projective map must retain the
degree of its pole divisor.  The latter comparison has a square-root
error coming from the moving degree-zero class in `Pic(E)`.

This result applies only to bounded-degree rational formulas with
subquadratic coefficient height.  It does not control arbitrary
pointwise frames and does not prove a uniform lattice-arc bound.

## 1. Bounded representations and the statement

Fix a Weierstrass model of `E/Q`, a rational origin `O`, a nontorsion
point `P`, and the twenty-five rational contact points of
[the disjoint elliptic contact profile](disjoint_mixed_elliptic_integer_contact_profile.md).
Set

```text
Q_N=NP,       h_N=hat h_[O](P)N^2.
```

For each contact `B`, its good-prime contact integer satisfies

```text
log d_B(N)=h_N+O(N).                                  (1)
```

A rational function has representation degree at most `D` and
coefficient height at most `H` if it is the ratio of two forms of
degree at most `D` in the fixed Weierstrass coordinates, restricted to
`E`, and the **joint** projective logarithmic height of the concatenated
numerator and denominator coefficient vectors is at most `H`.  The
joint normalization retains their relative scale: the presentation
`M/1` for an integer `M` has height `log max(1,|M|)`, rather than height
zero as it would if the two one-term vectors were normalized
separately.  For rational `M`, this is the usual joint
numerator-denominator height.  This condition is essential for the
moving-resultant estimate below.  Homogenizing or
performing the bounded additions and multiplications used below keeps
the degree bounded in terms of `D` and changes the joint coefficient
height to at most `O_D(H+1)`.  Local numerator-denominator cancellation
at a contact is handled by jets in Section 2; no unique factorization
or principal-gcd assertion in the coordinate ring of `E` is assumed.

For every relevant `N`, let

```text
z_N=a_N+i b_N,       w_N=c_N+i d_N,
a_N,b_N,c_N,d_N in Q(E),
a_N d_N-b_N c_N=1.                                  (2)
```

Assume all four entries have representation degree at most one fixed
`D`, coefficient height at most `H_N`, and

```text
H_N=o(h_N).                                           (3)
```

Allow a varying positive common denominator `q_N` at `Q_N` and put

```text
A_N=q_N a_N(Q_N),   B_N=q_N b_N(Q_N),
C_N=q_N c_N(Q_N),   D_N=q_N d_N(Q_N).                 (4)
```

Suppose these are integers and

```text
log q_N=o(h_N),
log max(|q_N|,|A_N|,|B_N|,|C_N|,|D_N|)
   <=4h_N+o(h_N).                                     (5)
```

Assign to each contact the same one of the four combinations as in the
fixed-frame theorem:

```text
contact direction       number       u_(B,N)
infinity                   7          z_N
0                          7          w_N
1                          7          z_N+w_N
rho=-3/2                   4          -3z_N+2w_N.      (6)
```

Finally, allow positive error integers `e_B(N)` with

```text
max_B log e_B(N)=o(h_N),
d_B(N) divides e_B(N) Norm(q_N u_(B,N)(Q_N)).          (7)
```

**Theorem.** For fixed `D`, conditions (1)--(7) hold for only finitely
many `N`.

## 2. Uniform resultant bound at a fixed contact

Fix one of the rational points `B` and put

```text
f_(B,N)=Norm(u_(B,N)) in Q(E).                         (8)
```

Its representation degree is bounded in terms of `D`, and its
joint numerator-denominator coefficient height is `O_D(H_N+1)`.

Suppose `ord_B(f_(B,N))<=0`.  Fix once and for all a local uniformizer
`t` at the rational point `B` and local trivializations for the finite
list of ambient line bundles allowed by the degree bound.  Expand any
bounded-degree numerator and denominator presentation as

```text
p=p_r t^r+p_(r+1)t^(r+1)+...,
q=q_s t^s+q_(s+1)t^(s+1)+...,
p_r q_s !=0.
```

The possible orders `r,s` are uniformly bounded.  The jet coefficients
are fixed linear evaluations of the global coefficient vectors, so
their joint heights are `O_D(H_N+1)`.  The inequality
`ord_B(f_(B,N))=r-s<=0` says `r<=s`.  Cancelling the common local power
`t^r` leaves a numerator whose nonzero constant coefficient is `p_r`.
A bounded local Bezout identity between `t` and this reduced numerator
therefore has coefficient height `O_D(H_N+1)`.  Equivalently, their
local resultant is the nonzero leading coefficient `p_r`, with the
fixed trivialization denominators included.  It gives an integer
`R_(B,N)` satisfying

```text
log R_(B,N)<=C_D(H_N+1),                               (9)
```

and, for every rational point `Q`,

```text
sum_p min(k_(B,p)(Q),max(0,v_p(f_(B,N)(Q)))) log p
   <=log R_(B,N).                                     (10)
```

All primes, including the finitely many primes of bad reduction,
trivialization denominators, and vertical coefficient factors, are
included here.  The local Bezout identity bounds the simultaneous
contact and numerator valuation, so no new avoidance condition is
imposed on the recurrence subsequence.  This argument uses only a
bounded Taylor jet at the fixed point `B`; it does not require a global
gcd of divisors on `E`.

Because `q_N` clears both real components of `u_(B,N)(Q_N)`, its norm is
an integer and

```text
v_p(Norm(q_Nu_(B,N)(Q_N)))
 =2v_p(q_N)+v_p(f_(B,N)(Q_N)).                         (11)
```

Taking valuations in (7), replacing the positive last term in (11) by
its minimum with the contact depth, and applying (10) yields

```text
log d_B(N)
 <=log e_B(N)+2log q_N+C_D(H_N+1)
 =o(h_N).                                             (12)
```

This contradicts (1).  Hence every contact satisfies

```text
ord_B Norm(u_(B,N))>0.                                (13)
```

The real components of `u_(B,N)` belong to `Q(E)`.  At the rational
point `B`, their leading coefficients are real, and a sum of their
squares cannot have leading cancellation.  If neither is identically
zero, then

```text
ord_B Norm(u_(B,N))
 =2min(ord_B Re(u_(B,N)),ord_B Im(u_(B,N))).           (14)
```

The evident one-component version applies if one component is zero;
both cannot be identically zero because of (2).  It follows from (13)
that both components of every assigned `u_(B,N)` vanish at `B`.

Choose a complementary rational combination `v_(B,N)` of `z_N,w_N`
whose coefficient determinant with `u_(B,N)` is a fixed nonzero
rational number.  For the four rows of (6), one may use `w_N,z_N,w_N,z_N`.
Equation (2) gives

```text
det(u_(B,N),v_(B,N))=nonzero constant.                 (15)
```

Since both components of `u_(B,N)` vanish at `B`, at least one
component of `v_(B,N)`, and hence at least one of
`a_N,b_N,c_N,d_N`, has a pole there.

Let

```text
Delta_N=sum_X max(0,-ord_X(a_N),-ord_X(b_N),
                      -ord_X(c_N),-ord_X(d_N))[X].     (16)
```

The twenty-five contacts are distinct, so

```text
delta_N=deg Delta_N>=25.                               (17)
```

The representation-degree bound also gives a uniform upper bound for
`delta_N`.

## 3. Uniform evaluation height for a moving map

We record the height estimate that prevents the moving coefficients
from cancelling (17).  A standalone version, including the matching
upper bound and explicit treatment of all small degrees, is
[the uniform moving-map height lemma](moving_frame_uniform_height_lemma.md).

**Uniform height lemma.** Fix `E`, its Weierstrass model, `D`, and the
number `r` of functions.
Let `f_1,...,f_r in Q(E)` have representation degrees at most `D` and
coefficient heights at most `H`.  Let `Delta` be their least common pole
divisor and define

```text
Phi=[1:f_1:...:f_r]:E -> P^r.
```

Then, for every rational `Q` at which the displayed affine values are
defined,

```text
h(Phi(Q))
 >=(deg Delta) hat h_[O](Q)
   -C_(D,r)(sqrt(hat h_[O](Q)(H+1))+H+1).             (18)
```

Here and below the constant may also depend on the fixed curve and
embedding, but not on the functions, their coefficients, or `Q`.

Here is a concrete proof that also handles moving base divisors.  Write

```text
Delta=sum_j m_j[B_j],       delta=deg Delta.
```

If `delta=0`, all the functions are constant and (18) is immediate.
The case `delta=1` cannot occur: a degree-one line bundle on a
genus-one curve has only one independent global section, whereas a
nonconstant projective map requires at least two.  Assume `delta>=2`
for the rest of the proof.

Each pole, with at least its pole multiplicity, lies in the zero divisor
of one of the bounded-degree denominators.  Arithmetic Bezout on the
fixed Weierstrass model gives

```text
sum_j m_j hat h_[O](B_j)<=C_(D,r)(H+1).                (19)
```

This includes complete Galois orbits of nonrational pole points.  Put
`S=sum_j m_jB_j` in the elliptic group.  The divisor is rational, so
`S in E(Q)`, and Cauchy--Schwarz gives

```text
hat h_[O](S)
 <=delta sum_j m_jhat h_[O](B_j)
 <=C_(D,r)(H+1).                                      (20)
```

Choose a point `T` with `[delta]T=S`.  Since `delta` is bounded in
terms of `D` and `r`, `T` has uniformly bounded degree over `Q`, and

```text
hat h_[O](T)=hat h_[O](S)/delta^2<=C_(D,r)(H+1).       (21)
```

Translate the functions by putting `f'_i(X)=f_i(X+T)`.  Their common
pole divisor is

```text
Delta'=sum_j m_j[B_j-T].
```

The sum of its points in the elliptic group is zero, so
`Delta'` is linearly equivalent to `delta[O]`.  There is therefore a
rational function `g`, defined over the bounded-degree field of `T`,
with

```text
div(g)=Delta'-delta[O].                                (22)
```

All of these objects have uniformly bounded algebraic degree.  Their
coefficient heights are `O_(D,r)(H+1)`: the addition law controls the
translated functions, and `g` is obtained by imposing the prescribed
bounded-order zeros on the fixed vector space `L(delta[O])`.  In a
fixed basis, that is a bounded-size linear system whose entries have
height `O_(D,r)(H+1)`; a nonzero maximal minor and Cramer's rule give the
same bound for a solution.

The tuple

```text
(g,gf'_1,...,gf'_r) in L(delta[O])^(r+1)               (23)
```

is base-point-free.  Away from `Delta'`, its first entry is nonzero.
At a point of `Delta'`, at least one function attaining the maximal
pole order in `Delta'` leaves a nonzero product with `g`.  For
`delta=2`, use the usual rank-two basis of `L(2[O])`; for `delta>=3`,
use any fixed basis of the `delta`-dimensional space `L(delta[O])`.
In that basis, the tuple (23) has coefficient height
`O_(D,r)(H+1)`.  The effective projective Nullstellensatz for this
base-point-free tuple bounds all nonarchimedean specialization gcds by
`C_(D,r)(H+1)` logarithmically; the archimedean coefficient comparison has
the same bound.  Since the standard height for `delta[O]` is
`delta hat h_[O]+O_(D,r)(1)`, it follows that

```text
h([g(X):gf'_1(X):...:gf'_r(X)])
 >=delta hat h_[O](X)-C_(D,r)(H+1).                    (24)
```

Take `X=Q-T`.  The projective point on the left of (24) is `Phi(Q)`.
The Neron--Tate pairing and (21) give

```text
hat h_[O](Q-T)
 >=hat h_[O](Q)-2sqrt(hat h_[O](Q)hat h_[O](T))
 >=hat h_[O](Q)-C_(D,r) sqrt(hat h_[O](Q)(H+1)).       (25)
```

Equations (24)--(25) prove (18).  This also explains why a fixed-map
`O(N)` error cannot simply be quoted when the map moves: the translating
point `T`, equivalently the degree-zero class of the pole divisor,
moves with the map, and (19)--(21) provide the needed uniform control.

## 4. The contradiction

Apply the uniform height lemma to

```text
Phi_N=[1:a_N:b_N:c_N:d_N].
```

Equations (3), (17), and (18) give

```text
h(Phi_N(Q_N))
 >=delta_N h_N
   -C_D(sqrt(h_N(H_N+1))+H_N+1)
 >=25h_N-o(h_N).                                      (26)
```

But (4) gives integral homogeneous coordinates

```text
[q_N:A_N:B_N:C_N:D_N]
```

for the same point.  Primitive reduction only decreases their maximum,
so (5) gives

```text
h(Phi_N(Q_N))<=4h_N+o(h_N),                            (27)
```

contradicting (26).

Thus slowly varying coefficients do not evade the fixed-frame
obstruction at bounded representation degree.  No unbounded-degree
conclusion is asserted without controlling the constants that depend
on the presentation degree.  More conditionally, choose a
nondecreasing function `C(D)>=1` dominating the resultant and uniform
height constants above.  The same proof permits degrees `D_N` if

```text
C(D_N)^2(H_N+1)=o(h_N),
```

in addition to the denominator and error bounds already stated.  This
condition does not give such control for arbitrary growing degrees.  A
coefficient height comparable to `h_N`, or a frame chosen without a
complexity-controlled rational formula, can likewise carry a full
contact-scale contribution.  No claim is made about those cases or
about the global uniform `sqrt(R)` problem.
