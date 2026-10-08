# Five-point CM contacts as an oriented norm-conic lift

The two rational cross-ratios of five directions do not retain the
integral frame, the orientations of the split contact factors, or the
determinant-one condition.  This note gives an exact simultaneous system
which retains all three.  It is a concrete lift problem over the rational
five-point moduli surface, not a proof that such a lift exists and not a
new uniform short-arc bound.

## 1. The determinant-one frame

Apply the median rotation and primitive normalization of
`median_rotation_integer_frame.md` to three anchor rows.  Write the five
resulting Gaussian rows as

```text
Q_i=U v_i,       U=(z w) in SL_2(Z),
v_i=(r_i,s_i) in Z^2,       v_1=(1,0).                 (1)
```

Here a real column `(x,y)` is identified with `x+iy`.  Consequently

```text
Im(bar(z)w)=1.
```

Put

```text
A=N(z),   B=Re(bar(z)w),   C=N(w).
```

Then, exactly,

```text
bar(z)w=B+i,       AC-B^2=1.                           (2)
```

For every row define

```text
X_i=Re(bar(z)Q_i),       s_i=Im(bar(z)Q_i).
```

The second equality agrees with the second coordinate in (1), since
`det(z,Q_i)=s_i`.  Multiplying (1) by `bar(z)` gives the stronger
oriented identity

```text
bar(z)Q_i=X_i+i s_i=A r_i+(B+i)s_i.                   (3)
```

Its norm and real part give

```text
X_i^2+s_i^2=A N(Q_i),
X_i == B s_i (mod A),
r_i=(X_i-Bs_i)/A.                                     (4)
```

Thus the five row norms lie on five norm conics with the same values of
`A,B`, where `B^2 == -1 (mod A)`.  The common oriented factorization
`bar(z)w=B+i` is stronger than the Gram equation `AC-B^2=1`: the latter
forgets which Gaussian primes of `B^2+1` occur in `z` and in `w`.

For two rows, (3) also gives

```text
X_i s_j-s_i X_j=A det(v_i,v_j)=A det(Q_i,Q_j).         (5)
```

Equations (2)--(5) are the determinant-one part of the lift.

## 2. The prescribed cross-ratios with their integral scales

Assume the directions are distinct and `s_2,s_3,s_4,s_5` are nonzero.
With the normalization `(infinity,0,1,a,b)`, the cross-ratios are

```text
a=det(v_1,v_3)det(v_2,v_4)
    /(det(v_1,v_4)det(v_2,v_3)),
b=det(v_1,v_3)det(v_2,v_5)
    /(det(v_1,v_5)det(v_2,v_3)).                       (6)
```

They are equivalent to the integral-scale barycentric identities

```text
s_2 s_3 Q_4
 =s_4((1-a)s_3 Q_2+a s_2 Q_3),
s_2 s_3 Q_5
 =s_5((1-b)s_3 Q_2+b s_2 Q_3).                        (7)
```

For example, if `a=p/q` is reduced, the first equation becomes the
Gaussian-integer equality

```text
q s_2 s_3 Q_4
 =s_4((q-p)s_3 Q_2+p s_2 Q_3).                        (8)
```

Unlike (6), equation (8) visibly retains the row orientations and the
four integral row scales.  Substitution of (5) gives an equivalent
norm-conic form:

```text
a s_4(X_2s_3-s_2X_3)=s_3(X_2s_4-s_2X_4),
b s_5(X_2s_3-s_2X_3)=s_3(X_2s_5-s_2X_5).              (9)
```

Equations (2)--(4), (8), and its `b` analogue parameterize every common
determinant-one lift once the integral scale variables are included.

## 3. Substitution of the oriented CM blocks

Let `chi_i(S)` be the membership bit of row `i` in a Boolean cut and put

```text
mu_S=median(chi_1(S),chi_2(S),chi_3(S)).
```

After the median rotation, the surviving core in row `i` is

```text
A_i' = product_(chi_i(S)-mu_S=1) H_S
       product_(chi_i(S)-mu_S=-1) bar(H_S).            (10)
```

This factorization needs a small correction-dependent trimming.  Before
the median rotation write the corrections as `K_i`, and put

```text
E=product_i K_i,
tilde(H_S)=H_S/gcd_G(H_S,(E bar(E))^4).                 (11)
```

At a Gaussian prime let `epsilon_i` be the signed valuation contributed
by `K_i`.  The median identity and the Lipschitz property of the median
give

```text
|epsilon_i-(median_a(chi_a(S)h+epsilon_a)-mu_S h)|
 <= |epsilon_i|+sum_(a=1)^3 |epsilon_a|,
```

where `h=v_pi(H_S)`.  The fourth power in (11) therefore removes more
than the possible local error.  The remaining `tilde(H_S)` divides the
rotated primitive row in exactly the orientation predicted by
`chi_i(S)-mu_S`.  Since the original blocks have disjoint support, the
product of all removed parts divides `(E bar(E))^4`; its total logarithmic
size is `O(sum_i log max(1,|K_i|))`.

Replace `H_S` below by its trimmed value and set

```text
Q_i=K_i' A_i',       log max(1,|K_i'|)=o(w).           (12)
```

The same local inequality bounds the residual factor `K_i'` by a fixed
multiple of the total correction height.  Thus (12) is an actual
Gaussian divisibility statement, rather than a formal factorization
inferred only from the median height estimate.

The exact oriented norm-conic equations are therefore

```text
bar(z) K_i' A_i'=X_i+i s_i,
X_i^2+s_i^2
 =A N(K_i') product_(chi_i(S)!=mu_S)N(H_S),
X_i == B s_i (mod A).                                 (13)
```

These equations contain information absent from the aggregate rational
contact norms.  In particular, all rows use the same square root `B` of
`-1` modulo `A`, and (12) fixes the actual Gaussian orientation of each
factor rather than only its norm.

For a boundary contact indexed by the complementary cuts `e,e^c`, put

```text
T_e={i:chi_i(e)!=mu_e}.
```

This collision cluster contains at most one of the three anchors.  The
single oriented divisor surviving the median rotation is

```text
G_e=H_e bar(H_(e^c))       if mu_e=0,
G_e=bar(H_e) H_(e^c)       if mu_e=1.                  (14)
```

It has `N(G_e)=N(H_e)N(H_(e^c))`, the retained rational contact norm,
and

```text
r_i z+s_i w == 0 (mod G_e)       for every i in T_e.  (15)
```

Equation (15) is the orientation-sensitive contact condition.  Taking
norms or projective cross-ratios loses it.

## 4. The simultaneous congruence lattice

For each contact choose one primitive coefficient vector
`lambda_e=(r_e,s_e)` from its collision cluster.  All vectors in the
same cluster define the same kernel modulo `G_e`: their determinants are
divisible by `N(G_e)`, and a primitive integer pair is unimodular over
every Gaussian prime dividing `G_e`.  Define

```text
L={ (z,w) in Z[i]^2 : r_e z+s_e w ==0 (mod G_e)
                         for all ten contacts e }.     (16)
```

For one contact the displayed map onto `Z[i]/(G_e)` is surjective, so
its kernel has index `N(G_e)`.  The ten `G_e` have pairwise disjoint
Gaussian support.  The Chinese remainder theorem therefore gives

```text
[Z[i]^2:L]=product_e N(G_e).                           (17)
```

Equivalently, `L` is a rank-two Gaussian module whose determinant ideal
is generated, up to a unit, by `Delta=product_e G_e`.  The integral
Hermitian matrix

```text
J=[[0,-i],[i,0]]
```

represents `2q(z,w)=2Im(bar(z)w)` and has determinant `-1`.  If `M` is a
Gaussian basis matrix for `L`, then

```text
det(bar(M)^t J M)=-N(det M)=-N(Delta)
                 =-product_e N(G_e).                  (18)
```

Equivalently, the determinant of the associated four-variable real
quadratic form scales by `[Z[i]^2:L]^2`.  The required common lift is a
vector of `L` representing `q=1`, together with (7) and (13).

There is also a local reciprocity constraint.  At an odd prime
`pi|G_e`, (15) makes `(z,w)` a multiple of `(-s_e,r_e)` modulo `pi`.
Using

```text
bar(z)w-z bar(w)=2i
```

shows that `r_e bar(z)+s_e bar(w)` is a unit modulo `pi`.  Thus the
contact vanishes in its chosen orientation and is automatically
nonvanishing in the conjugate orientation.  This is the local content
of determinant one which an unoriented norm allocation does not record.

## 5. The five-row height budget and the remaining problem

For the full five-row profile, `A_m=4`.  With profile error, correction
height, and residue height all `o(w)`, the median-frame estimates give

```text
log|Q_1|,log|Q_2|,log|Q_3| =4w+o(w),
log||U||=4w+o(w),
log||v_4||,log||v_5||=4w+o(w),
log|Q_4|,log|Q_5|=8w+o(w).                             (19)
```

Consequently

```text
log A=8w+o(w),
log|X_i|=8w+o(w)       for i=2,3,
log|X_i|=12w+o(w)      for i=4,5,
log|s_2|,log|s_3|=o(w),
log|s_4|,log|s_5|=4w+o(w).                             (20)
```

Each of the ten retained contacts has
`log N(G_e)=2w+o(w)`.  Hence (17) has logarithmic index

```text
log[Z[i]^2:L]=20w+o(w).                               (21)
```

The four real coordinates of `(z,w)` lie in a box of logarithmic volume
at most `16w+o(w)`.  This index comparison does not give a contradiction:
a lattice of large determinant can have one exceptionally short vector,
and the equation `q=1` does not by itself control the other minima.
Indeed, take `z=1,w=i`, choose arbitrarily large pairwise coprime
Gaussian integers `G_e=r_e+i s_e` with primitive coordinate pairs, and
use `lambda_e=(r_e,s_e)`.  Such a coprime sequence can be constructed
recursively by taking `s_e=1` and making `r_e` divisible by all earlier
norms.  Then every congruence in (16) holds and `q(z,w)=1`, although the
product in (17) is arbitrarily large.  This abstract example need not
satisfy the barycentric equations or come from one moduli point; it
shows only that the covolume and `q=1` conditions cannot supply the
missing obstruction by themselves.

The new exact target is therefore narrower than the two cross-ratios:
for the prescribed rational point `(a,b)` on the fixed elliptic moduli
section, find integral scales satisfying (7), the oriented block
equations (13), and a `q=1` vector in the contact lattice (16), all with
the heights (19)--(21).  Projecting to `(A,B,C,X_i,s_i)` produces a
family of norm conics over the moduli point, but that projection forgets
the prime orientations in `bar(z)w=B+i` and (15).  No conic-bundle or
universal-torsor theorem currently supplies the required integral point;
calling the norm projection a universal torsor would therefore be
premature.

The exact identities are checked in
[check_five_point_cm_norm_lift.py](check_five_point_cm_norm_lift.py).
They identify the common-moduli lift and its height requirements.
The subsequent [contact-lattice theorem](oriented_contact_lattice_rigidity.md)
proves unrestricted integral existence for the congruences and determinant
one, even with the shared moduli retained. It does not bound the frame or
the additional row factors at the required heights.

## 6. Grouping the ten contacts through the five rows

The ten congruences in (16) are not independent random conditions.  For
an edge `ij`, write `G_ij` for (14).  The exceptional row clusters are

```text
T_14={1,4},       T_15={1,5},       T_23={1,4,5},
T_24={2,4},       T_25={2,5},       T_13={2,4,5},
T_34={3,4},       T_35={3,5},       T_12={3,4,5},
T_45={4,5}.                                             (22)
```

Consequently define

```text
J_1=G_14 G_15 G_23,
J_2=G_24 G_25 G_13,
J_3=G_34 G_35 G_12,

J_4=G_14 G_24 G_34 G_23 G_13 G_12 G_45,
J_5=G_15 G_25 G_35 G_23 G_13 G_12 G_45.                (23)
```

After the correction trimming in (11), these products divide their
corresponding rows in `Z[i]`.  The five complementary singleton cuts
give one further private oriented block in each row.  Absorbing those
private blocks and the trimmed corrections gives the exact grouped
factorization

```text
Q_i=J_i R_i,       i=1,...,5.                           (24)
```

When every original cut block has logarithmic norm `w+o(w)`, each
paired block `G_ij` has logarithmic norm `2w+o(w)`.  Hence

```text
log|J_1|=log|J_2|=log|J_3|=3w+o(w),
log|J_4|=log|J_5|=7w+o(w),
log|R_i|=w+o(w)       for every i.                      (25)
```

Thus all five residual Gaussian variables have the same remaining
height.  This is substantially more structure than the covolume count
in (17).

Write

```text
v_2=(r_2,b),       v_3=(r_3,d),
Delta=det(v_2,v_3)=r_2d-br_3.
```

The three anchor determinants `b,d,Delta` have size at most the residue
bound `T`.  The exact planar relation among the anchors becomes

```text
Delta J_1R_1-d J_2R_2+b J_3R_3=0.                     (26)
```

The subset-content identity implies

```text
gcd(N(Q_1),bd)=gcd(N(Q_2),b Delta)
 =gcd(N(Q_3),d Delta)=1.                               (27)
```

In particular the two coefficients occurring opposite each `J_a` in
(26) are units at every prime of `J_a`; no whole prime-power block is
discarded to obtain this fact.  Since `J_1,J_2,J_3` have disjoint
Gaussian support, (26) gives the three full-depth congruences

```text
dJ_2R_2 == bJ_3R_3                         (mod J_1),
Delta J_1R_1 == -bJ_3R_3                   (mod J_2),
Delta J_1R_1 == dJ_2R_2                    (mod J_3).   (28)
```

They are coupled congruences for variables of size `exp(w+o(w))`
modulo Gaussian moduli of size `exp(3w+o(w))`.  The multipliers contain
the other large `J` factors, so this size disparity alone is not a
contradiction.

Substituting (24) in the two integral-scale moduli equations (7) gives
the rest of the exact grouped system:

```text
bd J_4R_4
 =s_4((1-a)d J_2R_2+a b J_3R_3),
bd J_5R_5
 =s_5((1-b_0)d J_2R_2+b_0 b J_3R_3),                  (29)
```

where `b_0` denotes the second moduli coordinate, to distinguish it
from the anchor determinant `b`.  Rational values of `a,b_0` are
cleared exactly as in (8).  Before clearing, the row factors in every
summand have logarithmic size `8w+o(w)` when the real moduli values are
subexponential.  Clearing adds at most the logarithmic height of the
prescribed rational coefficient to every term.  Equations (26) and
(29), together with the oriented factorizations (23)--(24), retain
every row-sharing condition that the independent-CRT model forgets.

There is a useful integral norm shadow.  Since `Q_1=z`, put

```text
A=N(J_1R_1),       B=Re(bar(z)w),       AC-B^2=1,
X_i=Ar_i+Bs_i,     n_i=N(J_iR_i).
```

Then

```text
X_i^2+s_i^2=A n_i,             X_i == Bs_i (mod A),
dX_2-bX_3=A Delta,
bdX_4=s_4(dX_2-aA Delta),
bdX_5=s_5(dX_2-b_0A Delta).                       (30)
```

If `D_e=N(G_e)` and `i in T_e`, the full contact depth gives

```text
(Ar_i+Bs_i)^2+s_i^2 ==0 (mod A D_e).                   (31)
```

The weaker reduction modulo `D_e`, together with `B^2==-1 (mod A)`, is
the advertised system of simultaneous norm congruences.  It is
necessary but not sufficient: taking norms in (24) loses which member
of each conjugate prime pair occurs in `G_e`.

In the rational coordinate

```text
tau_i=d det(v_2,v_i)/(s_i Delta),       i=4,5,
```

the first nine contacts in (22) are precisely the three boundary
directions `infinity,0,1`, each repeated three times.  The last contact
is `tau_4=tau_5`.  On a fixed elliptic moduli section its intersection
with this boundary has a fixed rational common value `c_E`; after
clearing its fixed numerator and denominator, the two congruences
`tau_4==tau_5==c_E` retain the full `D_45` contact depth.  This
rephrasing does not add an inequality: coefficients such as `b,d,Delta`
must remain in the cleared equations when they are nonunits, and the
Gaussian orientations are still carried only by (24), (26), and (29).

Norm comparison in (26) is exactly balanced: every summand has
logarithmic modulus `4w+o(w)`.  The row factors in (29) are balanced at
`8w+o(w)`, before the common cost of clearing the moduli coefficients.
Thus neither the grouped norms nor the fixed boundary values presently
force an obstruction.  The concrete remaining problem is the
simultaneous solvability of (26), (29), and the oriented block
factorizations with five residual variables of height `w+o(w)`.
