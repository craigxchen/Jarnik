# A matching-norm gap for four thin signed witnesses

Four signed products in the pattern below satisfy the unconditional bound
`min(N(A),N(B),N(C))<=8 M^2 max(1,Y)^4`, where `M` is the largest
correction modulus and `Y` is the largest imaginary coordinate.
In particular, balanced growing groups cannot all have subpower
imaginary coordinates after subpower Gaussian corrections.
The obstruction uses the repeated norm divisors on opposite edges of the
four witnesses. Their simultaneous three-row determinant identities force
two independent integer relations on three coprime norms.

This concerns actual thin witnesses, not merely an algebraic relation
between four signed products. In particular, it does not contradict the
bounded-coefficient family in
[the four-term circuit note](four_term_signed_circuit_height.md).
It also does not show that a general endpoint profile supplies these four
witnesses on one common support.

## 1. An exact lemma for repeated opposite-edge divisors

Let `U_j=X_j+iY_j` be four Gaussian integers, indexed by `0,1,2,3`,
with every `Y_j` nonzero and every pair determinant

```text
Delta_ij=X_i Y_j-X_j Y_i,       i<j,
```

nonzero. Suppose positive integers `s,nA,nB,nC` satisfy
`gcd(nA,nB,nC)=1`, and the six determinants have the exact form

```text
Delta_01=s nA t_01,       Delta_23=s nA t_23,
Delta_02=s nB t_02,       Delta_13=s nB t_13,
Delta_03=s nC t_03,       Delta_12=s nC t_12,
```

where all six `t_ij` are nonzero integers. Put

```text
Y=max_j |Y_j|,       T=max_(i<j) |t_ij|.
```

Then

```text
max(nA,nB,nC) <= 2 Y^2 T^2.                         (1)
```

**Proof.** The four three-row determinant identities, divided by `s`,
give an integer matrix annihilating `(nA,nB,nC)^t`:

```text
  [ Y_2 t_01   -Y_1 t_02    Y_0 t_12 ]
  [ Y_3 t_01    Y_0 t_13   -Y_1 t_03 ]
  [ Y_0 t_23    Y_3 t_02   -Y_2 t_03 ]
  [ Y_1 t_23   -Y_2 t_13    Y_3 t_12 ].             (2)
```

This matrix has rank exactly two over the rationals. For the lower
bound on its rank, divide each row by the product of the three `Y`'s
whose labels occur in that triangle, and set

```text
a=t_01/(Y_0Y_1), b=t_02/(Y_0Y_2), c=t_03/(Y_0Y_3),
d=t_12/(Y_1Y_2), e=t_13/(Y_1Y_3), f=t_23/(Y_2Y_3).
```

The resulting rows are

```text
(a,-b,d),       (a,e,-c),       (f,b,-c),       (f,-e,d).
```

All six letters are nonzero. If the first two rows are proportional,
their equal nonzero first coordinates make them equal, so `e=-b` and
`c=-d`. If the third row is also proportional to the first, its second
coordinate forces the proportionality factor to be `-1`. Its third
coordinate is then simultaneously `d` and `-d`, a contradiction.
The nonzero null vector already gives the upper rank bound of two.

Choose two independent rows of (2). Their nonzero cross product is an
integer vector parallel to `(nA,nB,nC)`. Since the latter vector has
coordinate gcd one, the cross product is a nonzero integer multiple of
it. Every matrix entry has absolute value at most `YT`; hence every
coordinate of the cross product has absolute value at most `2Y^2T^2`.
This proves (1). Notice that pairwise coprimality of the three norms is
stronger than needed for this lemma.

The same argument also retains the three different column heights.
Define

```text
T_A=max(|t_01|,|t_23|), T_B=max(|t_02|,|t_13|),
T_C=max(|t_03|,|t_12|).
```

The three columns of (2) have heights at most `YT_A,YT_B,YT_C`.
The first coordinate of any nonzero cross product is a nonzero integer
multiple of `nA`. Consequently

```text
nA <= 2Y^2 T_B T_C,                                (1a)
```

with the analogous two bounds obtained by permuting the letters.

## 2. Four corrected signed products

Let `A,B,C` be nonunit conjugate-primitive Gaussian integers, and let
`D` be a nonzero conjugate-primitive Gaussian integer, possibly a unit.
Assume their norms are pairwise coprime. Write

```text
nA=N(A), nB=N(B), nC=N(C), nD=N(D),
P=nA nB nC nD,
B_0=DABC,              B_1=DA bar(B) bar(C),
B_2=D bar(A) B bar(C), B_3=D bar(A) bar(B) C.
```

For arbitrary nonzero Gaussian integers `K_j`, put

```text
U_j=K_j B_j,
M=max_j |K_j|,                 Y=max_j |Im U_j|,
Q_*=min(sqrt(nA nB),sqrt(nA nC),sqrt(nB nC)).
```

The corrections may overlap core primes, and the `U_j` need not be
conjugate-primitive. There is the following exact alternative:

```text
M >= Q_*
```

or

```text
nD <= 8 M^2 Y^4.                                   (3)
```

**Nonvanishing when `M<Q_*`.** Every core `B_j` is
conjugate-primitive and has modulus `sqrt(P)`. If `U_j` were real,
conjugating it and cancelling the coprime factors `B_j,bar(B_j)`
would give `bar(B_j)|K_j`. Thus `|K_j|>=sqrt(P)>=Q_*`, impossible.
Consequently every imaginary coordinate is nonzero.

For two different core signs, write their products as `G E` and
`G bar(E)`, where the norm of `E` is one of
`nA nB,nA nC,nB nC`. If the corrected witnesses were collinear,
the equality of a product and its conjugate would give

```text
E^2 K_i bar(K_j)=bar(E)^2 bar(K_i) K_j.
```

Since `E` and `bar(E)` are coprime,
`E^2 | bar(K_i)K_j`. Taking moduli yields
`N(E)<=|K_i||K_j|<=M^2`, again contradicting `M<Q_*`.
Thus all six determinants are nonzero.

**The repeated norm divisors.** The common core factors of pairs
`01,23` have norm `nD nA`; those of `02,13` have norm `nD nB`;
those of `03,12` have norm `nD nC`. For example, the common factors
of `B_0,B_1` and `B_2,B_3` are respectively `DA` and `D bar(A)`.
If a Gaussian integer divides two witnesses, its norm divides their
integer determinant. Therefore the six determinant factorizations in
Section 1 hold with `s=nD`. Extra common factors from the corrections
do not affect this divisibility or require deleting any core primes.

Since `|U_j|<=M sqrt(P)`, the resulting nonzero integer residuals obey

```text
T_A <= 2 M sqrt(P) Y/(nD nA),
T_B <= 2 M sqrt(P) Y/(nD nB),
T_C <= 2 M sqrt(P) Y/(nD nC).                       (4)
```

Apply the column-specific bound (1a), using the coprimality of the
core norms, and then (4):

```text
nA <= 2Y^2 T_B T_C <= 8 P M^2 Y^4/(nD^2 nB nC)
                   = 8 M^2 Y^4 nA/nD.
```

Cancel `nA` to prove (3). More generally, this calculation proves (3)
whenever all imaginary coordinates and all determinants are nonzero;
the cutoff `M<Q_*` was only used to guarantee those conditions.
The calculation also permits unit groups under the same nonvanishing
conditions.

## 3. A bound without balanced group sizes

There is also a uniform consequence with no balance assumption. In the
notation of Section 2, put

```text
Yhat=max(1,Y),       H=max(M,Yhat).
```

Then

```text
min(nA,nB,nC) <= 8 M^2 Yhat^4 <= 8 H^6.             (5)
```

This bounds the smallest of the three changing groups. It is not a
bound for their maximum in all configurations.

**Proof.** Put `a=min(nA,nB,nC)`. If `a<=M^2`, (5) follows
immediately from `Yhat>=1`. Otherwise every one of the six products

```text
nA nB, nA nC, nB nC, nD nA, nD nB, nD nC
```

exceeds `M^2`, because `nD>=1`. Also
`sqrt(P)>=a^(3/2)>M`, so the argument in Section 2 excludes zero
imaginary coordinates.

One may make any of the four groups the common positive group by
conjugating the witnesses where that group has negative orientation
and permuting the group and witness labels. For `A,B,C`, respectively,
the conjugated witness pairs are `23,13,12`. This preserves `M` and
all absolute imaginary coordinates and leaves the same four-sign
pattern. Between any two different witnesses of a reoriented quartet,
the flip norm is one of the six products displayed above. The
collinearity calculation in Section 2 would require that norm to be
at most `M^2`. Thus all determinants remain nonzero after each of
these reorientations.

The common-group bound (3) therefore applies with any group chosen
as common. In this branch it gives the stronger conclusion

```text
max(nD,nA,nB,nC) <= 8 M^2 Y^4,
                     if min(nA,nB,nC)>M^2.          (6)
```

This proves (5). The other branch only bounded the smallest changing
group; no global maximum-norm assertion follows. If `D` is a unit,
it can become a changing group during reorientation: the nonvanishing
conditions just proved suffice for the common-group calculation.

## 4. A positive exponent in the balanced case

Suppose, along a sequence with `w` tending to infinity,

```text
log nA=log nB=log nC=log nD=w+o(w).
```

Then `log Q_*=w+o(w)`. If the corrections are subpower,
`log M=o(w)`, the first alternative in Section 2 is unavailable for
large `w`, and (3) gives

```text
2 log M+4 log Y >= w+o(w),
Y >= exp(w/4-o(w)).                                (7)
```

In particular, all four imaginary coordinates cannot be subpower.
Without assuming subpower corrections, a common upper bound `H>=1`
for both `M` and `Y` must satisfy

```text
H >= exp(w/6-o(w)).                                (8)
```

Indeed the first alternative gives the stronger `H>=Q_*`, while
the second gives `H^6>=exp(w+o(w))`.

The arithmetic step here is the rank-two system (2) on the same three
coprime norms. One isolated four-term signed relation does not provide
that system. Nor does conjugating one relation give a second
independent relation automatically. The circuit examples in the
four-term note therefore remain consistent with (7).

## 5. The limitation for general endpoint profiles

The theorem requires a common signed support, with the *same* three
large norm divisors repeated on opposite edges. The full four-row
Boolean profile and its five-star lift have different support data.
In particular, the star lift has six distinct, pairwise-coprime edge
norms on any four of its stars. Those six edge norms are not the three
repeated quantities in Section 1.

The available cheap imaginary-coordinate identities on original
endpoint rows also have differing supports. The common-monomial cost
in [the four-term circuit note](four_term_signed_circuit_height.md)
prevents identifying them with four subpower-corrected witnesses here
by that restricted normalization. A new construction or arithmetic
reduction is needed to make (5) or (7) apply to a general endpoint cluster.
No radius-uniform lattice-point count follows from this theorem alone.

There is a related distinction for the genuine thin pair products
`Z_ij=P_i bar(P_j)/G_ij`. Their support indicators obey Boolean
parallelogram identities, but their Gaussian orientation data have
three states: positive, negative, and absent. Multiplying pair products
adds signed valuations instead of taking their sum modulo two.
On a subset containing two positive endpoints and neither negative
endpoint, for example, the product contains `H_T^2`. Dividing out
one Gaussian block changes the phase and has no proved thinness
bound. The support identity alone is therefore insufficient to
construct the four witnesses required here.

The [exact checker](check_four_signed_witness_matching_norm_gap.py)
checks the repeated-divisor matrix, its primitive kernel, signed-product
divisibilities and the quantitative bounds. These finite checks supplement
the proof for arbitrary Gaussian integers given above; they do not certify
the missing reduction from general endpoint profiles.
