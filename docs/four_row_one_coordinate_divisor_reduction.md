# Four-row small-imaginary reduction to one Gaussian divisor

The balanced four-row norm and determinant data need only **one** small
imaginary coordinate once an integral metric realization is fixed. The
remaining three smallness bounds follow from the determinant scale. Thus the
unresolved arithmetic test can be stated as a one-coordinate condition on a
Gaussian divisor of the first Gram row. This is an exact reduction, not a
height gap or an exclusion of the full profile.

## 1. Quantitative one-row propagation

Let `P_i=X_i+iY_i` be Gaussian integers, with

```text
N_i=Norm(P_i)>0,
S_ij=Re(bar(P_i)P_j),
Delta_ij=Im(bar(P_i)P_j).
```

For every `j`, the identity `P_j=(S_1j+i Delta_1j)P_1/N_1` gives

```text
Y_j=(S_1j Y_1+Delta_1j X_1)/N_1,
|Y_j| <= sqrt(N_j/N_1) |Y_1|+|Delta_1j|/sqrt(N_1).       (1)
```

The inequality uses `|S_1j|<=sqrt(N_1 N_j)` and
`|X_1|<=sqrt(N_1)`; it does not require a sign choice or a primitive
row. It applies to real metric realizations too.

In the full four-row Boolean profile, let

```text
N_i=k_i product_(T contains i)n_T,
Delta_ij=t_ij product_(T contains i,j)n_T !=0,
log n_T=w+o(w),
log k_i, log max(1,|t_ij|)=o(w).
```

Each row contains eight blocks and each pair shares four. Therefore
`sqrt(N_i)=exp(4w+o(w))`, `N_j/N_1=exp(o(w))`, and
`|Delta_1j|/sqrt(N_1)=exp(o(w))`. Formula (1) proves

```text
max_j |Y_j| <= exp(o(w)) (|Y_1|+1).                       (2)
```

In particular, `log max(1,|Y_1|)=o(w)` is equivalent to the same
condition for all four imaginary coordinates. More quantitatively, if
`|log n_T-w|<=eta w`, `log k_i<=kappa w`, and
`log max(1,|t_ij|)<=tau w`, then

```text
max_j |Y_j|
 <= exp((4eta+kappa/2)w)|Y_1|
    +exp((8eta+tau)w).                                  (3)
```

For the first term, the four blocks incident only to `j` and the four
incident only to `1` give `log(N_j/N_1)<=8eta w+kappa w`.
The second term compares four shared block upper bounds with eight
incident block lower bounds under the square root. Formula (3) is
deliberately conservative; (2) is the asymptotic statement needed here.

## 2. Exact integral divisor test

Suppose positive integers `N_i` and nonzero oriented integers
`Delta_ij` satisfy the rank-two Pluecker equations and the exact
Euclidean metric equations of
[the norm-metric realization criterion](boolean_norm_metric_realization_criterion.md).
Those equations determine integral dot products `S_ij` when the norm
data and brackets are integral. Set

```text
Z_1j=S_1j+i Delta_1j,       Z_11=N_1,
H=gcd_Z[i](Z_11,Z_12,Z_13,Z_14).
```

If `N_1` is a Gaussian norm, the same criterion proves that the set

```text
G in Z[i],     G|H,     Norm(G)=N_1                       (4)
```

is nonempty and parametrizes every integral realization of the fixed
oriented Gram matrix by

```text
P_1=bar(G),       P_j=Z_1j/G.                            (5)
```

The **exact** first-row imaginary coordinate is `Y_1=-Im(G)`.
Restrict (4) to divisors for which every row in (5) is
conjugate-primitive; call this finite set `D_prim`. Then the full
balanced primitive profile has subpower imaginary coordinates exactly
when `D_prim` contains a divisor with

```text
log max(1,|Im(G)|)=o(w).                                (6)
```

The forward implication takes `G=bar(P_1)`; the reverse implication
uses (2). Positivity and a specified order of real coordinates can
be checked on (5) if the endpoint application requires them. They
are not hidden in (6).

The divisor test is finite for each data set, but its cardinality can
grow with `w`; finiteness alone is not a height inequality. A strict
four-row obstruction would require a uniform lower bound such as

```text
min_(G in D_prim) max(1,|Im(G)|) >= exp(epsilon w)       (7)
```

for some fixed `epsilon>0`, over the admissible balanced data with
subpower `k_i,t_ij`, when `D_prim` is nonempty. Proving (7), or a weaker alternative forcing one
correction or determinant residual to have height `Omega(w)`, is the
missing arithmetic step. The present argument does not establish it.

## 3. Prescribed block orientations after primitive realization

Assume in addition that each `n_T` is an odd norm supported on split
primes, the fifteen norms have pairwise disjoint rational-prime
supports, and the norm and determinant factorizations above hold
exactly. If a realization (5) has conjugate-primitive rows, then
Gaussian factors `H_T` of norms `n_T` can be chosen so that
`H_T|P_i` for every `i in T`; the residual quotient has norm `k_i`.

Here is the primewise proof. For `p^e|n_T`, every incident row has
`p^e|N_i`. Conjugate-primitivity places all its `p`-valuation in
exactly one of the two Gaussian primes above `p`. If `i,j in T`
chose opposite orientations, then `bar(P_i)P_j` is a unit at one
Gaussian prime above `p`, while its conjugate is divisible there.
Their difference divided by `2i` is a unit there, so
`Im(bar(P_i)P_j)` is a `p`-adic unit. This contradicts
`p^e|Delta_ij`. Thus all incident rows choose one common orientation.
Collect those oriented powers over the primes of `n_T` to form `H_T`.
Singleton blocks need no pair comparison. Disjoint support permits
all fifteen divisions simultaneously, and taking norms leaves
exactly `k_i` in row `i`.

This orientation observation is implicit in
[primitive completion](boolean_metric_primitive_completion.md), which
also handles nonprimitive initial realizations with subpower losses.
It explains the scope of (7): after the metric equations and primitive
integral lift, the four simultaneous angular inequalities reduce to a
single short-divisor inequality. The real algebraic polynomial profile
does not settle this divisor question because it does not give an
integral Gram realization with rational Gaussian factors.

## 4. The endpoint chamber gives a nearest-square congruence test

In the positive-real endpoint chamber, the first-row divisor search
has at most two candidates at subpower imaginary height. Put

```text
a=floor(sqrt(N_1)),       r_0=N_1-a^2,       0<=r_0<2a+1.
```

If `P_1=X_1+iY_1`, `X_1>0`, and `Y_1^2<2a-1`, then `X_1=a`.
Indeed `X_1<=a-1` would give
`Y_1^2=N_1-X_1^2>=N_1-(a-1)^2>=2a-1`.
For the balanced four-row profile, `a=exp(4w+o(w))`, so every
subpower `Y_1` satisfies this condition for large `w`.

Consequently such a realization requires `r_0=b^2` for a nonnegative
integer `b=exp(o(w))`, and

```text
P_1=a+i sigma b,       G=bar(P_1)=a-i sigma b,
sigma in {+1,-1}.                                      (8)
```

For each sign, the remaining rows are forced by (5). Their integrality
is exactly the three Gaussian congruences

```text
a-i sigma b divides S_1j+i Delta_1j,       j=2,3,4.  (9)
```

Equivalently, both of the following ordinary integer congruences must
hold for each `j`:

```text
a S_1j-sigma b Delta_1j = 0 mod N_1,
sigma b S_1j+a Delta_1j = 0 mod N_1.                  (10)
```

If `gcd(a,b)=1`, then `gcd(a,N_1)=1`, and the first congruence in
(10) implies the second. This is the case for a conjugate-primitive
candidate of odd norm. After (9), test conjugate-primitivity of the
three forced rows; Section 3 then supplies all fifteen oriented
blocks. Formula (1) supplies their subpower imaginary bounds.

Thus the exact missing arithmetic inequality can be posed without a
search over arbitrary Gaussian rotations: for every admissible
balanced four-row metric data set, prove that either the nearest-square
remainder `r_0` is not a subpower square, or both sign choices in (8)
fail (9) or primitivity. A fixed positive height gap would need this
failure uniformly whenever `b<exp(epsilon w)` for some fixed
`epsilon>0`. Neither the metric equations nor the real polynomial
profile proves such a gap.

The [integer-star checker](check_four_row_integer_five_star_lift.py)
independently verifies the threshold, reconstruction identities and both
congruence signs on the existing three-row family. The argument above is
valid in any row count with the stated metric scales; the test does not
assert existence of a full four-row small-residue profile.

## 5. An exact private-norm gcd and a unique candidate sign

The reduced pair quotients give an additional necessary arithmetic
test on the first-row candidate. This test uses actual row primitivity
and retains all correction-factor overlaps.

Let `I={1,2,3,4}`, `k_i=Norm(K_i)`, and now write

```text
Z_ij=bar(P_i)P_j/G_ij=x_ij+i t_ij,
e_i=n_(I\{i}),       h_i=k_i n_{ {i} },
A_i=h_i e_i.
```

Here `A_i` is a grouped private norm, not the nearest-square integer
`a=floor(sqrt(N_1))` in Section 4. The triangle identities in
[the frame note](endpoint_triangle_frame_holonomy.md#8-seven-norm-variables-and-an-exact-additive-gcd)
recover it by the exact additive gcd

```text
A_i=gcd_({j,k} subset I\{i})
       |(t_ik x_ij-t_ij x_ik)/t_jk|.                   (11)
```

All `t_jk` are assumed nonzero. Each quotient is an integer, and
the three quotients are `A_i` times the three pairwise coprime
matching norm products.

For every `j!=i`, the stronger actual-coordinate identity is

```text
gcd(A_i, X_i t_ij+Y_i x_ij)=h_i,
e_i=A_i/gcd(A_i,X_i t_ij+Y_i x_ij).                    (12)
```

The gcd is the positive ordinary integer gcd, also when its second
argument is negative. To prove (12), multiply the reduced quotient
by `P_i` and take imaginary parts:

```text
X_i t_ij+Y_i x_ij=(N_i/G_ij)Y_j=h_i c_ij Y_j,
c_ij=product_(T contains i but not j, T!={i}) n_T.     (13)
```

There are three core norms in `c_ij`. None is `e_i`, so disjoint
core support gives `gcd(c_ij,e_i)=1`. Also `e_i|N_j`; since `P_j`
has coprime integer coordinates, `gcd(e_i,Y_j)=1`. Consequently

```text
gcd(h_i e_i,h_i c_ij Y_j)=h_i.
```

No coprimality between `h_i` and `e_i` is needed: `K_i` may contain
arbitrary powers of primes dividing `e_i`. Thus (12) includes their
full valuation depth.

There is a useful finite consequence for the two signs in (8).
Fix `X_i>0`, `b>0` with `X_i^2+b^2=N_i` and `gcd(X_i,b)=1`.
If both candidates `Y_i=b` and `Y_i=-b` pass the necessary test
on the same edge `ij`, then

```text
gcd(A_i,X_i t_ij+b x_ij)
 =gcd(A_i,X_i t_ij-b x_ij)=h_i
   implies n_{ {i} } | t_ij.                          (14)
```

Indeed `n_{ {i} }|h_i` then divides both expressions, hence divides
`2X_i t_ij`. This core norm is odd, and
`gcd(X_i,n_{ {i} })=1`, since `n_{ {i} }|N_i` and the candidate
is primitive. Therefore it divides `t_ij` itself, including every
prime-power depth. In particular, if

```text
0<|t_ij|<n_{ {i} },
```

at most one sign can pass (12). This inequality holds eventually
in the balanced subpower-residual setting. One can therefore apply
(12) as a sign filter before the full Gaussian divisibilities (9).
Passing it does not establish those divisibilities or the remaining
row primitivity.

Equations (11)--(14) are additional exact tests on the full profile.
They do not prove a positive lower bound on `b` or on the residuals.
The [integer-star checker](check_four_row_integer_five_star_lift.py)
verifies the additive and coordinate gcds, exact co-singleton
recovery, and the sign implication on full fifteen-block fixtures
with correction factors overlapping the co-singleton core at
several valuation depths. These finite fixtures are not asserted
to satisfy endpoint smallness.

## 6. The three gcd outputs differ only at residual scale

The tests (12) must not be counted as three independent large-modulus
constraints. Fix a row `i` and arbitrary candidate integers `X,Y`, with
no norm or primitivity assumption on this candidate, and set

```text
A=A_i,       L_j=X t_ij+Y x_ij,
d_j=gcd(A,L_j),       d=gcd(A,(L_j)_(j!=i)).
```

All gcds are positive. The triangle identities behind (11) imply
`A | t_ik x_ij-t_ij x_ik`, and hence

```text
A | t_ik L_j-t_ij L_k.
```

It follows that `d_j | t_ij L_k` for every `k`. Since also
`d_j | A`, Bezout applied to `A` and all the `L_k` gives
`d_j | t_ij d`. As `d | d_j`, this proves

```text
d_j/d divides |t_ij|.                                (15)
```

In particular, the three outputs have the same logarithmic size up to
`log max_j |t_ij|`. Their valuations agree at every prime not dividing
the corresponding residual. They need not be identical at residual
primes. No coprimality between corrections and core blocks was used.
This limits independent divisibility counting; it neither removes the
sign obstruction nor proves a height gap.

There is also a depth-safe version of the root-of-minus-one description.
For one edge put `u=gcd(x_ij,t_ij)` and

```text
A'=A/gcd(A,u^2).
```

The norm factorization gives `A | x_ij^2+t_ij^2`, so primewise
cancellation proves `A' | (x_ij/u)^2+(t_ij/u)^2`. The two reduced
coordinates are coprime, and therefore `t_ij/u` is a unit modulo `A'`.
For `A'>1`,

```text
r=(x_ij/u)(t_ij/u)^(-1) mod A',       r^2=-1 mod A',
A' >= A/|t_ij|^2.                                   (16)
```

The case `A'=1` is vacuous. This operation removes at most the
valuation depth present in `u^2`; it does not discard an entire large
core prime power merely because the residual contains that prime.
The checker now verifies (15) on arbitrary candidate coordinates and
(16) on the same full-profile fixtures with overlapping corrections.

## 7. Exact imaginary-content filter, and its normalization limit

There is a small but depth-exact constraint on the actual pair residuals.
For any two conjugate-primitive rows in the full profile, with nonzero
imaginary coordinates and `Delta_ij=G_ij t_ij`, one has

```text
gcd(|Y_i|,|t_ij|)=gcd(|Y_i|,|Y_j|)
                   =gcd(|Y_j|,|t_ij|).                 (17)
```

This holds even when the correcting factors overlap core primes. Indeed,
integer-coordinate primitivity gives `gcd(X_i,Y_i)=1`, hence
`gcd(Y_i,N_i)=1`. Since `G_ij|N_i,N_j`, the integer
`d=gcd(|Y_i|,|Y_j|)` is coprime to `G_ij`. It divides
`Delta_ij=X_iY_j-X_jY_i=G_ij t_ij`, so `d|t_ij`.
Conversely, if `d'|Y_i,t_ij`, then `d'|Delta_ij`, whence
`d'|X_iY_j`. Coprimality of `X_i,Y_i` forces `d'|Y_j`.
Interchanging the rows proves (17).

At any rational prime `p`, (17) says

```text
min(v_p(Y_i),v_p(t_ij))=min(v_p(Y_i),v_p(Y_j)).       (18)
```

In particular, if the two `Y`-valuations differ, then
`v_p(t_ij)` equals their minimum; extra residual depth is possible
only when those valuations agree. Equivalently,
`gcd(|Y_i|,(|t_ij|)_(j!=i))=gcd_i |Y_i|`.

Put `g=gcd_i |Y_i|`, `y_i=Y_i/g`, and `r_ij=t_ij/g`.
Equation (17) gives `g|t_ij` for every pair and `gcd(g,G_ij)=1`.
The virtual integer vectors `X_i+i y_i` still satisfy

```text
X_i y_j-X_j y_i=G_ij r_ij,
gcd(|y_i|,|r_ij|)=gcd(|y_i|,|y_j|).                   (19)
```

Every triangle circuit in Section 1 of
[the Boolean circuit note](full_boolean_small_imaginary_circuits.md)
loses exactly the common factor `g^2` when written in `y_i,r_ij`.
Its large private norm factors and remaining coefficient ratios do
not improve. More decisively, the virtual vectors need not have the
original row norms or Gaussian block divisibilities:

```text
Norm(X_i+i y_i)=N_i-(g^2-1)y_i^2.                    (20)
```

For example, take `H=2+i`, `P_1=H(-1+2i)=-4+3i`, and
`P_2=H(5+2i)=8+9i`. They are conjugate-primitive, share the
oriented factor `H` of norm `5`, and have `g=3`.
Their determinant is `-60=5*(-12)`; the virtual determinant is still
`-20`, divisible by `5`, but `Norm(-4+i)=17` is not. Thus content
normalization preserves the determinant congruences, not the full
Gaussian profile. Equation (17) is a local filter, not a positive
residual-height gap.
