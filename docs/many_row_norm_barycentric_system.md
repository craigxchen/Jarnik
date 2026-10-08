# Full-subset norm interpolation and two limits of a many-row descent

The actual Gaussian norm equations give a family of barycentric moments
involving the singleton blocks, which do not occur in the ordinary
three- and four-row determinant circuits. This note gives the exact
normalization for arbitrary nonzero, nonconstant small imaginary
coordinates. It also gives an unbounded, prime-disjoint countermodel to
the isolated four-row norm circuit, preserving the coupling of its six
edge residuals and conjugate-primitivity of its correcting factors.
Finally, the five-row Gaussian circuit loses two blocks of rational
content when made primitive. These are limits of particular proposed
steps, not an exclusion or construction of the full endpoint profile.

## 1. Norm moments with every subset retained

Fix `m>=4` and suppose

```text
P_i=X_i+iY_i=K_i product_(T containing i) H_T,
n_T=Norm(H_T),        k_i=Norm(K_i),
Delta_ij=X_i Y_j-X_j Y_i
        =t_ij product_(T containing i,j) n_T != 0.
```

Extend `t_ji=-t_ij`. Assume `Y_i!=0`; no equality among the
`Y_i` is needed. Set

```text
tau_i=product_(j!=i) t_ij,
A_i=n_{ {i} } product_(T not containing i, |T|>=3) n_T^(|T|-2),
C_i=k_i A_i/tau_i,
L=product_(|T|>=3) n_T^(|T|-2),
Y=product_i Y_i,       x_i=X_i/Y_i.
```

The product defining `A_i` includes every eligible subset of the other
`m-1` rows. In particular, it has no full-set factor. Direct cancellation
of each incident block gives

```text
(x_i^2+1)/product_(j!=i)(x_i-x_j)
                  = (Y/L) C_i Y_i^(m-4).                 (1)
```

Indeed, an incident block occurs once in `Norm(P_i)` and `|T|-1`
times in the product of incident determinants. Its surviving exponent
is `2-|T|`; multiplying by `L` leaves exactly the stated `A_i`.
This calculation uses the actual norm equation, not only the rank-two
minor identities.

The top `m-3` coefficients in Lagrange interpolation of `x^2+1`
therefore yield

```text
sum_i C_i X_i^j Y_i^(m-4-j)=0,       0<=j<=m-4.           (2)
```

All exponents in (2) are nonnegative. Multiplication by
`product_(a<b)|t_ab|` clears every denominator. Under the full Boolean
profile, `log|t_ab|=o(w)` and `log k_i=o(w)`, so this is an integral
moment system with only subpower clearing cost. It permits arbitrary
subpower, nonconstant `Y_i`. The norm-block part of each `A_i` has height

```text
log A_i = ((m-5)2^(m-2)+m+2)w+o(w).                       (3)
```

For four rows the sole lower moment is the new norm circuit

```text
sum_i (k_i/tau_i) n_{ {i} } n_[4]\{i}=0.                (4)
```

After replacing individual rows by their negatives if necessary, take
`Y_i>0`, and order the distinct real numbers `x_1<...<x_m`.
The positive norms in (1) imply

```text
sign(C_i)=(-1)^(m-i).
```

Thus (4) is an equality between two positive sums. Positivity fixes
these alternating signs; it does not make all the terms have one sign.

## 2. The full positive quadratic, including its discriminant

Assume now the original rows are conjugate-primitive, as supplied by
the central endpoint reduction. The common Gaussian block
`H_[m]=A+iB` is then conjugate-primitive and has odd norm `g`.
In particular `gcd(B,g)=1`. Choose an integer `r` with
`r=A B^(-1) mod g`, and put

```text
a_i=(X_i-rY_i)/g in Z,       s=(r^2+1)/g in Z,
b_i=a_i/Y_i,                f(t)=g t^2+2rt+s.
```

These statements follow from divisibility by the actual orientation
`A+iB`; they are not consequences of norm divisibility alone. They give

```text
Norm(P_i)=g f(a_i,Y_i),
f(a,y)=g a^2+2ray+s y^2,          gs-r^2=1.              (5)
```

Let `L_0` be `L` with the full-set factor `g^(m-2)` removed, and put

```text
M_j=sum_i C_i Y_i^(m-4) b_i^j,
e_1=sum_i b_i,        e_2=sum_(i<j) b_i b_j,
alpha=Y/L_0.
```

The full barycentric identity is

```text
f(t)=alpha sum_i C_i Y_i^(m-4) product_(j!=i)(t-b_j).
```

In addition to `M_0=...=M_(m-4)=0`, its remaining coefficients are

```text
g  =alpha M_(m-3),
2r =alpha (M_(m-2)-e_1 M_(m-3)),
s  =alpha (M_(m-1)-e_1 M_(m-2)+e_2 M_(m-3)).            (6)
```

Write the three parenthesized quantities, in order, as `U,V,W`.
The exact positive quadratic condition is

```text
alpha U>0,                alpha^2(4UW-V^2)=4.           (7)
```

Conversely, the zero moments and (7) make the displayed interpolant
a strictly positive quadratic of discriminant `-4`; (6) supplies its
coefficients. Their required integrality and the oriented Gaussian
divisibilities remain additional arithmetic conditions in a converse
construction. Equations (6)--(7) do not assert them automatically.

For `m=4,Y_i=1`, the formulas are especially short:

```text
M_0=0,
gL_0=M_1,
2rL_0=M_2-e_1M_1,
sL_0=M_3-e_1M_2+e_2M_1,
4M_1(M_3-e_1M_2+e_2M_1)-(M_2-e_1M_1)^2=4L_0^2.
```

Here `L_0` is the product of all four triple norms. Thus the overlapping
barycentric denominators are explicit. Replacing these rational
quotients by independent integers would discard four full blocks of
denominator height.

## 3. A balanced arithmetic countermodel to the isolated norm circuit

Define eight positive quadratics `Q_j(t)=(t+a_j)^2+b_j^2`, where

```text
(a_j,b_j)=(-5,2),(-5,3),(-2,5),(3,3),
          (-5,1),(-4,4),(-1,5),(1,5).
```

Exact polynomial multiplication gives

```text
Q_1 Q_2+Q_3 Q_4=Q_5 Q_6+Q_7 Q_8.                       (8)
```

For example, the four quartics in this order have coefficient tuples
in descending degree

```text
(1,-20,163,-630,986),     (1,2,23,102,522),
(1,-18,138,-528,832),     (1,0,48,0,676).
```

They use eight different norm quadratics. Put

```text
S={2,3,5,7,13,17,29,37,53,73,89,101,109},
M=348442128341958211928640,
t=M h,           h=1,2,...,
J_j(t)=(t+a_j+i b_j)/(a_j+i b_j),
v_j(t)=Norm(J_j(t))=Q_j(t)/Q_j(0).
```

Then all `J_j(t)` are Gaussian integers; the eight `v_j(t)` are odd,
conjugate-primitive Gaussian norms with pairwise disjoint rational-prime
support, and

```text
log v_j(t)=2 log t+O(1).
```

Here is a proof covering the whole progression. The set `S` contains
every prime dividing an imaginary constant `b_j`, a constant `Q_j(0)`,
or a resultant of two distinct `Q_j`. All those resultants are
nonzero. For each `p in S`, `M` contains
`p^(1+max_j v_p(Q_j(0)))`. Therefore

```text
v_p(Q_j(t))=v_p(Q_j(0)),        p in S,
```

and `v_j(t)` has no factor from `S`. Also `Q_j(0)` divides `M`,
which proves that `J_j(t)` is Gaussian integral. A common prime of
two `v_j(t)` would divide their resultant and hence lie in `S`,
impossible. An inert prime dividing a `Q_j(t)` must divide `b_j`,
and is likewise in `S`. A rational prime dividing both coordinates
of `J_j(t)` would either lie in `S`, or, after multiplication by
`a_j+i b_j`, divide both coordinates of `t+a_j+i b_j`; the latter
also forces it into `S`. Finally `2` is absent. Thus the claims include
exact conjugate-primitivity, not just absence of inert primes.

Equation (8) becomes

```text
986 v_1v_2+522 v_3v_4=832 v_5v_6+676 v_7v_8.             (9)
```

Assign the singleton/complementary-triple pairs of rows `1,2,3,4`
respectively to `(v_1,v_2),(v_5,v_6),(v_3,v_4),(v_7,v_8)`.
The edge-product coupling in (4) can be preserved exactly. Take

```text
(K_1,K_2,K_3,K_4)=(10+11i,1,1,1),
|t_12,t_13,t_14,t_23,t_24,t_34|=(39,1352,12,1,87,4),
t_ij<0 for i<j,      t_ji=-t_ij,      D=2822976.
```

Every `K_i` is conjugate-primitive. The absolute incident products are
`(632736,3393,5408,4176)`, so

```text
(k_i/tau_i)_i=(-986,832,-522,676)/D.
```

Consequently (4) holds exactly with the required alternating signs,
bounded coupled residuals, primitive correcting factors, and eight
prime-disjoint oriented Gaussian blocks of equal logarithmic height.
One may add any finite number of distinct primitive Gaussian linear
polynomials and enlarge `S,M` by the same construction to supply the
remaining seven block positions at the same height.

This is a countermodel to obtaining a strict block-height gain from
the isolated norm circuit plus those listed conditions. It does not
realize the pair determinants, their Pluecker equations, small row
imaginary coordinates, or the positive quadratic constraints (6)--(7).
In particular it is not an endpoint counterexample. The same caution
applies to interpreting actual orientations of the individual blocks
as compatible orientations of the required row products.

## 4. The five-row Gaussian relation has large rational content

For `m=5`, (2) gives both `sum C_iY_i=0` and `sum C_iX_i=0`.
Combining them and removing the common Gaussian factor gives

```text
sum_i (k_i K_i/tau_i) U_i=0,
U_i=A_i product_(T containing i)H_T / product_(|T|>=3)H_T.
```

Prime by prime, this is the exact factorization

```text
U_i=H_i^2 bar(H_i) H_[5]\{i} bar(H_[5]\{i})^2
       product_(j!=i) H_{ij} bar(H_[5]\{i,j})
   =n_i n_[5]\{i} V_i,
V_i=(H_i bar(H_[5]\{i}))
       product_(j!=i)(H_{ij} bar(H_[5]\{i,j})).           (10)
```

Under the original prime-disjoint, conjugate-primitive hypotheses,
`V_i` is conjugate-primitive. The ordinary integer coordinate content
of `U_i` is therefore exactly `n_i n_[5]\{i}`. Its logarithm is
`2w+o(w)`, whereas

```text
log|U_i|=7w+o(w),             log|V_i|=5w+o(w).
```

The new primitive pair blocks
`H_{ij} bar(H_[5]\{i,j})` have norm height `2w+o(w)` and occur
in rows `i,j`. The private block `H_i bar(H_[5]\{i})` has the
same norm height. This does reduce the support pattern to singletons
and pairs, but its primitive relation has coefficients

```text
(k_i K_i/tau_i) n_i n_[5]\{i},
```

of height `2w+o(w)`, not subpower height. It is important here to
retain the angular information correctly. Set `Z_i=K_i V_i`. Then

```text
Z_i=(A_i/(n_i n_[5]\{i})) P_i / product_(|T|>=3)H_T,
Z_i/bar(Z_i)=(P_i/bar(P_i))
             product_(|T|>=3)bar(H_T)/H_T.              (11)
```

The first prefactor is positive real. Thus the `Z_i` have exactly the
same angular differences as the original rows, with one common rotation.
Any rational content created by `K_i V_i` is bounded by `Norm(K_i)`
and is subpower. Small imaginary coordinates relative to the original
integer real axis need not survive this common rotation, but angular
clustering does survive.

Moreover, the fifteen factors

```text
J_T=H_T bar(H_[5]\T),           |T|=1 or 2,
```

are exactly the fifteen nontrivial unoriented cuts on five labels,
each of norm height `2w+o(w)`. Their incidence in `V_i` chooses the
side `T` containing `i`. This is the familiar complement-cut grouping
obtained by forgetting the original anchor; (11) represents the same
five directions. Reanchoring this system yields the usual four-row
full Boolean normalization at block scale `2w`. Eight of these fifteen
cuts separate any fixed pair of labels, so each reanchored numerator
has logarithmic modulus `8w+o(w)`. This equals the original five-row
numerator modulus, from sixteen incident norm blocks of height `w`.
The circuit alone
has not improved the inherited angular precision or supplied a new
strict height inequality. It cannot be regarded as a new height
descent merely because its written support pattern has changed.

This is a different normalization from the fixed-polynomial character
tests in [the additive generic-content note](additive_balanced_character_generic_content.md).
Their general caution about rational content agrees with (10), while
the present identity is obtained from actual norm interpolation.

The [exact checker](check_many_row_norm_barycentric_system.py) verifies
the coefficient formulas with nonconstant imaginary coordinates, the
quartic identity, the resultant-prime progression certificate, its
primitive Gaussian factors and edge residual coupling, and all Boolean
exponents in the five-row factorization. The general proofs above are
independent of the sampled values.
