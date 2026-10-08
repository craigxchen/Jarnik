# Second-order vertical phases: a sharp conditional gain and no forced relation

Let

```text
w_j=-1+i t_j,
```

and measure a signed monomial by the ten-block cost

```text
h(r) = (sum_(i<j) |r_i+r_j|
          + sum_(i<j<k) |r_i+r_j+r_k|)/6.              (1)
```

This note separates two facts.  A short three-term arithmetic progression
of heights gives a genuine second-order phase estimate, including exact
control of every Gaussian gcd and ordinary-content loss.  But four actual
fixed-real-part cofactors do not force any reciprocal cancellation of cost
less than three.  The latter remains false even for full primitive circle
tuples with four equal chord contents and mutually comparable heights.

Thus second-order phases can improve the first-order monomial barrier on
structured height configurations.  No exclusion of the formal five-thirds
profile follows from the arguments here: compatibility between that sharing
profile and the endpoint hypotheses remains open.

## 1. The useful second difference

Put

```text
(t_1,t_2,t_3)=(t,t+q,t+2q),       r=(1,-2,1,0).
```

Then `h(r)=2`, and the exact reciprocal identity is

```text
1/t - 2/(t+q) + 1/(t+2q)
  = 2q^2/[t(t+q)(t+2q)].                              (2)
```

Consequently

```text
eta = atan(1/t)-2 atan(1/(t+q))+atan(1/(t+2q))
    = O(q^2/t^3)                                      (3)
```

when `q=O(t)`.  In particular this is a cubic phase when `q=O(1)`.
The same assertion with constants depending on a fixed real part `d`
holds for `w_j=-d+i t_j`.

There is no hidden content loss in (3).  For general `d`, set

```text
U=w(t)w(t+2q),        V=w(t+q)^2,       g=gcd_G(U,V).
```

Direct multiplication gives

```text
Im(U conjugate(V)) = 2d(t+q)q^2.                       (4)
```

Therefore

```text
Z=(U/g) conjugate(V/g)
```

has nonzero integral imaginary part, and

```text
N(g) | 2d(t+q)q^2.                                    (5)
```

If the ordinary content `c=gcd(Re Z,Im Z)` is removed as well, then

```text
c N(g) | 2d(t+q)q^2.                                  (6)
```

The primitive value has the same phase.  Hence Gaussian integrality and
(3) give

```text
|Z/c| >= constant_d t^3/q^2                           (7)
```

for `q=o(t)`, with the evident harmless adjustment on a fixed initial
range.  If `q=t^(beta+o(1))`, this has exponent `3-2 beta`; it beats the
cost-two profile whenever `beta<1/2`.  This is a real power improvement,
not only a constant improvement.

For three arbitrary positive heights `x<x+p<x+p+q`, the primitive
integer second-difference vector is

```text
r=(q/g,-(p+q)/g,p/g,0),       g=gcd(p,q),               (8)
```

and

```text
sum_j r_j/t_j = p q(p+q)/(g t_1 t_2 t_3),
h(r)=(p+q)/g.                                           (9)
```

This coefficient is forced, not merely convenient.  The integer kernel
of

```text
[[1,1,1],[t_1,t_2,t_3]]
```

has rank one, and (8), up to sign, is its primitive generator.  For four
distinct heights the corresponding affine-moment kernel has rank two and
is generated over the rationals by such three-point circuits.  Rank alone
does not bound an integral generator: its circuit coefficients are the
height gaps divided by their gcd.

Four points can have additional short circuits under an exact gap
condition.  For example,

```text
r=(1,-1,-1,1),       h(r)=4/3,
```

annihilates the constant and linear moments precisely when
`t_1+t_4=t_2+t_3`.  For a translated cluster `t_j=T+a_j` with fixed
offsets, any vector in this affine-moment kernel gives

```text
sum_j r_j/t_j=O_r(T^-3).
```

More generally, if `max_j |a_j|<=H=o(T)`, the bound is
`O(||r||_1 H^2/T^3)`.

Equal pair sums, like three-term arithmetic progressions, are useful exact
gap structure.  Neither condition follows from ordering and comparability.

Thus this standard divided-difference route has cost below three exactly
when `p=q`, in which case it is (2).  A short arithmetic progression is a
real extra hypothesis; four ordered comparable heights do not supply one.

## 2. A primitive endpoint family with no low-cost cancellation

Let

```text
K=3*7*13*31=8463,       M=K(2n+1),
t_j=5^j M  (j=0,1,2,3),       w_j=-1+i t_j.             (10)
```

All four moduli are comparable to `M`, with an absolute comparison ratio
`125`.  They are genuine fixed-real-part Gaussian integers, not formal
prime blocks.

For later use, write

```text
u_j=w_j/(1+i),       n_j=N(u_j)=(1+25^j M^2)/2.
```

The four odd integers `n_j` are pairwise coprime.  Indeed, if an odd prime
`p` divided `n_i` and `n_j`, with `i<j`, then

```text
p | 25^(j-i)-1.
```

Every odd prime divisor of `25^m-1`, for `m=1,2,3`, divides `K`, hence
divides `M`.  But then `2n_i=1+25^iM^2` is `1` modulo `p`, a contradiction.

It follows that the `u_j` are pairwise coprime in `Z[i]`, and

```text
L=lcm_G(w_0,w_1,w_2,w_3)=(1+i) product_j u_j           (11)
```

up to a unit.  This family even comes from primitive circle tuples.  Put

```text
B=(1+i)L,       P=B/2=i product_j u_j,
v_j=B/w_j=(1+i) product_(k != j) u_k,
Q_j=P+v_j.                                               (12)
```

Each `u_j` is ordinary-primitive and has odd norm.  Pairwise coprimality
of their norms prevents a split prime and its conjugate from occurring
in two different factors.  Hence every `v_j` is ordinary-primitive;
multiplication by `1+i` gives it two odd coordinates, not ordinary
content two.  Also the Gaussian gcd of `P,Q_0,...,Q_3` is a unit.  Finally,

```text
2P=v_j w_j,       Re(w_j)=-1,       |Q_j|=|P|.          (13)
```

Thus (12) is a full Gaussian-primitive five-point circle tuple whose four
chords `Q_j-P` all have ordinary content one.

This endpoint realization is deliberately a scope test, not a bounded-arc
counterexample.  Its Euclidean circle radius is

```text
R=|P| asymp M^4,
```

while the angular span from `P` is asymptotic to a constant times `M^-1`.
The containing arc therefore has length `asymp M^3`, whereas `sqrt(R)` is
`asymp M^2`.  Its required endpoint constant in `C sqrt(R)` grows like
`M`.  Thus the family refutes a cancellation lemma based only on four
actual fixed-real-part cofactors (even with primitive complements), but it
does not refute a stronger lemma that also assumes a uniform endpoint-scale
arc bound.

## 3. Balanced base five excludes every cost below three

Let `H(r)=6h(r)`.  If `p_ij=r_i+r_j` and
`q_j=sum_i r_i-r_j`, then, for example,

```text
6r_1 = p_12+p_13+p_14-p_23-p_24-p_34-q_1+q_2+q_3+q_4.
```

The analogous identity holds for every coordinate.  The triangle
inequality therefore gives

```text
h(r) >= max_j |r_j|.                                    (14)
```

If `h(r)<3`, all four entries lie in `{-2,-1,0,1,2}`.  On the family
(10), exact reciprocal cancellation would say

```text
0 = M sum_j r_j/t_j
  = r_0+r_1/5+r_2/25+r_3/125,

0 = 125r_0+25r_1+5r_2+r_3.                             (15)
```

Reduction modulo five successively forces `r_3,r_2,r_1,r_0` to vanish.
There is consequently no nonzero integer reciprocal relation of cost
less than three, even if the vector `r` is allowed to depend on `M`.

There is also no approximate second-order cancellation in this cost
range.  For nonzero `r` as above, the absolute value of the rational
number in the first line of (15) is at least `1/125`.  Since

```text
|atan x-x| <= x^3/3       (0<=x<=1),
```

the exact monomial phase satisfies, throughout (10),

```text
|sum_j r_j atan(1/t_j)| >= 1/(250M).                    (16)
```

Ordinary-content removal and Gaussian gcd cancellation rescale the
associated `Z` by positive real numbers and do not alter this phase.
Thus every cost-below-three monomial in this actual family remains on the
first-order scale.  The hoped-for `T^-3` phase is not forced.

## 4. The Pell endpoint family also has no forced low-cost phase

The geometric family above has arc constant growing with `M`.  The
four-cofactor construction in
[the Pell-family note](vertical_lcm_pell_counterfamily.md) supplies a
stronger scope test on bounded endpoint-scale arcs.  Retain its notation

```text
g=1+2i,        a_n=M_0^n i,        b_n=M_0^n(1+3i),
c_n=M_0^n(4-3i),                   e_n=M_0^n(2+i),
M_0=[[17,4],[4,1]].
```

On positive indices divisible by five, the four cofactors

```text
2g(a_nb_n, a_nc_n, b_nc_n, a_ne_n)                     (17)
```

all have real part `-10`, positive comparable imaginary parts, and
ordinary-primitive complements.  Their lcm has order `T^2`, their circle
radius has order `T^2`, and their containing arc has order `T=sqrt(R)`.

They nevertheless force no second-order relation with `h(r)<3`.  Put

```text
lambda=9+4sqrt(5),       s=sqrt(5)-2.
```

The expanding projections of the four initial vectors, after dropping a
common positive normalization, are

```text
A=s,       B=1+3s,       C=4-3s,       E=2+s.
```

The four imaginary parts in (17), divided by `lambda^(2n)`, tend to one
common nonzero constant times `AB,AC,BC,AE`.  Hence a phase of order
`o(T^-1)` would require

```text
r_1 CE+r_2 BE+r_3 AE+r_4 BC=0.                          (18)
```

Using `s^2+4s-1=0`, the rational and `s` coefficients in (18) are

```text
5r_1+5r_2+r_3-5r_4=0,
4r_1+r_2+7r_4=0.                                       (19)
```

Thus

```text
r_2=-4r_1-7r_4,       r_3=15r_1+40r_4.                (20)
```

For `h(r)<3`, (14) gives `|r_j|<=2`; (20) then forces `r_1=r_4=0`
and hence `r=0`.  There are only finitely many vectors in this cost range,
so convergence of the four ratios upgrades the pointwise statement to a
uniform one: for all sufficiently large Pell indices,

```text
|sum_j r_j/t_j| >= c/T                                (21)
```

for every nonzero `r` with `h(r)<3`, with one absolute `c>0`.  The
`O_r(T^-3)` arctangent remainder cannot change this first-order scale.

Here `h` remains the formal ten-block benchmark (1).  This argument does
not assert that the actual Pell monomial has modulus cost `T^h`; it says
that none of the coefficient vectors which would have formal cost below
three supplies the desired reciprocal cancellation.

This rules out a forcing lemma based on four endpoint-scale fixed-real-part
cofactors and primitive complements alone.  It does not test an additional
subquadratic-lcm or five-thirds-sharing hypothesis: this Pell family's lcm
already has the sharp quadratic order.  Its normalized arc constant tends
to one fixed positive value, so it also does not refute a lemma restricted
to a sufficiently small prescribed endpoint constant or to constants
tending to zero.

## 5. Why positive divided differences do not change the conclusion

For positive distinct nodes, the reciprocal function has the exact law

```text
[x_0,...,x_k](1/x)=(-1)^k/(x_0 x_1 ... x_k).            (22)
```

All divided differences of a fixed order therefore have one strict sign;
a positive combination of them cannot vanish.  Clearing their barycentric
denominators can produce an integer signed vector, but (8)--(9) show the
coefficient cost explicitly in order two.  On the geometric nodes (10),
the four possible triples have primitive costs `6,31,31,6`, respectively.
None enters the range `h<3`, and their reciprocal values remain of order
`M^-1`, rather than `M^-3`.

The same obstruction works at any prescribed bounded cost: choose an odd
base larger than twice that bound, replace `5^j` by its first four powers,
and enlarge `M` to contain the odd prime divisors of the three relevant
norm differences.  Balanced-base uniqueness rules out every relation in
the chosen coefficient box.  Hence real-rootedness or positivity alone
cannot force a uniform bounded-height reciprocal cancellation.  Any
successful uniform argument needs additional clustering or arithmetic
sharing beyond four fixed-real-part cofactors.

Run `python3 docs/check_vertical_second_order_phase_scope.py` for exact
checks of the low-cost lattice, divided differences, Gaussian lcm,
ordinary primitivity, and circle identities.  The arguments above prove
the infinite statements.
