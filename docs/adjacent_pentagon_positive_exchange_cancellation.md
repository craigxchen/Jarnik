# Exact cancellation in the positive exchange of two adjacent pentagons

Two adjacent five-point pentagons have a genuinely joint positive exchange.
Its primitive cancellation is explicit. At a prime supporting both gap
factors, all prescribed gap powers disappear from the common primitive
numerator of the two positive summands. Only the squared integer content
of the separately normalized shared interior chord remains.

This identifies a precise obstruction to using that common numerator as
a stronger gap divisor. It does not rule out other uses of the positive
exchange, and it gives no improved endpoint exponent.

## 1. The six-point positive exchange

Let `z_0,...,z_5` be six distinct cyclically ordered Gaussian integers of
common norm `N`, and write `l_ij=|z_i-z_j|`. Use the three consecutive
four-point ratios

```text
x=l_03 l_12/(l_02 l_13)=a/b,
y=l_14 l_23/(l_13 l_24)=c/d,
z=l_25 l_34/(l_24 l_35)=e/f,
```

each in lowest positive terms. The letter `z` without a subscript denotes
the third ratio only. All three ratios belong to `(0,1)` and are rational,
including for differing Gaussian row units. Define

```text
E_L=ad+bc-bd,                  E_R=cf+de-df,
T=ade+bcf-bdf.
```

Order of all six points gives the stronger compatibility

```text
xz+y>1,                       T>0.                    (1)
```

In particular both pentagon conditions `E_L>0,E_R>0` follow. One way to
see the additional condition is to put the projective coordinates of rows
`1,2,3,4` at `0,1,t,infinity`, where `t=1/(1-y)`. The other two
coordinates are `-s,-r`, with

```text
s=(1-x)/(x+y-1),
r=(z+y-1)/((1-y)(1-z)).
```

The full cyclic order is equivalent to `r>s>0`; after clearing the
positive denominators, `r>s` is exactly (1). Conversely (1), with
`0<x,y,z<1`, supplies that order.

Direct expansion now gives the positive integer exchange

```text
E_L E_R = d(d-c)(b-a)(f-e) + cT.                       (2)
```

For the separately normalized pentagons of
[the primitive-content note](positive_pentagon_primitive_content.md), put

```text
L_L=gcd(b,d) gcd(a,d-c) gcd(c,b-a),
L_R=gcd(d,f) gcd(c,f-e) gcd(e,d-c).
```

Thus `E_L=L_L F_L`, `E_R=L_R F_R`, and (2) couples the actual two
pentagon contents, with their full cancellation factors retained.

The geometric identity behind (2) is especially useful. Put

```text
q_L=l_04 l_23/(l_03 l_24)=E_L/(ad),
q_R=l_15 l_23/(l_13 l_25)=E_R/(de),
W=q_L q_R=E_L E_R/(ae d^2),
r_0=l_01 l_45/(l_04 l_15).
```

Ptolemy on rows `0,1,4,5` says

```text
W=A+B,
A=Wr_0   =l_01 l_45 l_23^2/(l_03 l_24 l_13 l_25),
B=W(1-r_0)=l_05 l_14 l_23^2/(l_03 l_24 l_13 l_25).   (3)
```

Both summands are positive rational numbers. Equation (2) is (3) with
its explicit common rational denominator cleared.

## 2. Complete evaluated primitive cancellation

Write `W=P/Q` and `r_0=h/k` in lowest positive terms. Let `num` denote
the positive reduced numerator. Then exactly

```text
J:=gcd(num(A),num(B))=P/gcd(P,k).                       (4)
```

Indeed let `g=gcd(P,k)`, `u=gcd(h,Q)`, and
`v=gcd(k-h,Q)`. Since `gcd(h,k-h)=1` and `gcd(P,Q)=1`, cancellation
in the two fractions gives

```text
num(A)=(P/g)(h/u),
num(B)=(P/g)((k-h)/v).
```

The two factors after `P/g` are coprime. This proves (4), including
all repeated prime powers.

There is a formula directly in the three reduced consecutive ratios:

```text
M=E_L E_R,
U=gcd(M,ae d^2),
S=gcd(M,d(d-c)(b-a)(f-e)),
P=M/U,                 k=M/S,
J=S/gcd(S,U).                                          (5)
```

No squarefree assumption or omitted cancellation factor occurs here.
Alternatively, if `q_L=n_L/d_L`, `q_R=n_R/d_R` are separately reduced,
the exact cancellation before (4) is

```text
G=gcd(n_L,d_R) gcd(n_R,d_L),
P=n_L n_R/G,            Q=d_L d_R/G.                    (6)
```

Thus even the first multiplication has a possible cancellation against
the other pentagon's denominator. At primes where both gaps are genuine,
that first cancellation does not remove either gap; the subsequent
factor `gcd(P,k)` in (4) does.

## 3. Exact loss of both gap powers

Fix an odd split source prime `p=pi conjugate(pi)` dividing `N`, and
put `t_i=v_pi(z_i)`. Separate normalization of each five-row subset
does not affect its gap. Let their positive gaps be `h_L,h_R>0`.
The crossing endpoint incidences force opposite orientations. After
possibly exchanging `pi` and its conjugate, the levels obey

```text
min(t_0,t_4) > max(t_2,t_3) >= min(t_2,t_3) > max(t_1,t_5).
```

Write

```text
A_0=min(t_0,t_4),       B_0=max(t_1,t_5),
m=min(t_2,t_3),         M_0=max(t_2,t_3),
w=M_0-m,
h_L=A_0-M_0,           h_R=m-B_0.
```

For a pair at equal allocation level put
`rho_ij=v_pi(z_i-z_j)-t_i`; put `rho_ij=0` for unequal levels.
Unequal-level differences have exactly their minimum valuation. At equal
levels the excess at the conjugate prime is also `rho_ij`, as explained
in the primitive-content note; arbitrary row units are included.
Direct valuation of the rational expressions in (3) yields

```text
v_p(q_L)=h_L+rho_04+rho_23,
v_p(q_R)=h_R+rho_15+rho_23,
v_p(P)=h_L+h_R+rho_04+rho_15+2rho_23,
v_p(k)=h_L+h_R+w+rho_04+rho_15.                         (7)
```

For the last line, the numerator of `r_0` has valuation `t_1+t_5`,
while its denominator has valuation
`A_0+rho_04+min(t_1,t_5)+rho_15`. Thus `r_0` has strictly negative
valuation, of absolute value `A_0-B_0+rho_04+rho_15`.

By (4), and because `rho_23=0` whenever `w>0`,

```text
v_p(J)=max(0,2rho_23-w)=2rho_23.                        (8)
```

If `G_23=gcd_G(z_2,z_3)` and `c_23` is the ordinary integer content
of `(z_3-z_2)/G_23`, then `v_p(c_23)=rho_23`. Consequently

```text
v_p(J)=2v_p(c_23).                                     (9)
```

The forced contribution `h_L+h_R` and both outer chord excesses have
disappeared. This is exact on the shared gap support, not a global
claim that `J=c_23^2`: primes supporting only one gap, or neither gap,
may also occur.

Because the left interior levels `t_1,t_2,t_3` are nonconstant, and so
are the right interior levels `t_2,t_3,t_4`, the earlier exact pentagon
formula gives here

```text
v_p(F_L F_R)=h_L+h_R+rho_04+rho_15,
v_p(k)=v_p(F_L F_R)+w.                                 (10)
```

This pinpoints the remaining denominator that prevents the primitive
common-summand content from retaining the product of gap factors.

## 4. Actual ordered fixtures and size accounting

The primitive ordered six-tuple already used in
[the oriented-overlap note](ordered_overlapping_gap_factors.md),

```text
N=14365,
((-107,-54),(-98,-69),(-91,-78),(-78,-91),(-69,-98),(-54,-107)),
```

has

```text
(x,y,z)=(5/8,29/34,5/8),    F_L=F_R=K_L=K_R=13,
q_L=q_R=13/17,             W=169/289,
r_0=153/1690,
A=9/170,                  B=1537/2890,
J=1.
```

Thus the product `169` is present in the reduced numerator of `W`
and is entirely absorbed by the outer Ptolemy denominator in (4).
This is an actual equal-radius integer realization, with all six rows
strictly ordered on a minor arc, not independent rational pentagons.

Two further primitive ordered first-quadrant fixtures test the two
different cases in (8), including nested source powers:

```text
N=105625=5^4 13^2,
((323,36),(312,91),(300,125),(280,165),(260,195),(253,204)),
pi=2+i: levels=(0,4,2,3,1,4),  h_L=h_R=1, w=1, rho_23=0;

N=1373125=5^4 13^3,
((1170,65),(1167,106),(1150,225),(650,975),(538,1041),(510,1055)),
pi=2+i: levels=(1,4,2,2,0,3),  h_L=h_R=1, w=0, rho_23=1.
```

The last fixture has `v_5(J)=2`; this comes solely from shared chord
content, not from the two gap exponents. All these fixtures are outside
the endpoint class `C=1/2`.

For the overlap-supported product
`K_cap=product_(h_L,h_R>0) p^(h_L+h_R)`, (7) does give `K_cap|k`.
If the global tuple is primitive, parity gives

```text
k <= l_04 l_15/2 <= N Delta_04 Delta_15/2.
```

On `Delta<=C N^(-1/4)` this is `O(N^(1/2))`. It does not improve
the exponent from multiplying the two existing interior-triangle
bounds. In fact the latter give, with `c_0=l_23`,

```text
K_L K_R <= c_0^2 l_12 l_13 l_24 l_34/(16N)
        <= (C^4/16)c_0^2,
```

as recorded in the oriented-overlap note. The weaker `k` estimate
applies only to `K_cap`, not automatically to all of `K_L K_R`.
On a minor arc `l_04,l_15>=c_0`, so the right side of the direct
bound `k<=l_04 l_15/2` is at least `c_0^2/2`. At `C=1/2`, the
triangle-product bound is `c_0^2/256`. Thus the uncanceled source
divisor in `k` comes with a weaker archimedean bound, even before
the loss from restricting to the common prime support.

The positive exchange still contains arithmetic cancellation between
its summands modulo large source powers. Turning that cancellation into
a radius gain requires an additional bound on their reduced residue
sizes or on the outer denominator `k`; positivity and the exact common
numerator (4) do not supply it.

The [exact checker](check_adjacent_pentagon_positive_exchange.py) verifies
the exchange and simultaneous cancellations over bounded rational
triples and verifies (7)--(10) on all three actual fixtures, together
with their order, primitivity, and exclusion from `C=1/2`.
