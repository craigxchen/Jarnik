# The least circle for an actual balanced three-row family

This is an independent radius and arc audit for the seven-factor
three-row family below. It constructs four actual integer points,
including the anchor, on their **least** common circle. The family
has one full Boolean norm block for every nonempty subset of three
rows, all at the same logarithmic height. It supplies four points
on an arc of order `rho^(3/7)` for a circle of radius `rho`.
The number of points is fixed at four; this is not a counterexample
to a uniform bound on the number of points.

## 1. Exact rows and a safe progression

For an integer `t`, define the Gaussian linear factors

```text
H_123=t-i,       H_12=t-6-5i,    H_13=t-6+3i,
H_23=t-6+4i,    H_1=t-4+3i,     H_2=t-3+2i,
H_3=t-7-6i.
```

Put `P_i=product_(T contains i) H_T`. Direct multiplication gives

```text
P_1=t^4-16t^3+106t^2-256t+105+240i,
P_2=t^4-15t^3+ 95t^2-195t+ 94+180i,
P_3=t^4-19t^3+151t^2-439t+150+420i.             (1)
```

In particular the constant terms are `105+240i`, `94+180i`,
and `150+420i`. Their ordinary coordinate gcds are `15,2,30`.
Let `g=(15,2,30)` and

```text
K_1=7+16i,       K_2=47+90i,       K_3=5+14i,
J_T(t)=H_T(t)/H_T(0),
R_i(t)=P_i(t)/g_i=K_i product_(T contains i) J_T(t).    (2)
```

For a concrete infinite progression take `t=Mh`, `h>=1`, where

```text
M=8088851149545216141000
 =2^3 3^3 5^3 7 13^2 17^2 29 37 41 53 61^2 101.
```

Every `J_T(t)` is Gaussian integral because `Norm(H_T(0))|M`.
Write `Q_T(t)=Norm(H_T(t))` and `n_T=Norm(J_T(t))
=Q_T(t)/Q_T(0)`. The bad-prime set used in `M` is

```text
S={2,3,5,7,13,17,29,37,41,53,61,101}.
```

It contains all prime divisors of the seven constant norms,
the fixed imaginary coordinates, the three `Norm(K_i)`, and every
pair resultant of the distinct quadratics `Q_T`. For factors
`H_T=t-a_T+i b_T` and `H_U=t-a_U+i b_U`, that resultant is

```text
((a_T-a_U)^2+(b_T-b_U)^2)
((a_T-a_U)^2+(b_T+b_U)^2),
```

which makes the finite prime check explicit. The exponent of `M` at
each `p` is greater than every `v_p(Q_T(0))`. Hence
`v_p(Q_T(Mh))=v_p(Q_T(0))` for `p in S`, so every `n_T` is
`S`-free. A common prime of two `n_T` would divide the corresponding
quadratic resultant and belong to `S`, which is impossible. Thus the
seven `n_T` have pairwise disjoint rational-prime support and are
coprime to every `Norm(K_i)`. Each `J_T` has odd norm and is
conjugate-primitive: a rational prime dividing both its coordinates
would divide the fixed imaginary coordinate of `H_T(t)` and hence
belong to `S`. Consequently the `R_i` are conjugate-primitive and

```text
Im(R_1,R_2,R_3)=(16,90,14),
log n_T=2 log t+O(1)       for all seven T.           (3)
```

The three bracket identities, checked by polynomial multiplication,
are

```text
Im(bar(P_1)P_2)=-60 Q_123 Q_12,
Im(bar(P_1)P_3)=180 Q_123 Q_13,
Im(bar(P_2)P_3)=240 Q_123 Q_23.
```

After (2) they become

```text
Delta_12=-122 n_123 n_12,
Delta_13=  18 n_123 n_13,
Delta_23= 208 n_123 n_23.                         (4)
```

The pair Gaussian gcd norms of the constant corrections are
`Norm gcd(K_1,K_2)=61`, `Norm gcd(K_1,K_3)=1`, and
`Norm gcd(K_2,K_3)=13`. Therefore the primitive pair numerators
have constant correcting factors

```text
K_12=(K_2 bar(K_1))/61=29-2i,
K_13=K_3 bar(K_1)=259+18i,
K_23=(K_3 bar(K_2))/13=115+16i.                     (5)
```

Their imaginary residues are `-2,18,16`, respectively, after the
pair core is included; equivalently (4) has the stated residuals.
Because the constant correction norms avoid every core norm,
the core numerator divides each actual primitive pair numerator.

## 2. The least common circle is explicit

Let `B=product_(T!=empty) J_T` and, for `a=1,2,3`, set

```text
B_0=bar(B),
B_a=product_(T contains a) J_T
    product_(T not containing a) bar(J_T).
```

The Gaussian lcm of the conjugate constant corrections is, up to a
unit,

```text
E=lcm_G(bar(K_1),bar(K_2),bar(K_3))=912-211i,
Norm(E)=876265=5*13^2*17*61.                         (6)
```

Indeed each `bar(K_i)` divides the displayed `E`, and the pair gcd
norms above give the same lcm norm. Put

```text
E_0=E,
E_1=E K_1/bar(K_1)=-464+813i,
E_2=E K_2/bar(K_2)=-348+869i,
E_3=E K_3/bar(K_3)=-572+741i,
z_a=B_a E_a       (0<=a<=3).                         (7)
```

All four `z_a` are Gaussian integers, with exactly equal norm

```text
rho(t)^2=Norm(z_a)=Norm(E) product_(T!=empty)n_T(t)
        =(1/4500) product_(T!=empty)Q_T(t).           (8)
```

The last equality uses the exact arithmetic
`product_T Q_T(0)=3943192500=4500*876265`.
Moreover

```text
z_i/z_0=R_i/bar(R_i).                                (9)
```

This is the **least** integral circle for these four phases. Since
the Gaussian prime supports of the `J_T` are mutually disjoint and
avoid those of the `K_i`, the denominator lcm at anchor zero is

```text
lcm_G(bar(R_1),bar(R_2),bar(R_3))=bar(B)E=z_0
```

up to a unit. Any integral starting point whose rotations by the
three reduced phases in (9) remain integral is divisible by this
lcm. The four points in (7) have Gaussian gcd one, so changing the
anchor gives the same least radius, as in the general
[all-anchor radius identity](all_anchor_rotation_radius_cocycle.md).

## 3. Arc scale

The real parts of the `P_i` in (1) are `t^4+O(t^3)`, positive for
large positive `t`. Therefore the relative angles of the four
points in (7), with `z_0` at angle zero, are

```text
theta_1=2 arctan(240/Re(P_1))=480 t^(-4)+O(t^(-5)),
theta_2=2 arctan(180/Re(P_2))=360 t^(-4)+O(t^(-5)),
theta_3=2 arctan(420/Re(P_3))=840 t^(-4)+O(t^(-5)).
```

Thus the cyclic order on the short arc is `z_0,z_2,z_1,z_3` for
large `h`, and all four points are distinct. From (8),

```text
rho(t)=t^7/(30 sqrt(5))*(1+O(t^(-1))),
arc_length=rho(t) theta_3
          =(28/sqrt(5)) t^3*(1+O(t^(-1)))
          =Theta(rho(t)^(3/7))=o(sqrt(rho(t))).         (10)
```

Consequently every fixed positive multiple of `sqrt(rho)`
eventually contains these four points. Their number stays four.

## Verification

The [family checker](check_balanced_three_row_gaussian_polynomial_family.py)
multiplies the Gaussian factors, verifies the bracket polynomial
identities, and checks the progression, correction gcds and lcm,
equal norms, phase ratios, Gaussian gcd one, and the exact radius
identity on eight positive multiples of `M`. The proofs above apply
to every positive multiple of `M`.
