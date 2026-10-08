# An unbounded balanced three-row Gaussian family

There is an exact unbounded family realizing the entire seven-block
Boolean Gaussian profile on three rows, with conjugate-primitive rows,
pairwise disjoint rational-prime block supports, bounded correcting
factors, bounded nonzero imaginary coordinates, and bounded nonzero
determinant residuals. All seven block norms have logarithm
`w+O(1)` as `w` tends to infinity.

Thus a universal positive residual-height gap cannot hold for three
rows in the endpoint reduction's general, nonconstant-`Y_i` scope.
This does not settle the uniform endpoint problem or produce clusters
of unbounded cardinality. It does not settle the special `Y_i=1`
interpolation problem either. The family has three fixed distinct
imaginary coordinates.

## 1. Seven factors and exact polynomial identities

For a real or integer variable `t`, put

```text
H_g=t-i,              H_12=t-6-5i,
H_13=t-6+3i,          H_23=t-6+4i,
H_1=t-4+3i,           H_2=t-3+2i,
H_3=t-7-6i.
```

The label `g` denotes the full subset `{1,2,3}`. Define the rows by
their exact incident factors:

```text
P_1=H_g H_12 H_13 H_1,
P_2=H_g H_12 H_23 H_2,
P_3=H_g H_13 H_23 H_3.
```

Expansion gives

```text
P_1=t^4-16t^3+106t^2-256t+105+240i,
P_2=t^4-15t^3+ 95t^2-195t+ 94+180i,
P_3=t^4-19t^3+151t^2-439t+150+420i.                     (1)
```

In particular the three imaginary parts are nonzero constants. Writing
`Q_T=Norm(H_T)`, direct multiplication also gives

```text
Im(bar(P_1)P_2)=-60 Q_g Q_12,
Im(bar(P_1)P_3)=180 Q_g Q_13,
Im(bar(P_2)P_3)=240 Q_g Q_23.                            (2)
```

Every `Q_T(t)` is a positive monic quadratic. Consequently all three
determinants in (2) are nonzero for every real `t`.

The constant values and their coordinate contents are

```text
P_1(0)=105+240i,       g_1=15,
P_2(0)= 94+180i,       g_2= 2,
P_3(0)=150+420i,       g_3=30.
```

The resulting fixed correcting factors are

```text
K_1=7+16i,       Norm(K_1)=305=5*61,
K_2=47+90i,     Norm(K_2)=10309=13^2*61,
K_3=5+14i,      Norm(K_3)=221=13*17.                    (3)
```

Each `K_i` has coprime integer coordinates of opposite parity, so
each is conjugate-primitive.

## 2. One progression supplies integral independent primitive blocks

Normalize every block by its constant Gaussian value:

```text
J_T(t)=H_T(t)/H_T(0),       n_T(t)=Norm(J_T(t)).
```

Use the fixed positive integer

```text
M=8088851149545216141000
 =2^3*3^3*5^3*7*13^2*17^2*29*37*41*53*61^2*101,
t=Mh,                      h=1,2,3,... .                (4)
```

For every such `t`, all `J_T(t)` are Gaussian integers. Their
norms are odd, exceed one, and have pairwise disjoint rational-prime
support. Each `J_T(t)` is conjugate-primitive, and every `n_T(t)`
is coprime to all three fixed norms in (3).

Here is a proof valid for the entire progression. Write
`H_T(t)=t+a_T+i b_T`, and let `c_T=a_T^2+b_T^2=Q_T(0)`.
In the order `g,12,13,23,1,2,3`, the constants are

```text
c_T=(1,61,45,52,25,13,85).                              (5)
```

Let

```text
S={2,3,5,7,13,17,29,37,41,53,61,101}.                   (6)
```

This set contains every prime dividing one of the following:

- an imaginary constant `b_T`;
- a constant norm `c_T`;
- a fixed correcting norm `Norm(K_i)`;
- the resultant of two different norm quadratics `Q_T,Q_U`.

All twenty-one resultants are nonzero. For clarity, if the quadratics
are `t^2+b t+c` and `t^2+d t+e`, their resultant is exactly

```text
c(d-b)^2-b(d-b)(e-c)+(e-c)^2.                            (7)
```

Applying (7) to the displayed quadratics gives the set (6). The
[exact checker](check_balanced_three_row_gaussian_polynomial_family.py)
verifies this finite integer certificate independently of evaluation
at any `t`.

At each `p in S`, the integer `M` in (4) contains the power

```text
p^(1+max_T v_p(c_T)).                                   (8)
```

Thus `Q_T(t)=c_T+2a_T t+t^2` satisfies

```text
v_p(Q_T(t))=v_p(c_T),       p in S,       t=Mh.           (9)
```

Also every `c_T` divides `M`. The exact coordinates

```text
J_T(t)=1+(a_T t/c_T)-i(b_T t/c_T)                        (10)
```

are therefore integers, and

```text
n_T(t)=Q_T(t)/c_T
```

has no prime divisor in `S` by (9). If two such norms had a common
prime outside `S`, their original quadratics would vanish at `t`
modulo that prime; it would divide their nonzero resultant, contrary
to (6). This proves disjoint rational-prime support.

To prove primitivity, suppose a rational prime `p` divided both
coordinates of `J_T(t)`. It cannot lie in `S`, since its norm is
`S`-free. Multiplication by `H_T(0)` would make `p` divide both
coordinates of `H_T(t)`, including its fixed imaginary coordinate
`b_T`, which again forces `p in S`. Thus the integer coordinates
are coprime. Their norm is odd because `2 in S`, so the Gaussian
integer is conjugate-primitive. This also shows directly that the
norm contains no inert prime. Equation (6) finally proves disjointness
from every correcting norm in (3).

All constants `a_T` have absolute value at most seven, and `M>14`.
Hence `Q_T(t)-Q_T(0)=t(t+2a_T)>0` throughout (4), proving
`n_T(t)>1` with no excluded initial values.

## 3. The normalized rows retain every equation

Define

```text
R_1=K_1 J_g J_12 J_13 J_1,
R_2=K_2 J_g J_12 J_23 J_2,
R_3=K_3 J_g J_13 J_23 J_3.
```

The normalization constants cancel exactly, giving

```text
R_i(t)=P_i(t)/g_i,
(Im R_1,Im R_2,Im R_3)=(16,90,14).                     (11)
```

These are Gaussian integers by their product definitions. Within each
row, all block norms are pairwise coprime and coprime to the correcting
norm. Each factor is conjugate-primitive. Therefore the entire row
is conjugate-primitive, and the listed imaginary coordinates are
nonzero for every `h>=1`.

The exact row norm equations are

```text
Norm(R_i)=Norm(K_i) product_(T containing i)n_T.           (12)
```

Dividing (2) by `g_i g_j` and using (5) yields all three actual
determinant equations:

```text
Im(bar(R_1)R_2)=-122 n_g n_12,
Im(bar(R_1)R_3)=  18 n_g n_13,
Im(bar(R_2)R_3)= 208 n_g n_23.                            (13)
```

For example, the first coefficient is
`-60 c_g c_12/(g_1 g_2)=-60*61/30=-122`. The other two are
`180*45/450=18` and `240*52/60=208`. Thus the determinant
residuals are the fixed nonzero integers

```text
(t_12,t_13,t_23)=(-122,18,208).
```

Equations (11)--(13) hold for the same three rows and the same seven
blocks. In particular the rank-two identities and all triple norm-metric
identities are automatically satisfied. This is an actual Gaussian
factorization, rather than a model for isolated necessary circuits.

## 4. The balanced limit and its scope

Put `w=2 log t`. For each of the seven fixed quadratics,

```text
log n_T(t)=2 log t-log c_T+O(1/t)=w+O(1),
log|R_i(t)|=4 log t-log g_i+O(1/t)=2w+O(1).              (14)
```

Meanwhile `K_i`, `Im R_i`, and every `t_ij` are fixed nonzero
constants. Their logarithmic heights are `O(1)=o(w)`. Thus for
every fixed positive `delta` and every height threshold, sufficiently
large members satisfy the complete three-row profile within relative
block error `delta` and with all residual heights below `delta w`.

This rules out a universal three-row exclusion theorem of the form
stated in the [endpoint quantifier audit](endpoint_full_profile_quantifiers.md),
even with exact row primitivity and coherent Gaussian orientations.
Any universal exclusion using that full-profile formulation must use
at least four nonanchor rows or additional hypotheses beyond this
three-row profile. A claim restricted to equal imaginary coordinates,
such as the `Y_i=1` setup in the
[three-row interpolation note](three_row_discriminant_interpolation_scope.md),
has different hypotheses and is not refuted by this family.

The number of rows here is fixed. The construction does not establish
unbounded endpoint counts, solve the four-row fifteen-block system, or
disprove the radius-uniform conjecture.

The checker verifies the polynomial products and determinant identities,
the exact resultant and progression certificates, and eight entire
normalized configurations. These samples check implementation; equations
(7)--(10) prove the prime support and primitivity claims for every
positive `h`.

## 5. Four points on their least integral circle

The three conjugate-primitive phases `R_i/bar(R_i)` have an explicit
least common circle. Let `B=product_T J_T`, over all seven nonempty
cuts, and set

```text
B_0=bar(B),
B_i=product_(T contains i) J_T product_(T not containing i) bar(J_T).
```

The Gaussian lcm of `bar(K_1),bar(K_2),bar(K_3)` is, up to a unit,
`E=912-211i`, of norm `876265`. Put

```text
E_0=E,
E_1=E K_1/bar(K_1)=-464+813i,
E_2=E K_2/bar(K_2)=-348+869i,
E_3=E K_3/bar(K_3)=-572+741i,
z_a=B_a E_a                  (0<=a<=3).                 (15)
```

The four `z_a` are Gaussian integers with
`z_i/z_0=R_i/bar(R_i)` and one common squared radius

```text
rho(t)^2=876265 product_T n_T(t)
        =product_T Norm(H_T(t))/4500.                   (16)
```

This circle is least: disjointness of the block and correction prime
supports gives
`lcm_G(bar(R_1),bar(R_2),bar(R_3))=bar(B)E=z_0`.
Every integral anchor point realizing all three reduced phases is
divisible by this lcm. The four points have Gaussian gcd one, so all
reanchorings have the same least radius.

The real parts in (1) are `t^4+O(t^3)` and their imaginary parts are
`(240,180,420)`. Hence the relative point angles are, in row order,
`(480,360,840)t^(-4)+O(t^(-5))`. For large `h`, the four points
occur on one short arc in order `z_0,z_2,z_1,z_3`, with

```text
rho(t)=t^7/(30 sqrt(5))*(1+O(t^(-1))),
arc_length=(28/sqrt(5))t^3*(1+O(t^(-1)))
          =Theta(rho(t)^(3/7))=o(sqrt(rho(t))).          (17)
```

Thus any fixed positive multiple of `sqrt(rho)` eventually contains
these four points. Their number remains four. The exact Gaussian
lcm, pair corrections, and all-anchor accounting are given in the
[least-circle audit](three_row_balanced_family_least_circle.md).
The family checker also verifies the exact constant lcm, all four
integral points, their common norm and Gaussian gcd one, and the phase
identities (15)--(16) at eight progression values.
