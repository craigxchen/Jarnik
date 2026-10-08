# Moving-base Gram reduction and weighted quartic content

This note records the exact arithmetic when an actual six-node clique is
reanchored at infinity and two finite cotangents. The final scale statement
is conditional on the size of the remaining weighted evaluation matrix; no
uniform fair asymptotic is assumed.

## 1. The canonical primitive Gram form

Let the finite cotangents be X_i/L > 0, with L > 0, and choose two of
them as A/L and B/L, where

```text
delta=B-A>0.
```

The inverse normalization sending the three nodes

```text
infinity, A/L, B/L  ->  infinity, 0, 1
```

is

```text
t=(delta*s+A)/L,
M=[[delta,A],[0,L]].
```

The integer Gram matrix is

```text
G=M^T M=[[delta^2,delta*A],[delta*A,A^2+L^2]].       (1)
```

The edge condition for A/L,B/L says

```text
delta | A^2+L^2.
```

Put C=(A^2+L^2)/delta and

```text
d=gcd(delta,A,C).
```

Then

```text
gcd(all entries of G)=delta*d.                       (2)
```

Moreover d divides L. For completeness, let p^e be the exact p-part of d.
Since p^e divides delta,A,C, the right side of

```text
A^2+L^2=delta*C
```

has p-valuation at least 2e. If v_p(L)<e, then
v_p(A^2)>=2e>2v_p(L), so v_p(A^2+L^2)=2v_p(L)<2e, a contradiction.
This also covers p=2.

Thus the canonical primitive positive Gram form is exactly

```text
Q0=G/(delta*d)
  =[[delta/d,A/d],[A/d,C/d]],                         (3)
det(Q0)=(L/d)^2.
```

The entries of Q0 have gcd one, and (3) is therefore the primitive
integral representative of the projective Gram point. The normalized
remaining cotangents are

```text
x=(X-A)/delta,
```

with x=0,1 for the selected pair and three further distinct values
x,y,z.

## 2. Exact weighted evaluation content

Let T=(0,1,x,y,z), and for t in T define

```text
D_t=product_(r in T, r!=t) (t-r),
E[t,k]=t^k/D_t,  0<=k<=4.
```

Let ell be the least positive integer for which

```text
N=ell*E
```

has integral entries. Append the infinity row to obtain the integer 6-by-5
matrix

```text
Nhat = [ N ]
       [ 0  0  0  0 -ell ].                           (4)
```

Write Q0=(a,b,c), meaning the symmetric matrix with diagonal a,c, and set

```text
P0(s)=a*s^2+2*b*s+c,
f=(c^2,4*b*c,2*a*c+4*b^2,4*a*b,a^2)^T.              (5)
```

The six raw central weights after this normalization are

```text
(Nhat*f)/ell.                                         (6)
```

The last coordinate is -a^2, the negative leading coefficient of P0^2;
the Lagrange identity gives that it is the required infinity weight.
Consequently the exact primitive integer central vector is

```text
u=(Nhat*f)/C_eval,
C_eval=gcd of the six entries of Nhat*f,              (7)
U=||Nhat*f||_infinity/C_eval.                         (8)
```

This is the weighted quartic evaluation gcd. It includes every denominator
and every accidental specialization cancellation. The Gram matrix in (1)
is (delta*d)Q0, so replacing Q0 by G multiplies f by (delta*d)^2, which
cancels from the primitive vector (7). The ordinary coefficient content

```text
kappa=gcd(f)
```

is at most 4 for primitive (a,b,c); hence one may replace f by f/kappa in
(7) and regard the remaining gcd as purely weighted evaluation content.

For reference, if N5 denotes the first five rows of N, then its
invertibility gives the quantitative but tuple-dependent comparison

```text
||f||_infinity/||N5^(-1)||_infinity
 <= ||Nhat*f||_infinity
 <= ||Nhat||_infinity*||f||_1.                        (9)
```

Since ||f||_infinity is comparable to H(Q0)^2 with absolute constants,
(8)-(9) identify exactly where the fixed-base estimate loses information:
it is the integer C_eval, together with the explicit norms of the moving
weighted matrix N.

## 3. A conditional fair-scale consequence

Suppose along a family that

```text
log|A|=log|B|=log(delta)=8w+o(w),
log L=o(w),
log U=7w+o(w),                                    (10)
```

and A,B,delta > 0. Since d divides L, (3) gives

```text
log H(Q0)=8w+o(w).                                  (11)
```

Indeed delta/d and A/d already have this size, while
C/d=(A^2+L^2)/(delta*d) has the same logarithm.

Define the explicit evaluation scale

```text
R_eval=||Nhat*f||_infinity/H(Q0)^2.
```

The exact identity (8) then yields

```text
log C_eval=log R_eval+9w+o(w).                       (12)
```

Equivalently, without any assumption on the moving matrix,

```text
log(C_eval/R_eval)=9w+o(w).
```

This ratio is the exact moving-base content/height loss relative to the
weighted evaluation scale.

Thus, if the normalized cross-ratios satisfy the additional checked
condition log R_eval=o(w), the weighted evaluation content costs
9w+o(w) of the nominal H(Q0)^2 scale. Without that condition, (12) is
still the rigorous statement: the deficit is 9w+o(w) after subtracting
the explicitly computable moving-matrix evaluation scale.
Typical cotangents and pair differences of size exp(8w+o(w)), and
log L=o(w), establish (11), but do not by themselves establish
log R_eval=o(w). That last comparison must be verified from the three
remaining cross-ratios and their weighted denominators.
