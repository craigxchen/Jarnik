# Ordered terminal circuits do not force positivity

This exact six-row example rules out a sign argument based only on ordered
half-angle slopes, the terminal circuit equations, and numerical
reconstruction. It is an algebraic countermodel, not an endpoint
counterexample to the core-profile short-height theorem.

Take the ordered slopes `b_i=i`, `0<=i<=5`, and define

```text
L_abc(x)=(b_c-b_b)x_a-(b_c-b_a)x_b+(b_b-b_a)x_c,
f=9 L_012 L_345-L_013 L_245.
```

Every `L_abc` is killed by both `D_1=sum partial_i` and
`D_b=sum b_i partial_i`. Therefore `D_1 f=D_b f=0`. Writing
`f=sum_(i<j)c_ij x_i x_j`, these identities are exactly the terminal
integer circuit equations for the binary rows `(X_i,Y_i)=(1,b_i)`:

```text
sum_(j!=i) c_(ij) (X_j,Y_j)=(0,0) for each i.
```

The full nonzero coefficient list is

```text
c_02=-2, c_03=9, c_04=-12, c_05=5,
c_12=3, c_13=-18, c_14=27, c_15=-12,
c_23=8, c_24=-18, c_25=9, c_34=3, c_35=-2.
```

In particular the array is nonzero and has mixed signs. Reconstruction
nevertheless vanishes:

```text
f(b_0^2,...,b_5^2)
 =9 V(0,1,2)V(3,4,5)-V(0,1,3)V(2,4,5)
 =9*2*2-6*6=0.
```

By the terminal coherent-image theorem, `f` is the coherent array of
`Q=9 T_012 T_345-T_013 T_245` in `K_1(6)`, where
`T_abc=det((X_i^2,X_iY_i,Y_i^2))_(i=a,b,c)`. Direct evaluation gives the
same zero: at rows `(1,i)`, the two products of triangle determinants are
`4` and `36`. More generally, for rows `(N,i)` they are `4 N^6` and
`36 N^6`, respectively. Taking `N` large makes the slopes lie in an
arbitrarily short interval. The associated unit-circle points

```text
((N^2-i^2)/(N^2+i^2), 2Ni/(N^2+i^2))
```

all have equal radius, while the half-angle numerator rows `(N,i)` have
different norms `sqrt(N^2+i^2)`. Positive normalization preserves the
sign pattern, but the integer circuit equations here belong to the
unnormalized half-angle rows; equal-radius identities cannot be substituted
without tracking that rescaling.

This example shows why ordering and circuit signs alone are insufficient.
To force the actual coherent array to vanish for short numerical-zero
relations, one needs additional arithmetic input such as the full
core-profile and height exclusion in the actual-cut theorem. The example
does not have that core profile and does not challenge that theorem.

The standard-library checker
[`check_terminal_ordered_circuit_sign_countermodel.py`](check_terminal_ordered_circuit_sign_countermodel.py)
verifies the circuit equations, mixed signs, reconstruction, coherent
triangle realization, and the clustered-slope family exactly.
