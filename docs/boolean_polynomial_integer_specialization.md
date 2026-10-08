# From polynomial Boolean profiles to primitive integer families

Every rational polynomial profile of the following kind specializes to
an unbounded family satisfying the full integer profile. In particular,
integrality, parity and disjoint prime supports are not additional
obstructions once such a polynomial family has been found. This does
not assert the converse, or construct polynomial profiles in four or
more rows.

The later [Kummer-cover theorem](boolean_section_kummer_uniform_bound.md)
proves that these exact hypotheses force `m<=35`, independently of the
block degree, including repeated roots. Thus this specialization scheme
cannot produce a counterexample by constructing such profiles for
arbitrarily large `m`. The specialization theorem remains valid for
each existing profile; there is still no converse from arbitrary
integer configurations to polynomial families.

## 1. Hypotheses and conclusion

Fix `m>=2` and `e>=1`. For every nonempty subset `T` of `[m]`, let
`H_T in Q(i)[t]` be monic of degree `e`. Assume all the polynomials
`H_T` and `bar(H_T)` are pairwise coprime, including each polynomial
with its own conjugate. Suppose there are nonzero rational constants
`c_i,Y_i` such that

```text
P_i=c_i product_(T contains i) H_T = X_i+iY_i,
X_i in Q[t].                                            (1)
```

The row degree is `d=2^(m-1)e`. A constant Gaussian coefficient
in place of `c_i` would automatically be real, since the monic
product in (1) has a real leading coefficient.

There are a rational `t_0`, a positive integer `M`, nonzero rational
row scalars `rho_i`, and fixed conjugate-primitive Gaussian integers
`K_i` such that, for `t=t_0+Mh` with all sufficiently large integers
`h`,

```text
J_T=H_T(t)/H_T(t_0) in Z[i],       n_T=Norm(J_T)>1,
R_i=rho_i P_i(t)=K_i product_(T contains i)J_T in Z[i].   (2)
```

Every `J_T` and `R_i` is conjugate-primitive. The `n_T` are odd,
pairwise coprime, and coprime to every `Norm(K_i)`. All imaginary
coordinates `Im(R_i)` and determinant residuals in

```text
Im(bar(R_i)R_j)=tau_ij product_(T contains i,j)n_T       (3)
```

are fixed nonzero integers. Writing `w=2e log h`, we have

```text
log n_T=w+O(1),       log|R_i|=2^(m-2)w+O(1).           (4)
```

All constants may depend on the fixed polynomial family.

## 2. Choose a base point with the required parity

Choose `t_0=2^(-L)` with `L` sufficiently large. For any fixed
nonconstant `X_i`, its leading term has strictly smaller 2-adic
valuation than every lower term at this input, once `L` is large.
Moreover `v_2(X_i(t_0))<v_2(Y_i)`. Thus the primitive integer vector
proportional over `Q` to `(X_i(t_0),Y_i)` has an odd first coordinate
and an even second coordinate. Choose a nonzero rational `rho_i`
so that

```text
K_i=rho_i P_i(t_0)
```

is that primitive Gaussian integer. Its norm is odd, so it is
coprime in `Z[i]` to its conjugate. No polynomial `H_T` vanishes
at any real input, since a real zero would also be a zero of its
conjugate. All normalization denominators are therefore nonzero.

## 3. A finite Bezout certificate handles every prime

Put `J_T(u)=H_T(t_0+u)/H_T(t_0)` and
`n_T(u)=J_T(u)bar(J_T(u))`. Then `J_T(0)=n_T(0)=1`.
The original coprimality hypotheses imply

```text
gcd_Q[t](n_T,n_U)=1              for T!=U,
gcd_Q[t](Re(J_T),Im(J_T))=1     for every T.             (5)
```

For each pair in (5), choose a rational polynomial Bezout identity
and clear the coefficients of its two multipliers. This gives
integer polynomials `A,B` and a nonzero integer `D` such that
`Af+Bg=D`; the polynomials `f,g` themselves can remain rational.
At a specialization where all `J_T` coordinates are integers,
this identity shows that a common prime divisor of `f,g` must
divide `D`.

Let `S` consist of 2, the prime factors of these finitely many
nonzero integers `D`, and those of every `Norm(K_i)`. Choose `M`
so that every nonconstant coefficient of every `J_T(Mu)` belongs
to `q Z[i]`, where `q` is the product of the primes in `S`.
For example, a common multiple of `q` and all coefficient
denominators, multiplied by `q`, suffices. Then for every integer
`h`,

```text
J_T(Mh) in Z[i],       J_T(Mh)=1 mod p for every p in S.
```

Hence no prime in `S` divides `n_T(Mh)`. A common prime of two
norms, or of the two coordinates of one `J_T`, would have to lie
in `S` by its Bezout identity. This proves disjointness and ordinary
primitivity. Odd norms give conjugate-primitivity and exclude inert
prime factors. The norms also avoid every correcting norm.
Products in (2) are consequently conjugate-primitive as well.

The `n_T` have degree `2e` and positive leading coefficient, so
they exceed one for all sufficiently large positive `h`.

## 4. Every determinant residual is already fixed

The product of the shared norm polynomials divides
`Im(bar(P_i)P_j)` in `Q[t]`. Indeed its Gaussian factor divides
both rows, and its conjugate divides both conjugate rows. Its
degree is `2*2^(m-2)e=d`, whereas

```text
Im(bar(P_i)P_j)=Y_j X_i-Y_i X_j
```

has degree at most `d`. The determinant is nonzero: its vanishing
would make the two rows proportional, contradicting their distinct
incident polynomial factors. The quotient is therefore a nonzero
rational constant. After the normalizations in (2), evaluation at
`h=0` in the polynomial identity gives the particularly useful formula

```text
tau_ij=Im(bar(K_i)K_j) in Z\{0}.                        (6)
```

The imaginary coordinates are fixed and nonzero by (1), and integral
by (2). Their heights, the `K_i` heights, and all residual heights
are `O(1)` as `h` grows. Equation (4) follows directly from degrees.

## 5. Consequence for actual circles, with its quantifiers

These families give `m+1` actual integer points on a common centered
circle without requiring that it be the least possible circle.
Let `B=product_T J_T` and choose any nonzero Gaussian integer `E`
divisible by every `bar(K_i)`. Define

```text
z_0=bar(B) E,
z_i=(product_(T contains i)J_T)
    (product_(T does not contain i)bar(J_T)) E K_i/bar(K_i).
```

All points are Gaussian integers, have the same absolute value

```text
R=|E| product_T |J_T| = Theta(h^(e(2^m-1))),
```

and satisfy `z_i/z_0=R_i/bar(R_i)`. Since `Im(R_i)` is fixed
and its real part has degree `d`, all these ratios tend to 1 with
arguments `O(h^(-d))`. Equation (3) makes them pairwise distinct;
the nonzero imaginary coordinates distinguish them from `z_0`.
They fit on an arc of length

```text
O(R^alpha_m),
alpha_m=(2^(m-1)-1)/(2^m-1)=1/2-1/(2(2^m-1)).          (7)
```

For each fixed family this is `o(sqrt(R))`: for a given fixed `C>0`,
take `h` sufficiently large to fit its `m+1` points in length
`C sqrt(R)`. The constants in (7) may depend on `m`, since `h` is
chosen afterward. Thus arbitrary-size families under these hypotheses
would have disproved the requested uniform bound. The Kummer-cover
theorem now excludes that possibility by proving `m<=35`.
Existence is established here over `Q(i)` only for `m=2,3`; no
unbounded point count follows from those fixed sizes.

A [full four-row solution over real algebraic coefficients](four_row_real_polynomial_profile_certificate.md)
is now certified. Its `Q(i)`-rationality is unproved, so it does not
meet this theorem's hypotheses. A field norm does not in general
preserve a row's constant nonzero imaginary part.

Conversely, failure of the uniform integer bound has not been
shown to produce a rational polynomial family. The polynomial
nonexistence theorem for `m>=36` therefore still needs an additional
argument to settle the original problem.

The specific seven-factor construction and its least radius are in
[the three-row family](balanced_three_row_gaussian_polynomial_family.md)
and [its circle realization](three_row_balanced_family_least_circle.md).
The [exact checker](check_boolean_polynomial_integer_specialization.py)
tests the rational base point and primitive normalization after composing
that family with `t`, `t^2`, and `t^3`. It verifies the finite prime
certificate for each whole progression and eighteen actual profile and
circle realizations. The general assertion rests on the Bezout and
2-adic arguments above, rather than extrapolation from those checks.
