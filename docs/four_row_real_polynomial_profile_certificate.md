# A certified full four-row real polynomial profile

There is an exact nondegenerate full four-row polynomial profile with all
fifteen linear factors over the complex algebraic numbers.
Each row product has a nonzero constant imaginary part, and each pair bracket
is a nonzero constant times its four shared quadratic norm factors.  The
coefficients are not proved rational.  Consequently this result rules out an
obstruction based only on the real polynomial equations and positivity; it
does **not** supply an integer profile or meet the hypotheses of
[the rational polynomial specialization theorem](boolean_polynomial_integer_specialization.md).

The proof uses a small exact rational contraction certificate, stored in
[four_row_real_polynomial_certificate.json](four_row_real_polynomial_certificate.json)
and verified by the dependency-free checker
[check_four_row_real_polynomial_certificate.py](check_four_row_real_polynomial_certificate.py).
The floating-point search is only the discovery step.

## 1. The equations and the certified root

Index the nonempty subsets of \([4]\) by masks \(1,\ldots,15\), with bit
\(i-1\) indicating membership of row \(i\).  Set \(z_{15}=i\).  The remaining
fourteen roots give twenty-eight real coordinates, ordered first by their
real parts and then by their imaginary parts.  Define the polynomial map
\(F:\mathbb R^{28}\longrightarrow\mathbb R^{28}\) by

\[
 F_{ik}(x)=\operatorname{Im}\sum_{T\ni i}z_T^k,
 \qquad 1\leq i\leq4,\quad1\leq k\leq7.
 \tag{1}
\]

The certificate stores a rational point \(x_0\), a rational matrix \(A\),
and the radius \(\rho=10^{-30}\).  It proves that the coordinate cube
\(\|x-x_0\|_\infty\leq\rho\) contains a unique zero \(x_*\) of (1).
For orientation, its root coordinates are approximately:

| Mask | Subset | Real part | Imaginary part |
|---:|:---|---:|---:|
| 1 | 1 | 0.4944488062001512 | 0.2308248408334662 |
| 2 | 2 | 0.4188656378298891 | 0.4492215096953632 |
| 3 | 12 | 0.2044399393979019 | 1.0074295723585633 |
| 4 | 3 | 0.1826695291154774 | 1.0268356347459103 |
| 5 | 13 | 0.1761569489398476 | 0.7824010045215722 |
| 6 | 23 | 0.3212402885633954 | -0.5807115383499899 |
| 7 | 123 | 0.0354399200976209 | -0.9943811171308613 |
| 8 | 4 | 0.2712857344955086 | 1.0847400160095055 |
| 9 | 14 | 0.1920072766720538 | -0.9447886210580642 |
| 10 | 24 | 0.3328492893790878 | 0.3093376984106465 |
| 11 | 124 | 0.5216144998258322 | -0.2151451095754565 |
| 12 | 34 | 0.3100773132537340 | 0.6079476015708537 |
| 13 | 134 | 0.1489545794959850 | -0.8663405699492197 |
| 14 | 234 | 0.1671262311948373 | -0.9757510154082653 |
| 15 | 1234 | 0 | 1 |

The table is not the certificate: the JSON file retains sixty-five decimal
places and the full rational approximate inverse.

## 2. Exact contraction proof

All the following inequalities are checked using Python `Fraction`, with
each stored decimal interpreted as an exact rational number:

\[
\begin{split}
 \alpha&=\|AF(x_0)\|_\infty<10^{-60},\\
 \beta&=\|I-AJ(x_0)\|_\infty<10^{-59},\\
 \|A\|_\infty&<2700,\\
 |\operatorname{Re}z_T(x_0)|+|\operatorname{Im}z_T(x_0)|+2\rho&<2.
\end{split}
\tag{2}
\]

Here \(J=DF\).  Every root throughout the cube has absolute value less
than two.  For a fixed equation of power \(k\), there are seven variable
roots, two coordinates per root, and two possibly nonzero second
derivatives per coordinate.  Each second derivative has absolute value at
most \(k(k-1)2^{k-2}\).  The mean value theorem therefore gives

\[
 \|J(x)-J(x_0)\|_\infty
 \leq 28\max_{1\leq k\leq7}\bigl(k(k-1)2^{k-2}\bigr)\rho
 =37632\rho.
 \tag{3}
\]

Thus \(T(x)=x-AF(x)\) has Lipschitz constant at most

\[
 q=\beta+\|A\|_\infty37632\rho<1.014\times10^{-22}<\tfrac12.
 \tag{4}
\]

Also \(\alpha+q\rho<\rho\), so \(T\) maps the closed cube into itself.
Banach's contraction theorem gives a unique fixed point.  Since
\(\beta<1\), the square matrix \(AJ(x_0)\), and hence \(A\), is invertible;
the fixed point is a zero of \(F\).  Inequality (4) also makes
\(AJ(x_*)\) invertible.  The zero is nonsingular for this polynomial system
over \(\mathbb Q\), hence its coordinates are real algebraic numbers.

The measured rational quantities, converted to floating point only after
the exact comparisons, are

\[
 \alpha\approx4.648886\times10^{-66},\quad
 \beta\approx1.497284\times10^{-62},\quad
 \|A\|_\infty\approx2693.122207.
\]

## 3. Nondegeneracy and the polynomial profile

The checker proves that every two members of
\(\{z_T(x_*),\overline{z_T(x_*)}:T\ne\varnothing\}\) have coordinate
\(\ell^\infty\) separation greater than \(1/50\).  It also proves
\(|\operatorname{Im}z_T(x_*)|>1/5\) for every \(T\).  In particular all
thirty roots are distinct, and no quadratic norm factor has a real zero.

Set

\[
 H_T(t)=t-z_T(x_*),\qquad
 n_T(t)=H_T(t)\overline{H_T(t)},\qquad
 P_i(t)=\prod_{T\ni i}H_T(t).
 \tag{5}
\]

Newton's identities applied to the eight roots in each row show from (1)
that all coefficients of \(P_i\), except possibly its constant term, are
real.  Hence

\[
 P_i(t)=X_i(t)+iY_i,
 \qquad X_i\in\mathbb R_{\rm alg}[t]\text{ monic of degree }8,
 \quad Y_i\in\mathbb R_{\rm alg}.
 \tag{6}
\]

The four imaginary constants are approximately

\[
 (Y_1,Y_2,Y_3,Y_4)=
 (0.0064814521683495,\ 0.0124041227613239,
  0.0038234779864941,\ 0.0295954655875636).
 \tag{7}
\]

The constant term is a product of eight roots, one fixed.  On the cube,
its change has absolute value at most
\(7(2\rho)2^7=1792\rho\).  Exact rational checks of the products at
\(x_0\) consequently prove

\[
 \min_iY_i>3/1000,\qquad
 \min_{i\ne j}|Y_i-Y_j|>1/500.
 \tag{8}
\]

For \(i\ne j\), the real polynomial

\[
 \Delta_{ij}=\operatorname{Im}(\overline{P_i}P_j)
             =X_iY_j-Y_iX_j
\]

is divisible by \(\prod_{T\supset\{i,j\}}n_T\).  Both polynomials have
degree eight; comparison of their leading coefficients gives the exact
identity

\[
 \boxed{\displaystyle
 \Delta_{ij}(t)=(Y_j-Y_i)
       \prod_{T\supset\{i,j\}}n_T(t).}
 \tag{9}
\]

Every residual in (9) is nonzero by (8).  All fifteen \(n_T\) are monic,
positive, pairwise coprime quadratics over \(\mathbb R\), and for real
\(t\longrightarrow+\infty\) they satisfy
\(\log n_T(t)=2\log t+o(1)\).  The rank-two bracket relations and all norm
identities hold automatically because (5) supplies the actual row products.
They are not being imposed independently of the roots.

## 4. Arithmetic scope

This is an exact profile over a finite extension of \(\mathbb Q(i)\),
with real parameter \(t\).  The proof does not establish
\(z_T\in\mathbb Q(i)\), integer values at any parameter, disjoint rational
prime supports, or conjugate-primitive Gaussian integer rows.  None of
these arithmetic conclusions follows from contraction certification.

Thus a proposed four-row exclusion must use arithmetic information absent
from a general real algebraic polynomial profile.  In particular, the joint
norm circuit, bracket relations, constant nonzero imaginary parts, and
strict root/conjugate separation do admit simultaneous nondegenerate
solutions over the real algebraic numbers.  Whether there is a rational
Gaussian profile remains open here.

Clearing denominators in a number field does not make its coefficients
rational.  Nor does taking a field norm preserve a constant imaginary
part: for \(d\geq1\), let
\(A=(t-\sqrt2)^d\), \(B=(t+\sqrt2)^d\).  Then \(A+i\) has constant
imaginary part, but its relative norm from
\(\mathbb Q(\sqrt2,i)\) to \(\mathbb Q(i)\) is
\[
 (A+i)(B+i)=(t^2-2)^d-1+i(A+B),
\]
whose imaginary part has degree \(d\).  A separate arithmetic construction
would therefore be needed to apply the rational specialization theorem.

The discovery search used NumPy with an analytic-Jacobian damped
least-squares iteration because SciPy was unavailable.  It ran for
approximately 53 seconds before locating this candidate; Decimal Newton
refinement then supplied the rational certificate.  Reproduction of the
proof needs only:

```sh
python3 docs/check_four_row_real_polynomial_certificate.py
```
