# Effective reciprocal-square tails for critical real-rooted supports

The non-effective compactness argument in the [uniform reciprocal-tail
note](critical_moment_uniform_reciprocal_tail.md) can be made quantitative.
There is a universal linear lower bound on the ordered magnitudes, and hence
an \(O(1/K)\) reciprocal-square tail.

Let \(m=2s+1\), \(s\geq1\), and let \(z_1,\ldots,z_m\) be nonzero real
numbers with distinct absolute values such that

\[
\sum_j z_j^{2r+1}=0,\qquad 0\leq r<s.
\]

Write \(a_1<a_2<\cdots<a_m\) for their absolute values. Then, with

\[
\lambda=\frac{\log 2}{4},
\]

one has

\[
 a_j\geq \lambda j\,a_1\quad(1\leq j\leq m),                 \tag{1}
\]

and, for every \(K\geq1\),

\[
 \sum_{j>K}a_j^{-2}
 \leq \frac{1}{\lambda^2K}\,a_1^{-2}
 =\frac{16}{(\log2)^2K}\,a_1^{-2}.                           \tag{2}
\]

The constants are deliberately simple rather than optimized.

## 1. The normalized odd polynomial

Set \(y_j=1/z_j\), and first scale so that \(a_1=1\). Define

\[
 E(t)=\prod_{j=1}^m(1-y_jt).
\]

Newton's identities applied to the given odd power sums show that the root
polynomial, after division by its nonzero constant term, has the form

\[
 E(t)=1+F(t),
\]

where \(F\) is an odd real polynomial of degree \(m\). Equivalently,

\[
 E(t)+E(-t)=2.                                           \tag{3}
\]

The roots of \(E\) are the \(z_j\), so \(F(z_j)=-1\). The roots of
\(E(-t)=1-F(t)\) are \(-z_j\), so \(F(t)=1\) also has \(m\) real simple
roots.

The coefficient of \(t\) in \(E\) is \(-p_1\), where

\[
 p_1=\sum_j y_j,
 \qquad F'(0)=-p_1.                                      \tag{4}
\]

The critical paired sign order from the [odd root-order note](critical_odd_moment_root_order.md)
gives, after ordering the \(y_j\) by decreasing absolute value,
the signs \(++--++--\cdots\), up to an overall sign. Since
\(|y_1|=1\), the alternating pair sums imply

\[
 |p_1|<2.                                                 \tag{5}
\]

The strict inequality also follows from \(p_1^2=p_2>0\) together with that
paired order; only the displayed upper bound is used below.

## 2. Every intermediate fiber is real

The derivative \(F'=E'\) has degree \(m-1=2s\). Rolle's theorem gives one
real derivative root in each interval between consecutive roots of \(E\), so
these are all the roots of \(F'\). Thus the real line is divided into
\(m\) monotonicity intervals for \(F\).

The polynomial \(F+1=E\) has one root in every such interval: it has \(m\)
real roots and a strictly monotone function has at most one root on each
interval. The same argument applies to \(F-1=-E(-t)\), whose \(m\) roots
are also real. Hence every monotonicity interval contains one point where
\(F=-1\) and one point where \(F=1\). Its image therefore contains the
whole interval \([-1,1]\). It follows that, for every \(u\in[-1,1]\),

\[
 F(t)=u
\quad\text{has exactly \(m\) real roots, one in each interval}. \tag{6}
\]

In particular, if \(z\) is nonreal then \(F(z)\notin[-1,1]\). Thus

\[
 F(\mathbb H)\subset\Omega:=\mathbb C\setminus[-1,1],               \tag{7}
\]

where \(\mathbb H\) is the upper half-plane.

## 3. Green function and boundary Harnack estimate

On \(\Omega\), choose the branch of \(\sqrt{w^2-1}\) for which
\(w+\sqrt{w^2-1}\sim2w\) at infinity, and put

\[
 g(w)=\log\left|w+\sqrt{w^2-1}\right|.
\]

This is the positive harmonic Green function of the slit plane, vanishing
continuously on \([-1,1]\). By (7),

\[
 u(z)=g(F(z))
\]

is positive harmonic on \(\mathbb H\).

For completeness, the needed sharp one-ray Harnack estimate follows from
the disk Harnack inequality. Fix \(0<\varepsilon<T\), and map \(\mathbb H\)
to the unit disk by

\[
 \psi_\varepsilon(z)=\frac{z-i\varepsilon}{z+i\varepsilon}.
\]

The point \(i\varepsilon\) maps to \(0\), while \(iT\) maps to

\[
 r=\frac{T-\varepsilon}{T+\varepsilon}\in(0,1).
\]

Applying \(v(r)\leq(1+r)(1-r)^{-1}v(0)\) to the positive harmonic function
\(v=u\circ\psi_\varepsilon^{-1}\) gives

\[
 u(iT)\leq \frac{T}{\varepsilon}u(i\varepsilon).              \tag{8}
\]

Since \(F\) is odd with real coefficients, \(F(i\varepsilon)\) is purely
imaginary. For real \(q\), the chosen Green function satisfies

\[
 g(iq)=\operatorname{arsinh}|q|.

\]

Therefore, using \(F(i\varepsilon)=iF'(0)\varepsilon+O(\varepsilon^3)\),

\[
 \lim_{\varepsilon\downarrow0}\frac{u(i\varepsilon)}{\varepsilon}
 =|F'(0)|=|p_1|<2.                                      \tag{9}
\]

Letting \(\varepsilon\downarrow0\) in (8) yields the degree-free estimates

\[
 u(iT)\leq |p_1|T<2T,\qquad
 |F(iT)|\leq\sinh(|p_1|T)\qquad(T>0).                    \tag{10}
\]

## 4. Counting and the tail

At \(iT\), equation (3) and oddness give

\[
 |E(iT)|^2
 =\prod_j(1+T^2y_j^2)
 =1+|F(iT)|^2.                                           \tag{11}
\]

If \(U=u(iT)=\operatorname{arsinh}|F(iT)|\), then the right side of (11)
is \(\cosh^2U\). Consequently, by (10),

\[
 \sum_j\log(1+T^2y_j^2)
 =2\log\cosh U
 \leq2U
 <4T.                                                     \tag{12}
\]

Let

\[
 N(T)=\#\{j:|y_j|\geq T^{-1}\}
     =\#\{j:a_j\leq T\}.
\]

Each counted term in (12) contributes at least \(\log2\), so

\[
 N(T)<\frac{4T}{\log2}.                                  \tag{13}
\]

At \(T=a_j\), distinct absolute values give \(N(a_j)=j\), proving (1).
The equivalent reciprocal estimate is \(|y_j|<4/(j\log2)\), and hence

\[
 \sum_{j>K}a_j^{-2}=\sum_{j>K}y_j^2
 \leq\frac{16}{(\log2)^2}\sum_{j>K}j^{-2}
 \leq\frac{16}{(\log2)^2K},

\]

which is (2) in the normalization \(a_1=1\). Rescaling restores the factor
\(a_1^{-2}\).

## 5. Consequence for the sixteen cross supports

Suppose the eight-row application is split into two groups of four and all
sixteen cross supports satisfy the corresponding odd critical moment
hypotheses, with globally distinct positive magnitudes. Discard columns
constant on all eight rows, calling the remaining columns active.
Let \(a_*\) be the smallest magnitude among the active columns,
and let \(N_D(X)\) count entries of a cross support \(D\) whose magnitudes
are at most \(X\). Applying (13) relative to the minimum magnitude in \(D\)
gives, for every \(T>0\),

\[
 N_D(Ta_*)
 \leq \frac{4}{\log2}\,\frac{Ta_*}{\min_{h\in D}a_h}
 \leq \frac{4T}{\log2}.                                  \tag{14}
\]

Each active column belongs to at least four of the sixteen cross supports:
if its positive-sign counts in the two groups are \(r,t\), its incidence is
\(r(4-t)+(4-r)t\), whose minimum over nonconstant columns is \(4\). Therefore,
if \(N_{\rm act}(Ta_*)\) counts active columns up to \(Ta_*\),

\[
 4N_{\rm act}(Ta_*)
 \leq\sum_DN_D(Ta_*)
 \leq\frac{64T}{\log2},
\]

so

\[
 N_{\rm act}(Ta_*)\leq\frac{16T}{\log2}.                  \tag{15}
\]

If \(A_j\) are the active magnitudes in increasing order, then

\[
 A_j\geq \frac{\log2}{16}j\,a_* ,\qquad
 \sum_{j>K}A_j^{-2}
 \leq \frac{256}{(\log2)^2K}\,a_*^{-2}.                  \tag{16}
\]

This is a statement about the critical supports and their active columns.
It does not by itself give a point-count bound in an original circle radius:
there is no degree bound, and the relation between these support magnitudes
and a lattice-circle radius is an additional problem. Chebyshev families of
unbounded degree are compatible with (1): their normalized support extends
linearly with the degree. The bound controls the number of support magnitudes
in each fixed normalized interval, but does not bound the total degree.
