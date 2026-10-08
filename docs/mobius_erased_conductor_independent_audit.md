# Independent audit of clipped-width radius retention

This note independently checks the global retention argument in
[the erased-conductor budget](mobius_erased_conductor_budget.md). The
prime-layer registration, optimized elementary inequality, arbitrary target
parity factor, and final exponent are all valid. In particular, for
\(m\ge5\), a primitive odd fixed-unit source and its full primitive image,
both with endpoint constants at most two, satisfy
\[
 \boxed{N'\ge\frac18N^{\,1-4/m}.}                         \tag{1}
\]
No target unit-class extraction is used.

The exact formula for the cut minimum is
\[
 c_m=
 \frac{k[(m-1)(m-2k-1)+1]}{4(m-2k-1)},\qquad
 k=\lfloor m/4\rfloor.                                  \tag{2}
\]
This corrects an earlier informal candidate for \(c_m\); the retention ratio
stated in the main note is unchanged.

## 1. Registering source cuts and clipped cuts

At an odd split prime \(p\), let the signed source valuations be
\(t_{1,p},\ldots,t_{m,p}\). Write
\[
 w_p=\max_i t_{i,p}-\min_i t_{i,p}.
\]
Every integer layer between the attained minimum and maximum defines a
nontrivial cut. If one side has \(n\) rows, that layer contributes
\[
 n(m-n)
\]
to the exponent of
\[
 P=\prod_{i<j}f_{ij}.
\]
Consequently, with \(W=\log N\),
\[
\begin{aligned}
 W&=\sum_{\text{source layers}}\log p,\\
 \log P&=\sum_{\text{source layers}}n(m-n)\log p,\\
 D:=\frac{m^2}{4}W-\log P
 &=\sum_{\text{source layers}}(n-m/2)^2\log p.            \tag{3}
\end{aligned}
\]
This is an equality for both even and odd \(m\). For odd \(m\), the
unavoidable \(1/4\) at a most-balanced cut is correctly retained in \(D\).

For every pair, the source endpoint hypothesis \(C\le2\) gives
\[
 f_{ij}\ge\frac4{C^2}\sqrt N\ge\sqrt N.
\]
There are \(\binom m2\) pairs, so
\[
 \log P\ge\frac{m(m-1)}4W,\qquad D\le\frac m4W.           \tag{4}
\]
No fairness or binary-prime assumption enters (3)--(4).

Let \(I_p\) be the integral coefficient interval and set
\[
 z_{i,p}=\operatorname{clip}_{I_p}(t_{i,p}),\qquad
 w_{p,\rm clip}=\max_i z_{i,p}-\min_i z_{i,p}.
\]
The source layers counted by
\(w_p-w_{p,\rm clip}\) are exactly the layers strictly beyond one endpoint
of \(I_p\). A row lying on an endpoint has distance zero and belongs to no
erased layer. If \(I_p\) and the source range are disjoint, every source
layer is erased; the additional gap between the intervals is nonnegative
overlap which is not needed below. If the source width is zero, including
at a prime outside \(N\), there is no source layer and hence no contribution
to either side of the retention inequality.

Define
\[
 W_{\rm clip}=\sum_p w_{p,\rm clip}\log p,\qquad
 E_0=W-W_{\rm clip}.
\]
The clipped-conductor divisibility
\(N_{\rm clip}\mid\gcd(N,N')\) gives
\[
 W_{\rm clip}\le W',\qquad W'=\log N',\qquad
 E_0\ge W-W'.                                            \tag{5}
\]

## 2. Every erased layer is registered in target residues

Consider one erased layer and orient its cut so that its \(n\) rows are the
rows farther beyond the relevant endpoint of \(I_p\). Every pair among
these \(n\) rows lies strictly on the same side of \(I_p\). The same-side
overlap theorem therefore contributes at least one \(p\)-valuation to that
target primitive pair residue. Summing the unit layers in one tail gives
the exact minimum-distance overlap; the other tail is disjoint, so there is
no double counting. Summing over primes gives
\[
 \log B'\ge
 \sum_{\text{erased layers}}\binom n2\log p,
 \qquad B'=\prod_{i<j}b'_{ij}.                            \tag{6}
\]
The side used for \(n\) in (6) can be the complement of the side used in
(3), but
\[
 (n-m/2)^2=((m-n)-m/2)^2.
\]
Thus the same deficit charge applies. Endpoint cancellation causes no gap:
positive overlap is used only for rows strictly outside the closed
coefficient interval.

For \(A\ge0\), put
\[
 c(A)=\min_{1\le n\le m-1}
 \left\{\binom n2+A(n-m/2)^2\right\}.                    \tag{7}
\]
Applying (7) to every erased layer and allowing the nonnegative deficit
from retained layers to remain on the right gives
\[
 \boxed{c(A)E_0\le\log B'+AD.}                           \tag{8}
\]
This also proves that boundary-zero layers and external primes require no
separate error term.

## 3. Exact optimization

Let
\[
 k=\lfloor m/4\rfloor,\qquad
 L=m-2k-1,\qquad A_m=\frac{k}{L}.
\]
For
\[
 F_A(n)=\binom n2+A(n-m/2)^2,
\]
the adjacent difference is
\[
 F_A(n+1)-F_A(n)=n+A(2n+1-m).                            \tag{9}
\]
At \(A=A_m\), the difference in (9) vanishes at \(n=k\). Since it is
strictly increasing in \(n\), the integer minimum is attained at
\(k,k+1\). Therefore
\[
\begin{aligned}
 c_m=c(A_m)
 &=\binom k2+\frac{k}{L}(k-m/2)^2\\
 &=\frac{k[(m-1)L+1]}{4L},                               \tag{10}\\
 G_m&=c_m+\frac m8
 =\frac{(m-1)[m+2k(m-2k-2)]}{8(m-2k-1)}.
\end{aligned}
\]

Suppose
\[
 \log B'\le\frac m8W'+\eta.                              \tag{11}
\]
Using (4)--(5) in (8) gives
\[
 (c_m-A_mm/4)W\le(c_m+m/8)W'+\eta.
\]
Direct substitution into this ratio yields
\[
 \boxed{
 W'\ge\rho_mW-\frac{\eta}{G_m},\qquad
 \rho_m=
 \frac{2k(m-2k-2)}
 {m+2k(m-2k-2)}.}                                       \tag{12}
\]

This is also the optimum over the one-parameter family (8). The lower
envelope in (7) is piecewise affine in \(A\), and on each piece the ratio
\((c(A)-Am/4)/(c(A)+m/8)\) has derivative of constant sign. It is therefore
enough to inspect adjacent-line breakpoints
\[
 A=\frac{j}{m-2j-1}.
\]
The resulting numerator is governed by \(j(m-2-2j)\), whose maximizing
integer is \(j=\lfloor m/4\rfloor\), with the expected tie when
\(m\equiv0\pmod4\).

To verify the simple lower bound, write \(m=4k+s\),
\(s\in\{0,1,2,3\}\). The inequality
\(\rho_m\ge1-4/m\) is equivalent to
\[
 4m+8k(m-2k-2)-m^2=s(4-s)\ge0.                           \tag{13}
\]
For \(m\ge5\), \(\rho_m>0\).

## 4. Full target parity budget

No fixed target unit class is needed for (11). Let \(r\) of the \(m\)
primitive target columns have two odd coordinates. At every odd split
prime, target pair norms obey the layer-cake Plotkin bound. At two, a
primitive target pair norm has exponent one exactly when its endpoints lie
on opposite sides of the parity cut, giving \(r(m-r)\) such pairs. The
individual target chord inequalities therefore give
\[
 \log B'\le
 \frac m8W'
 +\binom m2\log(C'/2)
 +\frac{r(m-r)}2\log2.                                  \tag{14}
\]
This uses all target rows and changes neither \(N'\) nor the target
configuration.

For \(C'\le2\), (14) permits
\[
 \eta\le\frac{m^2}{8}\log2.
\]
Also
\[
 G_m\ge\frac{m^2}{24}\qquad(m\ge5),
\]
because the bracket \(m+2k(m-2k-2)\) in (10) is at least
\(m^2/4\), while \(m-2k-1\le(m+1)/2\); hence
\[
 G_m\ge\frac{m^2(m-1)}{16(m+1)}\ge\frac{m^2}{24}.
\]
so \(\eta/G_m\le3\log2\). Equations (12)--(13) now give
\[
 W'\ge(1-4/m)W-3\log2.
\]
Exponentiation proves (1).

If the target parity factor is absent, one may use \(\eta=0\) and retain
the sharper exact conclusion
\[
 N'\ge N^{\rho_m}.
\]

## 5. Audit conclusion

The retention argument is rigorous for arbitrary nested source
prime powers and arbitrary locations and lengths of the coefficient
intervals. Clipped layers, erased layers, interval-boundary rows, disjoint
intervals, and primes outside the source conductor are all registered with
the correct sign. The only correction needed in the preliminary statement
was the closed form (2) for \(c_m\); it does not alter \(\rho_m\).

The result is a near-unit comparison of the two radii. It does not prohibit
\(N'=N\), give a radius-independent point bound, or improve the known
qualitative exponent from the separate projective-orbit argument.
