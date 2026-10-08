# Exact normalization under actual source reanchoring

This note computes every permitted source-anchor normalization of a rational
Möbius map. The result is an exact primewise distance formula and an explicit
obstruction to bounding normalized matrix height from the source radius
alone. It does not construct a radius-nonincreasing counterexample and does
not prove a uniform endpoint bound.

Let \(H_0=1,H_1,\dots,H_{m-1}\) be primitive integer Gaussian half-angle
columns of odd norms \(f_i\), for a primitive fixed-unit circle configuration
of least squared radius \(N\). Let
\[
 Mz=\alpha z+\beta\bar z,\qquad
 w=\beta/\alpha=v/u,
 \qquad u,v\in\mathbb Z[i],\quad (u,v)=1,
\]
where both coefficients are nonzero and \(Nw\ne1\). Write
\[
 Nu=U,\quad Nv=V,\qquad V/U=a/b,
 \quad a,b\in\mathbb Z_{>0},\quad(a,b)=1,
\]
and let
\[
 \epsilon=4/N\gcd(2,u-v)\in\{1,2,4\}.
\]

## 1. Actual reanchoring and exact Gaussian cancellation

After choosing source point \(j\) as the new anchor, represent the new phase
\(q_i/q_j\) by the integer-primitive normalization of \(H_i\bar H_j\).
Multiplication of this column by \(H_j\) is a nonzero real multiple of
\(H_i\). Therefore the corresponding map is
\(M\circ(z\mapsto H_jz)\), and its full output configuration is unchanged.
In particular its ratio of complex coefficients is
\[
 w_j=w\bar H_j/H_j.
\]
Put
\[
 g_j=\gcd(uH_j,v\bar H_j),\qquad G_j=Ng_j.
\]
Its exact reduced coefficient pair is
\[
 \boxed{u_j=uH_j/g_j,\qquad v_j=v\bar H_j/g_j.}
 \tag{1}
\]
Because each source column is primitive and has odd norm,
\((H_j,\bar H_j)=1\). Together with \((u,v)=1\), primewise comparison
gives the useful exact factorization
\[
 g_j=\gcd(u,\bar H_j)\gcd(v,H_j)
 \quad\text{up to a Gaussian unit}.
 \tag{2}
\]
Indeed, at any Gaussian prime the two pairs of valuations, those of
\((u,v)\) and of \((H_j,\bar H_j)\), each contain a zero. Thus
\(\min(v(u)+v(H_j),v(v)+v(\bar H_j))\) is precisely the sum of the two
cross minima in (2).

Consequently
\[
 Nu_j=U f_j/G_j,\qquad Nv_j=V f_j/G_j.
 \tag{3}
\]
There is a positive integer \(k_j\) such that
\[
 \boxed{Nu_j=bk_j,\qquad Nv_j=ak_j.}
 \tag{4}
\]
The factor \(k_j\), not the invariant rational ratio \(a/b\), is the part
that actual source reanchoring can change.

## 2. Exact local interval formula

For each odd split prime \(p=\pi\bar\pi\), fix one orientation and write
\[
 d_p=v_\pi(w),\qquad e_p=v_{\bar\pi}(w),\qquad
 t_{j,p}=v_\pi(H_j)-v_{\bar\pi}(H_j).
\]
Primitivity implies that the two source Gaussian valuations are
\((t_{j,p})_+\) and \((-t_{j,p})_+\). Under reanchoring,
\[
 v_\pi(w_j)=d_p-t_{j,p},\qquad
 v_{\bar\pi}(w_j)=e_p+t_{j,p}.
\]
Reduction of the numerator and denominator at the two Gaussian primes gives
\[
 v_p(Nu_j)=(t_{j,p}-d_p)_++(-t_{j,p}-e_p)_+,
\]
\[
 v_p(Nv_j)=(d_p-t_{j,p})_++(e_p+t_{j,p})_+.
\]
Their difference is the invariant \(d_p+e_p=v_p(a/b)\). Their minimum is
the exponent of \(k_j\), and therefore
\[
 \boxed{
 v_p(k_j)=\operatorname{dist}(t_{j,p},I_p),\qquad
 I_p=[\min(d_p,-e_p),\max(d_p,-e_p)].}
 \tag{5}
\]
Equivalently the distance equals
\[
 \frac{|t_{j,p}-d_p|+|t_{j,p}+e_p|-|d_p+e_p|}{2}.
\]
All endpoints and source coordinates are integers, so this is an integer.
The formula holds for arbitrary nested source valuations and requires no
common ordering of the rows at different primes.

Inert primes cannot divide a primitive source norm. They also cannot divide
both \(Nu_j\) and \(Nv_j\), because \(u_j,v_j\) are Gaussian coprime.
The same assertion holds at two, whose unique Gaussian prime cannot divide
both reduced coefficients. Hence
\[
 \boxed{k_j=\prod_{p\equiv1\pmod4}
 p^{\operatorname{dist}(t_{j,p},I_p)}.}
 \tag{6}
\]
Thus the permitted anchors are a finite collection of vectors of source
valuations, and their matrix heights are determined by weighted sums of
distances to a fixed collection of intervals.

## 3. The parity factor, trace, discriminant, and good determinant

The parity factor \(\epsilon\) is unchanged by actual reanchoring. Indeed,
\(H_j\equiv\bar H_j\pmod{2\mathbb Z[i]}\), and both \(H_j\) and \(g_j\)
are units modulo two. Therefore
\[
 u_j-v_j\equiv H_jg_j^{-1}(u-v)\pmod{2\mathbb Z[i]}.
\]
This preserves the valuation at \(1+i\) truncated at two, exactly the
information defining \(\epsilon\).

Applying the exact output-conformal normalization theorem to (1) now gives
the minimal integral representative at anchor \(j\):
\[
 \boxed{
 T_j=\frac\epsilon2(a+b)k_j,\qquad
 K_j=\epsilon^2ab\,k_j^2,\qquad
 \delta_j=\frac\epsilon4(b-a)k_j.}
 \tag{7}
\]
Set
\[
 T_c=\epsilon(a+b)/2,\qquad K_c=\epsilon^2ab,\qquad
 c=\epsilon(b-a)/4.
\]
These are integers. For \(c\), integrality of \(\delta_j\) and oddness of
\(k_j\) show that the possible denominator four divides the numerator.
The trace assertion follows in the same way. Also \(c\ne0\).

Define the good determinant part of the coefficient core by
\[
 D_c=\prod_{\substack{p\mid c\\p\nmid K_c}}p^{v_p(c)}.
\]
Since any prime dividing \(k_j\) also divides \(K_j\), the exact good
determinant factor is
\[
 \boxed{
 D_j=\prod_{\substack{p\mid c\\p\nmid K_c\,k_j}}p^{v_p(c)}
 =\frac{D_c}{\displaystyle\prod_{p\mid k_j}p^{v_p(D_c)}}.}
 \tag{8}
\]
The product in the denominator only uses primes dividing \(D_c\).
In particular, a good prime contribution is retained at anchor \(j\)
exactly when its source valuation lies inside \(I_p\). The two-part and
every inert-prime part of \(D_c\) are constant across anchors, because
\(k_j\) contains only split primes.

The determinant-free compressed matrix loss is exactly
\[
 T_j^mK_j^{m-2}=T_c^mK_c^{m-2}k_j^{3m-4}.
 \tag{9}
\]
The joint conductor refinement can multiply the resulting lower bound by
\[
 \max\{1,D_j^{2m-3}/(T_c^2k_j^2)\}.
 \tag{10}
\]
Thus both bounds can be compared across actual anchors using (5), (8),
and no further hidden coordinate normalization.

These minimal representatives also optimize the conductor-enhanced bound
within each output-conformal class. The normalization ideal shows that every
other integral representative is \(L_\lambda M_j\), where \(L_\lambda\)
is multiplication by a nonzero Gaussian integer \(\lambda\). Put
\(s=N\lambda\ge1\). Its invariants are
\[
 T=sT_j,\qquad K=s^2K_j,\qquad\delta=s\delta_j.
\]
Every prime dividing \(s\) now divides \(K\) and contributes nothing to
the good determinant part. At the other primes the determinant valuations
are unchanged. Hence \(D(L_\lambda M_j)\mid D_j\). Both the base lower
bound and the conductor-enhanced lower bound are therefore maximized at
the displayed minimal representative.

## 4. What source valuation widths do control

The primitive source radius has the exact factorization
\[
 N=\prod_{p\equiv1\pmod4}p^{W_p},\qquad
 W_p=\max_j t_{j,p}-\min_j t_{j,p}.
\]
For two anchors \(i,j\), let \(f_{ij}\) be the norm of the primitive
half-angle numerator of \(q_i/q_j\). Then
\[
 f_{ij}=\prod_p p^{|t_{i,p}-t_{j,p}|}\mid N.
\]
Distance to an interval is one-Lipschitz, so (5) yields exact divisibilities
\[
 \boxed{k_i\mid k_jf_{ij},\qquad k_j\mid k_if_{ij}.}
 \tag{11}
\]
Consequently
\[
 1/N\le k_i/k_j\le N,
 \qquad \min_j k_j\ge k_0/N.
 \tag{12}
\]
If every \(I_p\) intersects the full source interval
\([\min_jt_{j,p},\max_jt_{j,p}]\), then (5) additionally gives
\(k_j\mid N\) for every anchor. This condition includes primes not dividing
\(N\): their intervals must contain zero.

Without that condition, source widths say nothing about the distance of the
map intervals from the source. At every split prime absent from \(N\), all
source valuations are zero and the factor
\(p^{\operatorname{dist}(0,I_p)}\) is shared by every \(k_j\). Actual
reanchoring cannot remove any of it.

## 5. A uniform obstruction using a prime outside the source radius

For any source configuration as above, choose a split prime
\(p=\pi\bar\pi\) not dividing \(N\), and any integer \(L\ge1\). Consider
the rational map with reduced coefficient pair
\[
 u=2\bar\pi^L,\qquad v=\pi^L,
 \qquad w=\frac{\pi^L}{2\bar\pi^L}.
\]
Its norm ratio is the fixed number \(1/4\), independently of the source,
\(p\), and \(L\). Because \(p\nmid N\), no \(H_j\) is divisible by
\(\pi\) or \(\bar\pi\); odd source norms also exclude the prime over two.
Thus (2) gives \(g_j=1\) for every anchor. Here
\[
 a=1,\quad b=4,\quad\epsilon=4,\qquad
 \boxed{k_j=p^Lf_j.}
 \tag{13}
\]
In particular \(k_0=p^L\) and no permitted anchor has smaller \(k_j\).
All normalized matrix invariants are explicit:
\[
 \boxed{
 T_j=10p^Lf_j,\quad K_j=64p^{2L}f_j^2,\quad
 \delta_j=3p^Lf_j,\quad D_j=3.}
 \tag{14}
\]
The last assertion uses that three is inert and hence divides neither
\(p\) nor a primitive source norm \(f_j\). The conductor gain in (10)
therefore supplies no growing cancellation as \(L\) increases.

For fixed source radius, the minimum normalized trace over all actual
anchors tends to infinity with \(L\). This refutes any upper bound on that
minimum depending only on the source, even after fixing the full rational
conformal double-orbit invariant \(Nw=1/4\). It is an exact obstruction to
an unrestricted anchor-height argument, not an example satisfying \(N'\le N\).
An argument using actual reanchoring to settle the endpoint problem must
also use a restriction from the output configuration or the proposed
radius descent.
