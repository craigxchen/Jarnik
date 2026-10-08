# A joint row-content and phase-overlap conductor bound

This note proves a determinant gain by accounting jointly for ordinary row
contents and Gaussian phase overlaps. Both an all-pairs bound and a
compressed bound with exponent \(m-2\) are proved. Taking the larger
of the compressed bound and the earlier bound gives a uniformly valid
improvement factor, which need not exceed one strictly. These estimates
do not supply the missing uniform endpoint obstruction.

Let \(M\) be an integral real \(2\times2\) matrix of nonzero determinant
\(\delta\), and let
\[
 T=\operatorname{tr}(M^TM),\qquad K=T^2-4\delta^2>0.
\]
Let \(H_1,\dots,H_m\), \(m\ge3\), be pairwise nonparallel primitive
integer columns of odd squared norms \(f_i\). Put
\[
 r_i=\gcd_{\mathbb Z}(MH_i),\quad B_i=MH_i/r_i,\quad
 b_{ij}=\frac{|\det(H_i,H_j)|}{N\gcd_{\mathbb Z[i]}(H_i,H_j)}.
\]
Here each \(r_i\) is positive. The old Gaussian gcd norm divides the
integer determinant, so \(b_{ij}\) is a positive integer.

## 1. Local theorem at a good determinant prime

Fix an odd split prime \(p\mid\delta\) with \(p\nmid K\). Write
\[
 \ell=v_p(\delta),\quad c_i=v_p(r_i),\quad
 \beta_{ij}=v_p(b_{ij}).
\]
Let \(t_i\) be the signed Gaussian orientation valuation of the primitive
output \(B_i\). Thus \(|t_i|=v_p(NB_i)\), and the local exponent of the
least squared circle radius is
\(W=\max_i t_i-\min_i t_i\).

**Joint pair inequality.** One has
\[
 \boxed{
 \beta_{ij}\ge\min(c_i,c_j)
 +\mathbf1_{t_it_j>0}\min(|t_i|,|t_j|).}
 \tag{1}
\]
In particular,
\[
 \boxed{\sum_i c_i\le\ell+\sum_{i<j}\beta_{ij},}
 \tag{2}
\]
and the stronger joint loss estimate is
\[
 \boxed{
 2\sum_i c_i+\sum_i|t_i|-W
 \le3\ell+\sum_{i<j}\beta_{ij}.}
 \tag{3}
\]

For a proof of (1), work in isotropic coordinates over \(\mathbb Z_p\).
The matrix becomes
\[
 L=\begin{pmatrix}\alpha&\beta\\\gamma&\eta\end{pmatrix},
 \qquad K=16\alpha\beta\gamma\eta.
\]
All four entries are units. Also \(0\le c_i\le\ell\), by ordinary
primitivity and the adjugate identity.

If \(c_i>0\), the source lies modulo \(p\) in the kernel of \(L\).
That kernel has both isotropic coordinates nonzero, since both entries of
the first row of \(L\) are units. Thus the source norm is a \(p\)-unit.
In particular its Gaussian gcd with any other source row has norm prime to
\(p\).

If \(t_i\ne0\), the normalized output has exactly one isotropic coordinate
divisible by \(p\). Applying \(\operatorname{adj}(L)\), whose four entries
are units, gives two unit coordinates. Since
\[
 H_i=(r_i/\delta)\operatorname{adj}(L)B_i
\]
is primitive, necessarily \(c_i=\ell\). This proves the useful disjointness
statement
\[
 t_i\ne0\implies c_i=\ell,
 \qquad c_i<\ell\implies t_i=0.
 \tag{4}
\]

When both \(c_i,c_j>0\), the old Gaussian gcd contributes no \(p\)-factor.
The first row of \(L\) is a primitive linear functional, and takes values
divisible by \(p^{c_i},p^{c_j}\) on the source columns. Completing that
functional to an invertible coordinate system gives
\(\beta_{ij}=v_p\det(H_i,H_j)\ge\min(c_i,c_j)\).
If one content exponent is zero this inequality is automatic.

When the extra term in (1) is nonzero, (4) gives \(c_i=c_j=\ell\).
The determinant identity then gives
\[
 \beta_{ij}=\ell+v_p\det(B_i,B_j).
\]
Outputs of the same isotropic orientation have determinant valuation at
least \(\min(|t_i|,|t_j|)\). This proves (1) in every case.

For (2), pair every index other than a maximal-content index with that
index, and use (1). To prove (3), first note the two elementary inequalities
\[
 2\sum_i c_i-\sum_{i<j}\min(c_i,c_j)\le3\ell,
 \tag{5}
\]
\[
 \sum_i|t_i|-W\le
 \sum_{\substack{i<j\\t_it_j>0}}\min(|t_i|,|t_j|).
 \tag{6}
\]
For (5), order the \(c_i\) increasingly. The coefficients in the left side
are \(i-m+2\); only the last two are positive, and they are one and two.
Each \(c_i\le\ell\). Inequality (6) is the signed-width inequality for
\(m\ge3\), including zeros and the case when every sign is the same.
Adding (5)--(6) and applying (1) proves (3).

At an odd inert prime dividing \(\delta\) but not \(K\), all primitive
source and output norms are units, so the phase term vanishes. The matrix
has an entry that is a unit, since otherwise \(p\mid K\). The same primitive
linear-functional argument proves \(\beta_{ij}\ge\min(c_i,c_j)\).
Thus (2)--(3) remain valid, with all \(t_i=0\).

## 2. The coefficient three is sharp

Take \(p=17\), \(q=p^\ell\), and
\[
 M=\begin{pmatrix}1&q\\1&0\end{pmatrix},\qquad
 H_0=(1,0),\quad H_1=(q,2),\quad H_2=(2q,1).
\]
The source columns are primitive and have odd norms. The matrix is primitive,
\(\delta=-q\), and \(K=q^4+4\), which is prime to \(p\). The row content
exponents are \((0,\ell,\ell)\), while the normalized outputs are
\[
 (1,1),\quad(3,1),\quad(3,2).
\]
Their norms are \(2,10,13\), all prime to \(17\), so every \(t_i=0\).
The only positive \(\beta_{ij}\) is \(\beta_{12}=\ell\), since the
source norms are \(p\)-units and their pair determinant is \(-3q\).
Thus the left side of (3) equals \(4\ell\), and its pair-content term
equals \(\ell\). The remaining \(3\ell\) cannot be lowered uniformly.
This example even retains the unit anchor; it makes no endpoint-clustering
claim.

## 3. A global determinant gain

Define
\[
 D=\prod_{\substack{p\mid\delta\\p\nmid K}}
 p^{v_p(\delta)},\qquad E=\binom m2.
\]
Combining (3) at primes dividing \(D\) with determinant-free pair transfer
at the other primes gives
\[
 \boxed{
 N'\ge\frac{D^{2m-3}\prod_i f_i}
 {2T^mK^E\prod_{i<j}b_{ij}}.}
 \tag{7}
\]
Here \(N'\) is the least squared realization radius of the entire output
phase configuration; no source anchor is needed for (7).

For details, let \(A_i=N(MH_i)=r_i^2N(B_i)\). At an odd prime dividing
\(D\), the exponent of \(\prod A_i\) minus the radius exponent is bounded
by \(3v_p(\delta)+\sum\beta_{ij}\), by (3). At any other odd prime,
ordinary contents cost at most \(2m v_p(\delta)\), while the signed-width
inequality and pair transfer cost at most
\(E v_p(K)+\sum\beta_{ij}\). At two, if \(2\mid D\), use the joint
ramified inequality proved in Section 5 below; it has the same conductor
charge \(3v_2(\delta)\), together with a single extra factor two. If
\(2\nmid D\), retain the ordinary-content bound and the usual single extra
factor two: if exactly \(r\) normalized outputs
are odd-odd, their product norms contribute \(r\), and their Gaussian pair
gcd norms contribute \(\binom r2\), with
\(r-\binom r2\le1\). Pair transfer bounds these gcds by \(Kb_{ij}\).
Consequently
\[
 N'\ge\frac{\prod_i A_i}
 {2D^3(|\delta|/D)^{2m}K^E\prod_{i<j}b_{ij}}.
\]
The least eigenvalue of \(M^TM\) is at least \(\delta^2/T\), so
\(\prod A_i\ge(\delta^2/T)^m\prod f_i\). Cancelling proves (7).

Under the anchored endpoint hypotheses used in the main radius theorem,
this specializes to
\[
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{D^{2m-3}N^{(3m-4)/8}}{T^mK^E}.
 \tag{8}
\]
The factor \(D^{2m-3}\) is a genuine gain from the determinant at primes
away from \(K\), including two when \(K\) is odd. Its benefit must be
weighed against the larger exponent
\(E\) on \(K\). The previously established bound with denominator
\(T^mK^{m-2}\) remains available independently; one may take the larger
of the two lower bounds. No claim is made here that (8) dominates it for
all matrices, or that it controls the remaining moving-source problem.

## 4. A compressed-discriminant variant with an anchor

When \(H_0=1\), the good-prime gain can also be combined with the anchored
signed-width inequality. The result is
\[
 \boxed{
 N'\ge\frac{D^{2m-3}\prod_i f_i}
 {2T^{m+2}K^{m-2}\prod_{i<j}b_{ij}}.}
 \tag{9}
\]
At odd primes dividing \(D\), retain the joint estimate (3). At other odd
primes, use the anchored inequality
\[
 \sum_i|t_i|-W\le(m-2)v_p(K)
 +\sum_{i<j}\beta_{ij}+2|t_0|,
\]
and charge ordinary contents by \(2m v_p(\delta)\). Multiplying these local
estimates produces an extra anchor factor bounded by
\((N B_0)^2\); it is harmless to include all its primes, including primes
already treated by (3). Because \(H_0=1\),
\[
 N B_0=N(MH_0)/r_0^2\le T.
\]
Thus the anchor costs at most \(T^2\), and the same eigenvalue comparison
used for (7) proves (9).

Here is a separate parity check. If \(2\mid D\), apply the joint ramified
inequality of Section 5, which already includes the single extra factor
two. Otherwise the integer
\(K=(A-D)^2+4B^2\) has two-adic valuation either zero or at least two.
If \(2\mid K\), then the factor \(2K^{m-2}\) supplies exponent at least
\(1+2(m-2)\ge m\), covering all primitive output norm parities.
If \(K\) is odd, pair transfer forces \(b_{ij}\) even for every pair
of odd-odd primitive outputs. Hence the remaining exponent is at most
\(r-\binom r2\le1\). Ordinary row contents at two are already included
in the \((|\delta|/D)^{2m}\) charge. This proves that the single displayed
factor two still suffices.

For a primitive anchored endpoint cluster, equation (9) yields
\[
 N'\ge\frac12(2/C)^{2(m-1)+E}
 \frac{D^{2m-3}N^{(3m-4)/8}}{T^{m+2}K^{m-2}}.
 \tag{10}
\]
Relative to the existing determinant-free compressed bound, this improves
the lower bound by \(D^{2m-3}/T^2\). One can therefore multiply that
existing bound by
\[
 \max\{1,D^{2m-3}/T^2\}.
\]
There is no assertion that this gain is large for every nonconformal map.

## 5. Ramified determinant primes are also usable

Suppose \(\ell=v_2(\delta)>0\) and \(K\) is odd. Because
\(K\equiv T^2\pmod4\), the matrix has an odd number of odd entries.
It cannot have three odd entries: any such matrix is invertible modulo two.
Thus \(M\) has exactly one odd entry, and its image modulo two is a
coordinate axis.

Write \(c_i=v_2(r_i)\), and let \(o_i\in\{0,1\}\) indicate that the
primitive output \(B_i\) has two odd coordinates. Its norm has valuation
\(o_i\). If \(o_i=1\), then \(\operatorname{adj}(M)B_i\) has at least
one odd coordinate, since the adjugate likewise has exactly one odd entry.
The identity
\[
 H_i=(r_i/\delta)\operatorname{adj}(M)B_i
\]
and primitivity force \(c_i=\ell\).

The source norms are odd, so every old Gaussian pair gcd norm is odd.
Consequently \(\beta_{ij}=v_2\det(H_i,H_j)\). A matrix row containing
the odd entry is a primitive linear functional over \(\mathbb Z_2\),
so the same coordinate argument as before gives
\(\beta_{ij}\ge\min(c_i,c_j)\). If both outputs are odd-odd, then both
contents equal \(\ell\) and their output determinant is even. Hence
\[
 \beta_{ij}=\ell+v_2\det(B_i,B_j)\ge\ell+1.
\]
Together these prove
\[
 \boxed{\beta_{ij}\ge\min(c_i,c_j)+o_io_j.}
 \tag{11}
\]
Set \(r=\sum_i o_i\). The same content inequality (5), followed by
\(r-\binom r2\le1\), gives
\[
 \boxed{2\sum_i c_i+r\le3\ell+\sum_{i<j}\beta_{ij}+1.}
 \tag{12}
\]
The left side is exactly the two-adic exponent of \(\prod_iN(MH_i)\).
Thus (12) supplies the full required local estimate, even without subtracting
any two-adic contribution of the least radius. The extra one is precisely
the single factor two already present in (7)--(10). This justifies including
the complete two-part of \(\delta\) in \(D\) whenever \(K\) is odd.

## 6. Nearly conformal diagonal maps

For the particular family
\[
 M_n=\operatorname{diag}(n+1,n),\qquad n\ge1,
\]
one has
\[
 \delta=n(n+1),\qquad T=2n^2+2n+1,\qquad K=(2n+1)^2.
\]
Every prime divisor of \(n(n+1)\) is coprime to \(2n+1\), so the full
good determinant part is \(D=\delta\), including its two-part. The elementary
bounds
\[
 D\ge n^2,\qquad T\le5n^2,\qquad K\le9n^2
\]
give
\[
 \frac{T^{m+2}K^{m-2}}{D^{2m-3}}
 \le5^{m+2}9^{m-2}n^6.
\]
Thus under the anchored endpoint hypotheses,
\[
 \boxed{
 N'\ge
 \frac{(2/C)^{2(m-1)+E}}
 {2\,5^{m+2}9^{m-2}}
 \frac{N^{(3m-4)/8}}{n^6}.}
 \tag{13}
\]
For \(m\ge5\), radius-nonincreasing images in this family require
\[
 \boxed{n\gg_{m,C}N^{(m-4)/16}.}
 \tag{14}
\]
The power six in the map-size loss is independent of \(m\). This is a
stronger necessary condition for this family, whose geometric distortion
tends to one, but it is still not a bound on arbitrary rational maps or on
the size of an endpoint cluster.

## 7. Persistent exact checks

The root agent extended
[the conductor checker](check_mobius_conductor_transfer.py) to verify both
global determinant-gain inequalities in 7,808 nonconformal matrix/source
cases. It uses the complete part of the determinant coprime to K, including
two, and checks the endpoint consequence on an actual Pell cluster.
These finite tests supplement the local proofs above.
