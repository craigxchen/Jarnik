# Same-side interval overlap forces target pair residues

This note proves the proposed local residue bound for the coefficient
intervals in
[mobius_actual_anchor_normalization.md](mobius_actual_anchor_normalization.md).
It also gives its exact product consequence. The bound can be strong for a
specified balanced binary profile and a point interval; general coefficient
intervals can make it vacuous. No uniform endpoint conclusion is obtained.

## Setup and invariant pair residues

Write the real-linear rational map as
\(Mz=\alpha z+\beta\bar z\), with \(\alpha\beta\ne0\), and put
\(w=\beta/\alpha\), assuming \(Nw\ne1\). For primitive Gaussian source
columns \(H_i\), let \(t_i=v_\pi(H_i)-v_{\bar\pi}(H_i)\) at an odd split
prime \(p=\pi\bar\pi\). Write
\[
 d=v_\pi(w),\quad e=v_{\bar\pi}(w),\quad
 L=\min(d,-e),\quad R=\max(d,-e),\quad I=[L,R].
\]
The anchor cost is \(k_{i,p}=\operatorname{dist}(t_i,I)\).
Let \(B_i\) be the coordinate-primitive transformed columns and define
\[
 b'_{ij}=\frac{|\det(B_i,B_j)|}{N\gcd_{\mathbb Z[i]}(B_i,B_j)},
 \qquad \beta'_{ij}=v_p(b'_{ij}).
\]
The source and target columns are assumed pairwise nonparallel, so these
residues are nonzero integers. This residue is the absolute imaginary
coordinate of the primitive Gaussian numerator of the pair phase ratio.
It is invariant under a common nonzero Gaussian multiplier on all columns,
followed by coordinate normalization: the pair phase ratios do not change,
and their primitive Gaussian representatives are unique up to sign.
This proves invariance under source or target reanchoring and permits
factoring out the output-conformal coefficient \(\alpha\).

## Exact local formulas outside the interval

Use the split embedding into \(\mathbb Q_p\), writing
\(r_i=H_i/\bar H_i\) so that \(v_p(r_i)=t_i\). After removing \(\alpha\),
the map in isotropic coordinates is
\[
 (P,Q)\longmapsto(P+wQ,\ \bar w P+Q).
\]
Its determinant is \(1-w\bar w\). Set
\[
 \ell=v_p(1-w\bar w),\qquad
 \kappa=\ell-\min(0,d+e)\ge0.
\]
If \(d+e\ne0\), unequal valuations in the determinant difference give
\(\kappa=0\). A positive core term is possible only when \(d+e=0\).

If both \(t_i,t_j>R\), strict domination gives
\[
 v_p(P_i+wQ_i)=v_p(Q_i)+d,\qquad
 v_p(\bar wP_i+Q_i)=v_p(Q_i).
\]
Thus their normalized output signed valuations both equal \(d\). Removing
the output coordinate contents and then the pair Gaussian gcd gives exactly
\[
 \beta'_{ij}
 =\ell+v_p(r_i-r_j)-d
 =v_p(r_i-r_j)-R+\kappa.
 \tag{1}
\]
Here \(\min(0,d+e)-d=-R\). Since
\(v_p(r_i-r_j)\ge\min(t_i,t_j)\),
\[
 \beta'_{ij}\ge\min(t_i-R,t_j-R)+\kappa. \tag{2}
\]
If \(t_i\ne t_j\), this inequality is an equality. Equal source valuations
can give additional cancellation and hence a larger residue.

If both \(t_i,t_j<L\), the symmetric calculation uses \(1/r_i\).
Both normalized output signed valuations equal \(-e\), and
\[
 \beta'_{ij}
 =\ell+v_p(r_i^{-1}-r_j^{-1})-e
 =v_p(r_i^{-1}-r_j^{-1})+L+\kappa
 \ge\min(L-t_i,L-t_j)+\kappa.
 \tag{3}
\]
Again distinct source valuations give equality in the last inequality.
Consequently, whenever the two source valuations are strictly on the same
side of the coefficient interval,
\[
 \boxed{\beta'_{ij}\ge\min(k_{i,p},k_{j,p}).} \tag{4}
\]
There is no omitted determinant-core loss: its effect is the nonnegative
\(\kappa\) in (2)--(3). The argument is local at an odd split prime and
requires no fixed Gaussian-unit class for the target columns.

At an interval endpoint the strict domination used above can fail, so the
exact formulas (1) and (3) must not be asserted there. The desired lower
bound is nevertheless zero whenever one row lies in the interval, and
ordinary integrality of \(b'_{ij}\) suffices for that zero bound. There is
no claimed positive same-side overlap for rows on opposite sides.

## Clipping identity and global cost bound

Put \(c_i=\operatorname{clip}_I(t_i)\) and \(k_i=|t_i-c_i|\). For each
pair, let \(o_{ij}=\min(k_i,k_j)\) if both lie strictly on the same side,
and let \(o_{ij}=0\) otherwise. Directly in the three possible regions,
\[
 2o_{ij}=k_i+k_j+|c_i-c_j|-|t_i-t_j|.
\]
Summing gives the exact identity
\[
 \boxed{2\sum_{i<j}o_{ij}
 =(m-1)\sum_i k_i+\sum_{i<j}|c_i-c_j|
   -\sum_{i<j}|t_i-t_j|.} \tag{5}
\]
Apply this independently at each odd split prime. Define
\[
 k_i^{\rm glob}=\prod_p p^{k_{i,p}},\quad
 f_{ij}=\prod_p p^{|t_{i,p}-t_{j,p}|},\quad
 J_{ij}=\prod_p p^{|c_{i,p}-c_{j,p}|},\quad
 B'=\prod_{i<j}b'_{ij}.
\]
The factors \(k_i^{\rm glob}\) are exactly the anchor costs from the cited
normalization note. Equations (4)--(5) imply
\[
 \boxed{
 (\prod_i k_i^{\rm glob})^{m-1}\prod_{i<j}J_{ij}
 \le (B')^2\prod_{i<j}f_{ij}.} \tag{6}
\]
Indeed the right-minus-left prime exponent is nonnegative at every odd
split prime, and all remaining prime factors of \((B')^2\) can only help.
If the primitive source squared radius is \(N\), Plotkin's bound gives
\(\prod f_{ij}\le N^{m^2/4}\). In particular
\[
 \min_i k_i^{\rm glob}
 \le\left((B')^2N^{m^2/4}/\prod J_{ij}\right)^{1/[m(m-1)]}.
 \tag{7}
\]
Using a short-arc estimate for \(B'\) requires its own target normalization
and unit hypotheses; (6) itself does not assume such an estimate.

## Balanced point intervals and the failure boundary

For example, let \(m\) be even. At each source prime suppose the valuations
are two levels separated by \(n_p\), each occupied by \(m/2\) rows. Suppose
also that the coefficient interval is a point inside those levels. The two
same-side groups have distances \(u\) and \(n_p-u\). Hence
\[
 \sum_{i<j}o_{ij}
 =\binom{m/2}{2}n_p=\frac{m(m-2)}8n_p.
\]
For such a profile at every prime, (4) gives
\(B'\ge N^{m(m-2)/8}\). If an independently justified target estimate is
\(B'\le(C'/2)^{\binom m2}(N')^{m/8}\), then
\[
 N'\ge (2/C')^{4(m-1)}N^{m-2}.
\]
This is an obstruction for that specified profile. It is not a statement
for arbitrary nested source valuations or arbitrary target unit classes.

In general an interval may contain every source valuation. Then every
anchor cost and same-side overlap is zero at that prime; the residue bound
has no force there, even when the source pair distances are large. More
generally, large pair distances can be carried by the clipped distances
inside a long interval, or by opposite sides, without creating the same-side
quantity needed by (4). Therefore source minimum-distance information alone
does not yet supply the missing positive lower bound for this overlap.
Any application beyond a specified profile must control the coefficient
intervals and the source distribution relative to them.

An independent exact rational test at the split prime five checked 602
strict same-side instances with \(d,e,t_i,t_j\in[-3,3]\), including negative
coefficient valuations. A separate check verified 2,268 pairwise clipping
identities on integer intervals. All checks passed; the arguments above
establish the statements at arbitrary valuations.

Subsequent work combines this overlap with the exact source cut deficit
and the clipped conductor that divides the target. It bounds the aggregate
amount erased without prescribing a balanced binary profile; see
[the removed-conductor budget](mobius_erased_conductor_budget.md).
That comparison still permits unchanged radius and proves no uniform count.
