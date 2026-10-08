# Exact primitive chord residues and the balanced-layer bonus

## Status

This note does **not** prove D11 of the multipoint continuation.
Section 6 proves a power lower bound for the residue on an exact balanced
stratum, with exponent \(1/14\) in prime-layer weight normalization.
The earlier sections preserve
all finite-place residue information in an explicit positive integer and
isolates the precise additional growth statement. The identities hold with
arbitrary conductor exponents and inherited constant layers. The common
Gaussian-unit assumption is essential to the particularly simple factor 4.

## 1. Pair factorization with inherited layers

Write the eight points in one common Gaussian-unit class as

\[
z_i=\epsilon\prod_p\pi_p^{a_i(p)}
                  \overline{\pi_p}^{\,e_p-a_i(p)},\qquad
N=|z_i|^2=R^2,\qquad
W=\log N=\sum_pe_p\log p.
\]

The conductor is inherited from the ambient normalized cluster: neither
the minimum nor the maximum of \(a_i(p)\) over these eight selected rows
must be \(0,e_p\). Thus constant inherited layers are retained.

For a pair, put

\[
g_{ij}=\prod_p\pi_p^{\min(a_i,a_j)}
                  \overline{\pi_p}^{\,e_p-\max(a_i,a_j)},
\]

\[
h_{ij}=\prod_{a_i>a_j}\pi_p^{a_i-a_j}
       \prod_{a_i<a_j}\overline{\pi_p}^{\,a_j-a_i}.
\]

Then exactly

\[
z_i=\epsilon g_{ij}h_{ij},\qquad
z_j=\epsilon g_{ij}\overline{h_{ij}}.
\]

The two numbers \(h_{ij}\) and \(\overline{h_{ij}}\) are coprime.
Distinctness and the common unit imply

\[
t_{ij}=|\operatorname{Im}h_{ij}|\in\mathbf Z_{\geq1}.
\]

If

\[
d_{ij}=\sum_p|a_i(p)-a_j(p)|\log p,
\]

then

\[
|g_{ij}|^2=N e^{-d_{ij}},\qquad
\boxed{|z_i-z_j|^2=4t_{ij}^2N e^{-d_{ij}}.} \tag{1}
\]

This is an equality, including all extra prime divisibility of the chord.

## 2. One explicit integer contains the entire determinant residue

For each inherited prime layer \(a=(p,t)\), let \(r_a\) be its size among
the eight points. Put

\[
D=\sum_a(r_a-4)^2\log p,\qquad
W_4=\sum_{a:r_a=4}\log p,\qquad
T=\prod_{i<j}t_{ij}\in\mathbf Z_{\geq1}.
\]

The cut count gives

\[
\sum_{i<j}d_{ij}=16W-D.
\]

Multiplying (1) over the 28 pairs yields

\[
\prod_{i<j}|z_i-z_j|^2=4^{28}T^2N^{12}e^D. \tag{2}
\]

The eight-row Ramana determinant \(J\) has norm equal to this product
divided by \(N^{12}\). Thus

\[
\boxed{N_{\mathbf Q(i)/\mathbf Q}(J)=4^{28}T^2e^D.} \tag{3}
\]

In particular, the non-defect part is an explicit square, not merely an
unspecified collection of nonnegative prime-valuation remainders.

There is an exact Gaussian form as well. Set

\[
q_p=\sum_{i=0}^7a_i(p)-4e_p,\qquad
D_p=\sum_{t=1}^{e_p}(r_{p,t}-4)^2.
\]

The integers \(D_p\pm q_p\) are even and nonnegative: each square
\((r-4)^2\) has the same parity as \(r-4\), and dominates its absolute
value. Define

\[
H_D=\prod_p\pi_p^{(D_p+q_p)/2}
            \overline{\pi_p}^{(D_p-q_p)/2}.
\]

With \(V=\prod_{i<j}(z_j-z_i)\) and
\(J=N^6V/\prod_i z_i^3\), direct substitution in the preceding pair
factorizations gives

\[
J=\pm2^{28}H_DT. \tag{4}
\]

Indeed, the exponent of \(\pi_p\) outside \(T\) is

\[
6e_p+\sum_{i<j}\min(a_i,a_j)-3\sum_i a_i
=\frac{D_p+q_p}{2},
\]

and the conjugate exponent is \((D_p-q_p)/2\).
The factor \(\epsilon\) cancels because its remaining power is
\(\epsilon^{28-24}=\epsilon^4=1\); also \(i^{28}=1\).
The sign absorbs the signs of the primitive imaginary parts and the
choice of ordering of the pair differences.

## 3. Exact geometric slack and the proposed bonus

Suppose all chords are at most \(C\sqrt R\), and define the nonnegative
geometric slack

\[
S_{\rm arc}=\sum_{i<j}
 \log\frac{C\sqrt R}{|z_i-z_j|}.
\]

Comparing (2) with the exact product of these chord lengths gives

\[
\boxed{D+2\log T+2S_{\rm arc}
       =2W+56\log(C/2).} \tag{5}
\]

At \(C=1/2\), the constant is \(-56\log4\). In particular, one gets the
unconditional refinement

\[
D+2\log T\leq2W-56\log4. \tag{6}
\]

The proposed D11 inequality

\[
D+\delta W_4\leq2W+B
\]

is therefore equivalent, on these endpoint configurations, to

\[
\boxed{\delta W_4\leq
2\log T+2S_{\rm arc}+B+56\log4.} \tag{7}
\]

A sufficient purely arithmetic strengthening would be

\[
2\log T\geq\delta W_4-B_0 \tag{8}
\]

for fixed \(\delta>0,B_0\). Neither (7) nor (8) is proved here.
The geometric slack is important: (8) is stronger than the precise
requirement (7). Very unequal angular gaps can supply some or all of the
needed contribution through \(S_{\rm arc}\).

The exact-balanced phase exclusion does not directly furnish (8).
For completely balanced layers, \(D=0\) and \(W_4=W\), so D11 for
\(0<\delta\leq1\) would already be numerically automatic. The critical
regime is instead \(D\) close to \(2W\), with substantial \(W_4\), such
as the binomial profile distribution in the sampling criterion.

## 4. Limits of a bounded-residue argument

Bounded primitive imaginary parts do not by themselves imply bounded
radius, even for fixed-constant endpoint arcs with four points.

Fix four distinct integers \(c_j\), chosen so that there is an infinite
arithmetic progression of positive integers \(t\) for which

\[
A_j(t)=t+c_j+i
\]

are pairwise coprime and mutually conjugate-coprime. Such choices can be
made explicitly by first taking all \(c_j\) to be distinct multiples of
a sufficiently large fixed primorial, then avoiding the finitely many
prime divisors of their nonzero pairwise resultants by CRT. At primes
larger than twice the number of forms, the finitely many forbidden
residues cannot exhaust the residue field; at the smaller primes the
common primorial makes all forms congruent, and \(t\equiv0\) makes them
congruent to \(i\).

Put \(G=\prod_jA_j\), \(R=|G|\), and

\[
z_j=A_j\overline{G/A_j}.
\]

For each pair, the primitive quotient is
\(A_i\overline{A_j}\), whose imaginary part is \(c_j-c_i\).
Thus every \(t_{ij}\), and their product, is fixed while \(R\to\infty\).
Moreover, \(R\asymp t^4\), and the phase span is \(O(t^{-2})\), so these
are endpoint arcs with a fixed constant depending on the shifts.

The same construction with eight forms gives eight unbounded-radius
configurations with bounded \(T\), but radius \(R\asymp t^8\) and angular
span of order \(t^{-2}\); it does not meet the endpoint scale. Its
conductor cuts are singleton cuts, so \(W_4=0\).

These examples rule out dropping the balanced-layer hypothesis from the
residue-growth target. They do not falsify (7), (8), or D11.

## 5. Remaining arithmetic problem

After keeping the exact square residue and the exact angular-gap slack,
the missing statement is (7), or a sufficient stronger statement such
as (8). The existing integer lower bounds \(t_{ij}\geq1\), local
valuation balance, and the rational-ray obstruction for exactly
balanced profiles do not prove positive logarithmic growth in \(W_4\).
No fixed positive \(\delta\) has been established for the general
mixed-layer statement. Section 6 below proves such a power gain on a
specified exactly balanced stratum.


## 6. A proved power lower bound on an exact balanced stratum

Here there is a genuine positive residue-growth result, stronger than
integer nonvanishing. Assume every conductor layer is exactly balanced
among the eight rows, and that some chosen zero row has no complementary
valuation row among the other seven. Equivalently, all its complementary
products \(B_i\) in the balanced model are nonunit. The complete
seven-row, 35-profile model satisfies this condition. At each split prime,
the threshold cuts are nested; if every layer has size four, those cuts
are identical. Reorienting that prime away from the chosen zero row
therefore produces exactly the original-root complementary-product model
used by the rational-ray theorem.

Then

\[
\boxed{T\geq\frac1{100}R^{1/14},\qquad
2\log T\geq\frac1{14}W-2\log100.} \tag{9}
\]

No endpoint assumption is needed for this result.

First consider the graph whose edges are the pairs with \(d_{ij}\geq
W/2\). This graph is connected. Indeed, fix a nontrivial vertex subset
\(A\), of size \(a\), and consider a balanced cut with
\(r=|A\cap S|\). Its number of separating pairs across the partition
\(A,A^c\) equals

\[
r(4-a+r)+(a-r)(4-r)
=4a-2r(a-r)\geq\frac{a(8-a)}2.
\]

After weighting and summing the cuts, the mean crossing distance is at
least \(W/2\). Hence at least one crossing edge belongs to the graph,
for every nontrivial partition, proving connectedness.

Put \(u_i=z_i/R\). For a graph edge, (1) gives

\[
|u_i-u_j|=2t_{ij}e^{-d_{ij}/2}
\leq 2T R^{-1/2}. \tag{10}
\]

If \(T\geq\sqrt R/100\), the first inequality in (9) already follows
from \(R>1\). Otherwise put \(\tau=T/\sqrt R<1/100\). Each graph edge
has a principal angular difference at most
\(2\arcsin\tau\leq4\tau\). Lift phases along a spanning tree rooted at
the chosen zero row. Every path has at most seven edges, so all lifted
phases lie in \([-28\tau,28\tau]\), an interval of width

\[
\Delta\leq56T R^{-1/2}<0.56<\pi/2.
\]

The rational-ray theorem in phase_frontier_continuation.md applies to
this original root and its nonunit \(B_i\), and yields

\[
\frac4{3\sqrt2}R^{-3/7}
\leq\Delta\leq56T R^{-1/2}.
\]

Consequently

\[
T\geq\frac1{42\sqrt2}R^{1/14}
>\frac1{100}R^{1/14},
\]

finishing the proof.

Since \(D=0\) and \(W_4=W\) here, (9) proves the stronger arithmetic
target (8) with \(\delta=1/14\) on this stratum. This is more than the
original exact-balance phase exclusion: it forces a power of the
conductor into the explicit primitive chord residues without assuming
the points initially occupy a short arc.

The hypothesis excluding complementary valuation rows is substantive.
If every row has its complementary row, no root has all \(B_i\)
nonunit; the rational-ray product has a unit-complement degeneracy.
In the common-unit model this is also a geometric reflection symmetry
of the eight points. The present proof gives no power residue bound in
that case. Nor does the crossing-distance argument establish the same
threshold graph connectivity for general mixed cuts: it used every cut
having size exactly four. Thus (9) does not prove the mixed-layer D11
or the capped local test.
