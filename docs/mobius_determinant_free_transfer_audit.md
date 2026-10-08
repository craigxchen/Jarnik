# Determinant-free pair-content transfer

Let \(M\) be an integral \(2\times2\) matrix of nonzero determinant
\(\delta\), and put
\[
 Q=M^TM=\begin{pmatrix}A&B\\B&D\end{pmatrix},\qquad
 K=(A-D)^2+4B^2>0.
\]
Let \(H_i,H_j\) be nonparallel integer-primitive Gaussian integers of odd
norm. Define
\[
 e=N(\gcd_{\mathbb Z[i]}(H_i,H_j)),\qquad
 b=|\det(H_i,H_j)|/e.
\]
For each row divide \(MH_r\) by its ordinary coordinate gcd \(c_r>0\),
obtaining an integer-primitive Gaussian integer \(B_r\). Then
\[
 \boxed{N(\gcd_{\mathbb Z[i]}(B_i,B_j))\mid Kb.}
\]
In particular no factor \(|\delta|\) is required. The proof does not require
that the entries of \(M\) have gcd one.

## Odd primes

Inert primes divide no primitive row norm, so consider an odd split prime
\(p\). Choose a square root \(\iota\) of \(-1\) in \(\mathbb Z_p\), and
use coordinates \((z,\bar z)=(x+\iota y,x-\iota y)\). This coordinate
change and its inverse are integral over \(\mathbb Z_p\), with unit
determinant. In these coordinates \(M\) becomes
\[
 L=\begin{pmatrix}\alpha&\beta\\\gamma&\eta\end{pmatrix},
 \qquad \det L=\delta,
\]
where
\[
 \alpha=\frac{a+d+\iota(c-b)}2,\quad
 \beta=\frac{a-d+\iota(c+b)}2,\quad
 \gamma=\frac{a-d-\iota(c+b)}2,\quad
 \eta=\frac{a+d-\iota(c-b)}2.
\]
Direct multiplication gives
\[
 K=16\alpha\beta\gamma\eta.
\]
All four entries are nonzero because \(K>0\), and all have nonnegative
valuation. Write
\(a_0=v_p(\alpha)\), \(b_0=v_p(\beta)\),
\(g_0=v_p(\gamma)\), \(d_0=v_p(\eta)\), and
\(k=v_p(K)=a_0+b_0+g_0+d_0\).

Let \(s=v_p(N\gcd(B_i,B_j))\). If \(s\le k\), there is nothing to
prove. Otherwise the two normalized output rows share one isotropic
orientation. Exchange the two isotropic coordinates if necessary so that
both output rows \((U_r,V_r)\) satisfy
\[
 v_p(U_r)\ge s>k,\qquad v_p(V_r)=0.
\]
The inverse numerator is
\[
 \operatorname{adj}(L)(U_r,V_r)^T
 =(\eta U_r-\beta V_r,\ -\gamma U_r+\alpha V_r)^T.
\]
Its two valuations are exactly \(b_0,a_0\): the terms involving \(U_r\)
have strictly greater valuation than either potentially competing term.
This strict comparison is the reason for first disposing of \(s\le k\).

Put \(\ell=v_p(\delta)\) and \(u=\min(a_0,b_0)\). Since
\(H_r=(c_r/\delta)\operatorname{adj}(L)(U_r,V_r)^T\) and the source row
is integer primitive, its minimum isotropic-coordinate valuation is zero.
Thus for both \(r=i,j\),
\[
 v_p(c_r)=\ell-u,
\]
and the two source coordinate valuations are exactly
\((b_0-u,a_0-u)\). In particular their common Gaussian gcd norm exponent is
\[
 v_p(e)=|a_0-b_0|.
\]
Writing \(D'=v_p(\det(B_i,B_j))\), determinant invariance gives
\[
 \begin{aligned}
 v_p(b)
 &=v_p(c_i)+v_p(c_j)+D'-\ell-v_p(e)\\
 &=\ell-a_0-b_0+D'.
 \end{aligned}
\]
Both first output coordinates are divisible by \(p^s\), so \(D'\ge s\).
As \(a_0+b_0\le k\) and \(\ell\ge0\), this yields
\[
 v_p(b)\ge s-k,
\]
which proves the required divisibility at every odd split prime, including
primes dividing \(\delta\).

## The prime two

An integer-primitive Gaussian integer has valuation at \(1+i\) at most
one. Therefore the target gcd norm has two-adic exponent either zero or
one. If \(K\) is even, its valuation already covers this exponent.

Suppose \(K\) is odd and the target gcd norm is even. Both normalized
output rows then have two odd coordinates. We claim \(b\) is even.
Since the source norms are odd, \(e\) is odd. If \(b\) were odd, the two
source vectors would form a basis modulo two. Each unnormalized image is
its normalized odd-odd vector multiplied by an integer, so each image
modulo two lies in the line spanned by \((1,1)\). Consequently the entire
image of \(M\) modulo two lies in that line. Each column of \(M\) has even
squared norm, making \(A-D\) even and hence \(K\) even, a contradiction.
Thus \(b\) is even, proving the required final prime valuation.

This establishes the determinant-free transfer claim. Its use in a global
radius estimate still requires the separate signed-valuation argument,
source normalization, and bounds on the ordinary coordinate contents.
