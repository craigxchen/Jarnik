# Chan almost-square and axis-strip results: scale audit

This note audits whether Tsz Ho Chan, *Factors of almost squares and lattice
points on circles* (arXiv:1406.2230), gives a usable bound for the lifted
near-axis points in the current argument.  The source is the primary paper:
<https://arxiv.org/abs/1406.2230> (the PDF is
<https://arxiv.org/pdf/1406.2230>).

## What Chan actually proves

For sufficiently large `n=N^2`, Theorem 3 bounds the number of representations
`n=A^2+B^2`, `|B|<n^(1/4)(log n)^(1/7)`, by ten. Theorem 4 bounds it by
thirty-six with logarithmic exponent `1/14` for any sufficiently large n
having one representation with `|B_1|<=exp((log n)^(2/7))`.

Theorem 3 writes `Y^2=h(2N-h)` and `h=s t^2`, with s squarefree, then
uses simultaneous Pell estimates. Theorems 1 and 2 instead count ordinary
divisors near the square root of one square or suitably almost-square
integer. The following scale comparisons are our inferences from those
statements and the already-proved Gaussian identities.

## Comparison with the Gaussian lift

Here the lifted point associated with a Gaussian divisor
`A=u+iv`, `P=Norm(A)`, is

`Z=(N/P) A^2 = (N-h)+iY`,

so `|Z|=N` and `Z` lies on the integer circle
`X^2+Y^2=N^2`.  The squared radius in Chan's notation is therefore
`n=N^2`, while the current angular argument only gives

`|Y| <= C N^(3/4)`

for a fixed constant `C`.  Chan's Theorem 3 reaches only

`N^(1/2) (2 log N)^(1/7)`,

and Theorem 4 reaches
`N^(1/2) (2 log N)^(1/14)`.  The target strip is larger by a factor of
order `N^(1/4)` (up to logarithms), so neither theorem applies to the full
lifted set.  The square `n=N^2` has the trivial special representation
`(N,0)`, so Theorem 4's hypothesis is satisfied for this fixed circle, but
its conclusion still concerns only the much narrower strip.

The Pell parameter also does not transfer in the needed way.  In the lift,
putting `g=N/P` gives the exact identity

`h = 2 g v^2`.

Thus the squarefree coefficient in Chan's decomposition is
`sf(2g)`, not merely the small squarefree part of `v`.  For near-critical P, one only knows
`N^(1/2-1/M)<=g<=(C^2/4)sqrt(N)`. The displayed size information alone
gives no small bound on `sf(2g)` beyond `2g`. In particular the small v
estimate does not make the squarefree part of g small; no stronger claim
about the realized g-values has been proved.

## Candidate norms and the stretched-log hypothesis

Suppose a near-critical candidate has ordinary norm `d=Norm(A)` with
`d=N^(1/2+o(1))`, and the current pair-separation estimate gives, up to a
fixed factor,

`|v| <= N^(1/(2M))`.

Chan's exceptional-point condition has the stretched-exponential form
`|v| <= exp((log d)^(2/7))`.  Taking logarithms, and using
`log d=(1/2+o(1)) log N`, shows that this condition can only be inferred
from the displayed power bound once

`M >= (2^(-5/7)+o(1)) (log N)^(5/7)`,

with fixed multiplicative constants changing only the implicit threshold.
This verifies the `5/7` estimate.  It is a possible regime for invoking
Theorem 4 on one fixed norm, but it does not solve the counting problem:
the all-pair norm argument already makes the `binom(M,2)` ordinary norms
distinct.  Chan's at-most-36 conclusion counts representations of each
fixed `d`, not the number of distinct values `d` arising from the pairs.
For `M >= 4`, the same candidates do lie inside Chan's representation
strip for their individual radius (`d^(1/4)>=N^(1/8)`), but
this still supplies no aggregate bound over distinct norms.

The divisor theorem (Theorem 2) has the same mismatch: it concerns many
divisors of one fixed almost-square integer in one interval around that
integer's square root.  The present ordinary norms are indeed divisors of the common integer N.
However, the reduction gives neither the special almost-square
factorization required of N nor the required additive interval around
sqrt(N). Its interval extends up to `N^(1/2+1/M)`, whose excess over
sqrt(N) need not be of order `N^(1/4)` times a logarithm.

## Conclusion

Chan's paper gives valid special-radius, short-axis-strip results and
confirms the stretched-log threshold `M` of order `(log N)^(5/7)` for the
exceptional-point hypothesis.  It does not bound the current arbitrary
near-axis arc of length `O(N^(3/4))`, and even in that large-`M` regime its
fixed-norm representation count cannot be summed over the distinct pair
norms.  The slow-`M` case and the required uniform growth bound therefore
remain open by this route.
