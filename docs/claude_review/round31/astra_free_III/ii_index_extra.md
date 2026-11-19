**PROVED — pruning the high-index tuples.** Suppose a shortened one-unit cluster is known to have full-rank valuation lattices for all five-tuples and every prime dividing their indices is at most \(Q\). There are absolute \(\eta_0,c>0\) and a constant \(K_Q\) with the following properties.

1. Deleting at most \(K_Q\) marks eliminates every \(\eta_0\)-balanced five-tuple of index larger than one. Ordinary relative balance means that the fifteen nonconstant unoriented cut weights lie in \([(1-\eta_0)w,(1+\eta_0)w]\) for some tuple-specific \(w\).
2. There is a subcluster of size at least \(cm-K_Q\) in which **every** five-tuple has index one.

The arguments apply to genuine points and gap-admissible angle models whenever the stated bound on the possible primes holds.

*Proof.* For the second claim, write \(f_i(\ell)=\pm1\) for the side of row \(i\) at threshold layer \(\ell=(p,h)\), weighted by \(\log p/W\). Use \(G_{ij}=\mathbb E f_i f_j=1-2d_{ij}/W\). Shortening ensures \(G_{ij}<0\) for \(i\ne j\), and \(G\) is positive semidefinite. The symmetric nonnegative matrix \(A=I-G\) has largest eigenvalue at most one; its operator norm equals its largest eigenvalue. Indeed \(|x^t A x|\le |x|^t A|x|\) entrywise in \(x\). Thus \(\sum_{j\ne i}G_{ij}^2\le1\) for each \(i\). The graph of pairs with \(G_{ij}<-1/25\) has maximum degree at most \(625\). Greedy deletion of successive closed neighborhoods gives an independent set of size at least \(m/626\). Use this subset in the second claim.

The square-relation bound excludes \(q=2\) on the shortened arc. Fix an odd remaining \(q\) and a nonzero coefficient vector \(c\bmod q\) with \(\sum c_i=0\bmod q\). A tuple \(T\) having this affine relation satisfies \(b_T=\sum c_i a_{T_i}\equiv0\pmod q\) coordinatewise for any fixed choice of integral representatives. First multiply the class \(c\) by a scalar chosen to make its centered representatives satisfy
\[
 S_2=\sum c_i^2\le 5q(q+1)/12 .
\]
Indeed the average over the nonzero scalars on any nonzero coordinate equals \(q(q+1)/12\). Put \(L=\sum|c_i|\). In a single unit class the angles of \(\exp(2i\langle b_T/q,\phi\rangle)^q\) all lie in the same arc of width at most \(L\Delta\). Its center may depend on \(\sum c_i\) and the original common phase. Within each of its \(q\) inverse branches, two distinct sums therefore obey
\[
 \sum_p |b_{T,p}-b_{U,p}|\log p\ge qW/2-O(q)
\]
by the ordinary monomial gap for \((b_T-b_U)/q\), using \(L/q\le5/2\). For a branch choose one representative tuple per distinct sum and sample two representatives independently and uniformly, with replacement. For \(k\) representatives the mean separation is at least \((1-1/k)(qW/2-O(q))\).

Express each \(a_{T_i,p}\) as the sum of its threshold indicators \(B_{T_i,\ell}\). By the triangle inequality the mean separation divided by \(W\) is at most the weighted layer mean of \(\mathbb E|Z_\ell-Z'_\ell|\), where \(Z_\ell(T)=\sum c_i B_{T_i,\ell}\) and the same layer is used on the two independent tuples. Flipping all bits at a layer sends \(Z\) to \(\sum c_i-Z\) and preserves this expectation. Mixing over a common fair orientation flip and then over layers can only increase the Gini mean difference. Indeed it is concave in the underlying law, by \(\mathbb E|Z-Z'|=2\int F(t)(1-F(t))\,dt\). For any real law of finite variance it also satisfies
\[
 \mathbb E|Z-Z'|\le(2/\sqrt3)\sqrt{\operatorname{Var} Z}.
\]
For a proof, write \(Z=H(U)\) with \(H\) its quantile function and \(U\) uniform on \([0,1]\). Integration on the two triangles in the unit square gives \(\mathbb E|Z-Z'|=2\int_0^1(2u-1)(H(u)-\mathbb EZ)\,du\); now use Cauchy–Schwarz. This applies also to discrete laws. Under the mixed law, each representative tuple supplies a distribution of mean \(\tfrac12\sum c_i\) because of the common fair flip.

For the positive-density subcluster, every distinct pair in a tuple satisfies \(|G_{ij}|\le1/25\), so its contribution to the mixed variance is
\[
 \tfrac14c^tG_Tc\le\tfrac14(S_2+L^2/25)\le\tfrac14(1+5/25)S_2 .
\]
As all means agree, the same bound holds for the pooled variance. Its squared Gini is at most \((1+5/25)S_2/3\le q(q+1)/6\le2q^2/9<q^2/4\).

For the relative-balance claim, let \(\alpha=W_T/W\) for a given balanced tuple. Its conditional distribution on the thirty oriented nonconstant patterns after the common flip is within \(O(\eta_0)\) in total variation of the uniform law \(U_*\) there. Each pair is separated by eight of the fifteen nonconstant cuts and has separation \(>W/2\). Thus \(\alpha\ge 15/16-O(\eta_0)\), while \(\alpha\le1\). With the two constant patterns equiweighted, the full five-bit law is within \(O(\eta_0)\) of a mixture of \(U_*\) and the fair law \(U\) on all 32 patterns. These mixtures form a convex set, so the same holds after pooling. Under \(U\) and \(U_*\), the variances are \(S_2/4\) and \((16S_2-(\sum c_i)^2)/60\) respectively, by subtracting the constant-pattern contribution. Thus any such mixture has variance at most \(4S_2/15\). The perturbation changes variance by \(O(\eta_0 L^2)\) since all distributions are centered about the same mean. The squared Gini bound is now \(16S_2/45+O(\eta_0 L^2)\le 4q(q+1)/27+O(\eta_0 q^2)\). Because \(q\ge3\), a sufficiently small absolute \(\eta_0\) makes this less than \((1/4-\epsilon)q^2\) for an absolute \(\epsilon>0\).

In either case the gap compared with the earlier lower bound gives \(k=O(1)\) for large \(W\). Each exact sum \(b_T\) determines the coefficient assigned to each active row: subtract two realizations and use full affine independence of *all* exponent rows to force every coefficient in the resulting exact affine relation to vanish. Full affine independence follows from \(d_{ij}>W/2\): given an affine dependence, sample independently from its normalized positive and negative coefficient distributions. Their coordinates have the same mean \(\mu\in[0,1]\) at each normalized prime, and the expected coordinate distance is at most \(2\mu(1-\mu)\le1/2\), using \(|x-y|\le x+y-2xy\). This contradicts strict pair separation. All tuples considered have distinct marks, so knowing their weighted active-row measure determines their active mark set. Deleting these marks for the boundedly many sums in all branches and all finitely many \(q,c\) establishes both claims. Bounded \(W\) gives boundedly many possible exponent rows and can be absorbed in \(K_Q\).

**COMPUTED — independent check on the smallest odd primes.** For \(q=3,5,7\), the exact worst over projective affine-relation vectors, best over their scalings, worst over \(r\in[0,1]\) of the Gini under \((1-r)U+rU_*\) is respectively \(315/256,239/128,633/256\). The script enumerates every zero-sum vector up to permutation and scaling (4, 8 and 13 orbits), and checks the endpoints and interior stationary point of each rational quadratic. [Finite Gini certificate](ii_index_gini.py)
