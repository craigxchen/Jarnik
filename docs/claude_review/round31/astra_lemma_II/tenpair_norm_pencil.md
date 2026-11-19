# A finiteness theorem for a fixed rational pencil of block norms

**PROVED.** Fix a rational subspace \(U\subset\mathbb Q^{10}\) of dimension at most two and a constant \(T\). In the genuine ten-pair-cut model of Part I, §5, only finitely many private pairwise coprime block-norm vectors \(n=(n_{ij})\in U\) have \(0<|t_{ij}|\le T\) for all pairs. This allows the norms to grow at unrelated rates a priori. The finiteness assertion is individual to \(U,T\).

## From the cofactor bounds to a conic

**PROVED.** Write \(P_i=\prod_{j\ne i}\beta_{ij}\), \(n_{ij}=|\beta_{ij}|^2\), \(t_{ij}=\det(P_i,P_j)/n_{ij}\), and \(s_{ij}=\operatorname{Re}(\bar P_iP_j)/n_{ij}\). Both quotients are integral. The norm gives
\[
 s_{ij}^{\,2}+t_{ij}^{\,2}=\prod_{a\notin\{i,j\}}n_{ia}n_{ja}.
 \tag{1}
\]
For distinct \(i,j,k\) and complementary indices \(l,m\), define
\[
 A_{ijk}=t_{jk}^{\,2}n_{jk}n_{il}n_{im}
 +t_{ik}^{\,2}n_{ik}n_{jl}n_{jm}
 -t_{ij}^{\,2}n_{ij}n_{kl}n_{km}.
\]
Put \(m_{ij}=t_{ij}n_{ij}\). Taking squared absolute values in \(m_{jk}P_i-m_{ik}P_j=-m_{ij}P_k\), using \(|P_i|^2=\prod_{a\ne i}n_{ia}\), and dividing by \(n_{ij}n_{ik}n_{jk}\) gives \(2t_{ik}t_{jk}s_{ij}=A_{ijk}\). In particular \(s_{ij}\), with a fixed choice of \(k\), is a real homogeneous cubic in the \(n\)'s once the \(t\)'s are fixed.

Fix the signed \(t\)'s along a putative infinite subsequence. Dimension one gives finiteness by (1) and the homogeneous cubic, since their difference is a homogeneous sextic taking a nonzero constant value. In dimension two choose a lattice basis of \(U\cap\mathbb Z^{10}\). Write \(n_e=L_e(p,q)\) for rational integer linear forms and integer \(p,q\). Substituting for fixed \(i,j,k\) gives one homogeneous rational binary sextic \(F(p,q)=c\ne0\), with \(c=t_{ij}^2\) and \(F=\prod_{a\notin\{i,j\}}L_{ia}L_{ja}-s_{ij}(L)^2\). There are infinitely many projective parameter slopes since \(F\) takes a fixed nonzero value.

After clearing denominators each primitive integer irreducible factor of \(F\) has a fixed nonzero integer value on an infinite subsequence, by divisibility. Any irreducible factor of degree at least three is excluded by Thue's theorem on binary forms. Two independent rational linear factors fix \(p,q\). A linear factor with fixed value confines the parameters to an affine line. On this line all equations (1), with their cubic substitutions, hold identically by the infinite subsequence. At least one \(L_{ab}\) is nonconstant since the map into \(U\) is injective on parameters. Evaluating (1) for \(i=a\) and \(j\notin\{a,b\}\) at its real root is impossible since \(t_{ij}\ne0\). Finally, two nonproportional irreducible quadratic factors with fixed nonzero values would give a nonzero homogeneous quadratic vanishing on infinitely many slopes by taking a linear combination. A definite quadratic with fixed nonzero value has finitely many integer representations. The only remaining possibility is an indefinite nonsplit rational quadratic \(Q(p,q)=d\ne0\).

**PROVED.** Diagonalize the real conic as \(u v=1\) by a real linear change of the two coordinates. The point \((p,q)\) on the conic is a real linear combination of \(u\) and \(u^{-1}\). Pass to one of its two real branches and take \(u>0\). Every \(n_e\) now has the form \(a_eu+b_e/u\) with \(a_eb_e\ne0\), since the isotropic lines of a nonsplit rational quadratic have irrational slopes and \(L_e\) is rational and nonzero. All the norm identities (1) and cubic reconstruction identities extend over the conic by equality at infinitely many parameters. Every \(n_e\) stays positive on \(u>0\). Indeed it is positive at actual tuples, and (1) for a pair whose product includes \(n_e\) precludes a real zero. Hence \(a_e,b_e>0\). Distinct forms \(L_e,L_f\) are nonproportional. If they were proportional, write \(n_e/n_f=A/B\) in lowest terms. Coprimality forces \(n_e=A,n_f=B\) at all sample points, impossible for a nonconstant Laurent polynomial taking a fixed value infinitely often. Consequently the ten ratios \(b_e/a_e\) are distinct.

Define \(V_{ij}=(s_{ij}+it_{ij})/(s_{ij}-it_{ij})\), a rational function of \(u\). On every sample point it equals \((P_j/\bar P_j)/(P_i/\bar P_i)\), so \(V_{ij}V_{jk}=V_{ik}\) identically. At zero and infinity, (1) shows that \(s_{ij}\) has a pole of order three. Thus \(V_{ij}=1+O(u^3)\) at zero and \(1+O(u^{-3})\) at infinity. The only remaining possible zeros and poles of \(V_{ij}\) are the two roots \(\pm i\sqrt{b_e/a_e}\) for any of the six edges occurring in (1). At each such root exactly one of \(s_{ij}\pm it_{ij}\) vanishes simply, since the six edge forms occurring in (1) have disjoint simple roots and \(2it_{ij}\ne0\). Real conjugation exchanges the zero and the pole. At a fixed edge root the valuations of \(V_{ij}\) form a cocycle on the five indices. The valuation vanishes for pairs \(ij\) on the same side of the cut given by that edge and is \(+1\) or \(-1\) for pairs on opposite sides. A cocycle with these properties is a signed difference of the indicators of membership in the edge. We can thus choose \(y_e=\pm\sqrt{b_e/a_e}\) independently for each edge so that
\[
 V_{ij}=\prod_e\left(\frac{u-iy_e}{u+iy_e}\right)^{
                  {\bf1}_{j\in e}-{\bf1}_{i\in e}} .
\]
The divisor determines the right side up to its constant, which is fixed by \(V_{ij}(\infty)=1\). Comparing linear terms at infinity and logarithmic derivatives at zero shows that both \(\sum_{e\ni i}y_e\) and \(\sum_{e\ni i}1/y_e\) are independent of \(i\). The following real lemma excludes such a configuration.

## A real five-star lemma

**PROVED.** Suppose \(y_{ij}=y_{ji}\) for the ten edges of \(K_5\) are nonzero reals distinct up to sign. Their star sums and their reciprocal star sums cannot both be constant.

*Proof.* Suppose the star sums are \(s\) and the reciprocal sums \(r\). Then the five monic quartics of the star values are
\[
 P_i(z)=\prod_{j\ne i}(z-y_{ij})
       =z^4-sz^3+c_i z^2-r d_i z+d_i.
\]
If \(r=s=0\), each \(P_i\) is even and its roots cannot be distinct up to sign. Otherwise consider the real rational parametric plane curve
\[
 C(z)=[1-rz:z^4-sz^3:z^2].
\]
The degree of the parameterization as a linear system is four. The homogenized forms have no base point. On the open third-coordinate chart, two coordinates are \(Y=z^2-sz\) and \(X=z^{-2}-rz^{-1}\). Since \(Y\) has degree two, failure of birationality requires its involution \(z\mapsto s-z\) to fix \(X\). The only pole of \(X\) is zero, forcing \(s=0\), then \(r=0\) by invariance. Thus \(C\) is an irreducible plane quartic with normalization \(\mathbb P^1\). At the distinct points \(C(0)\) and \(C(\infty)\) the normalization has only one preimage, as the third coordinate can vanish only at these two parameters. Its differential vanishes at each because both vanishing local coordinates have order at least two. These are two distinct real unibranch singularities.

The five lines \(w_1+c_iw_2+d_iw_0=0\), in coordinates \([w_0:w_1:w_2]\), are distinct. Equality would identify two star quartics with four common roots although their star-value sets intersect in exactly one element. The intersection of lines \(i,j\) is \(C(y_{ij})\). No three lines concur, since that would give a point with three distinct finite nonzero preimages while \(Y=z^2-sz\) has degree two. A nonsingular real conic is tangent to all five lines. In the dual plane the five points with no three collinear impose five independent conditions on conics (four points are covered by two lines avoiding the fifth). The resulting conic cannot be reducible since two lines cannot cover five such points, even over \(\mathbb C\). The dual conic is real with real tangencies, hence is split over \(\mathbb R\).

Identify this plane with \({\rm Sym}^2\mathbb P^1\) over \(\mathbb R\) by assigning a point its unordered pair of tangents to the conic, counted with multiplicity over \(\mathbb C\). In coordinates on the dual conic this is the familiar identification with binary quadratics up to scale. The five lines become divisors of pairs containing one of five distinct real parameters \(h_1,\ldots,h_5\). Choose all \(h_i\) finite. Pull back a real equation of \(C\) to a symmetric polynomial \(D(x,y)\) of bidegree at most \((4,4)\) on the product. Since no one of the lines is contained in \(C\), writing \(f(x)=\prod_i(x-h_i)\) gives \(D(h_i,y)=k_i f(y)/(y-h_i)\) with \(k_i\in\mathbb R^\times\). Let \(g\) be the interpolation polynomial of degree at most four with \(g(h_i)=k_i\). Interpolation in the first coordinate shows
\[
 D(x,y)=\frac{f(x)g(y)-f(y)g(x)}{x-y}.
\]
Indeed both sides have degree at most four in each variable and agree at the five distinct first coordinates \(h_i\). Therefore \(C\) parametrizes unordered pairs in a common fiber of the real degree-five map \(h=f/g\), including a double point on the diagonal at every ramified preimage. The fiber of \(h\) at zero consists of five distinct real points.

We now classify the possible singularities of this plane curve over \(\mathbb C\). At an unordered pair of distinct points in a fiber with local degrees \(m,n\), the map taking ordered to unordered pairs is locally an isomorphism, and the completed analytic local equation is \(X^m-Y^n=0\), after a change of local parameters. This follows by taking a local coordinate on the target and its analytic \(m\)th and \(n\)th roots at the respective sources. The only singular unibranch possibility with \(m+n\le5\) is \(\{m,n\}=\{2,3\}\), an ordinary cusp of delta invariant one. At a diagonal point of local degree \(m\), the equation is \((X^m-Y^m)/(X-Y)\) in the symmetric parameters \(v=XY,w=X+Y\). For \(m=2,3\) the equations \(w=0,w^2-v=0\) are smooth. For \(m=4\) one gets \(w(w^2-2v)=0\). For \(m=5\) one gets \(w^4-3w^2v+v^2=0\), with two branches. The simple diagonal case does not belong to \(C\). It follows that the two unibranch points each force a different ramified fiber of type \((3,2)\). Both critical points in each of these fibers are real because the unordered pair is real and the local degrees distinguish its two members.

The two fibers account for six of the eight ramification units of a degree-five rational map. The genus formula for an irreducible rational plane quartic gives total delta invariant three. After the two cusps there must be another singularity. With ramification budget two, the sole choice that gives any singularity in the local equations above is one fiber of type \((2,2,1)\), giving an ordinary node. The other choices are one triple critical point or two simple critical points over distinct values, neither of which gives a singularity.

Apply real Möbius changes to put the two cusp values at zero and infinity, the triple zero at \(x=0\), the double zero at \(x=\infty\), and the triple pole at \(x=1\). The double pole is \(b\in\mathbb R\setminus\{0,1\}\). The new function, up to real scaling, is
\[
 H(x)=\frac{x^3}{(x-1)^3(x-b)^2}.
\]
The two remaining critical points are distinct roots of \(q(x)=2x^2+x-3b\). Their \(H\)-values coincide. Neither is a pole. Thus the coefficient vectors of the linear remainders of \(x^3\) and \((x-1)^3(x-b)^2\) modulo \(q\) are proportional. Their constant/linear determinant is \(b^2(9b+16)(24b+1)/16\), by polynomial division (the standalone exact verification is tenpair_belyi_check.py). The choice \(b=-1/24\) makes the two roots equal. Hence \(b=-16/9\), giving discriminant \(1+24b=-125/3<0\). The two remaining critical points are nonreal.

The only real critical values contributed by *real* source points are now zero and infinity. The number of real inverse images is constant on each interval of \(\mathbb P^1(\mathbb R)\) complementary to these two, since the restriction to the real circle is proper and a local diffeomorphism there. On either side, sufficiently close to zero there is one real inverse near the triple zero and either none or two near the double zero, by their real local degrees. Thus every regular real fiber has at most three real points. This contradicts the original fiber with five distinct real points, a property invariant under both real Möbius changes. This proves the lemma and the fixed-pencil finiteness theorem. \(\square\)

**PROVED — scope.** The theorem excludes bounded cofactors on each fixed rational plane of norm vectors. The proof gives neither a power of the parent norm for growing cofactors nor a bound independent of the plane. The polynomial obstruction for a two-dimensional real span of norms is more direct: a homogeneous binary sextic constant on a nonconstant polynomial parameterization forces a nonzero linear combination of the two parameter polynomials to be constant, by factorization over \(\mathbb C\). With real independent input take real and imaginary parts to obtain a real such combination. On the resulting affine line at least one norm form has a real root in the line parameter; equation (1) and cubic reconstruction give the same contradiction. This polynomial argument is independent of the real five-star lemma.
