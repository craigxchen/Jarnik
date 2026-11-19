### Exact degree-eight count of primitive constant-offset divisors

**PROVED.** For eight distinct roots \(\rho_i=x_i+iy_i\) in the upper half-plane, write
\[
 h_i=(y_i,x_i y_i,3x_i^2 y_i-y_i^3).
\]
The greatest possible number of ternary moment-kernel vectors of weight four **up to sign** is **12**, or 24 oriented vectors. Equality holds precisely for configurations supporting the eight-row cube sign code of Part I. This is uniqueness of the code, with room for more than one choice of roots. Each vector up to sign gives a different support and thus a different squarefree real norm divisor.

**PROVED.** The three moment equations are exactly \(\sum_i m_i h_i=0\). The balance lemma of Part I for total weight four forces positive and negative mass two each. The ternary vectors therefore correspond to four distinct indices with \(h_i+h_j=h_k+h_l\). These points are distinct since \(y>0\) and \(x=(xy)/y\). Two different pairings on one four-support cannot both give equal sums, by subtraction. If repeated factors are allowed, one must also count every midpoint \(2h_i=h_j+h_k\). A bound that the minimum weight is *at least* four does not guarantee that any vector of weight four exists.

**PROVED.** No four distinct moment points lie on an affine line. A line with height \(Y\) varying has its other coordinates \(V,W\) affine in \(Y\). Its intersections solve \(Y^4+YW-3V^2=0\), a monic quartic with no cubic term. Four positive roots would have positive sum, contradicting Vieta. On a line of fixed positive height the surface is the nondegenerate parabola \(3V^2=YW+Y^4\), with no three collinear points.

**PROVED — exhaustive reduction.** Choose a generic linear functional and label the moment points in increasing order. For \(i<j<k<l\) the only possible equality of pair sums on those four points is \(h_i+h_l=h_j+h_k\). There are 70 candidate rows \(e_i+e_l-e_j-e_k\). For repeated-factor minima add the 56 rows \(e_i+e_k-2e_j\), \(i<j<k\). Given a span \(E\) of these relations, the formal point for label \(i\) is the image of \(e_i\) in \(\mathbb Q^8/E\). Any actual moment configuration satisfying those relations is a linear image of the formal points. Discard a span that identifies two formal points or puts four on a line. These exclusions persist in every extension. The formal affine dimension is \(7-\dim E\) because the rows all sum to zero. Thus no span of rank six or above can occur.

**PROVED — enumeration certificate.** Start at the zero span. At each node reduce all candidate rows modulo a basis of the span. Count the zero images. Group nonzero images by proportionality and extend by the first index in each group if that index exceeds the last basis index used. Every desired span is visited: take its lexicographically greedy independent basis. If a candidate for the next selected basis vector is proportional modulo the current span to an earlier-indexed row, that row also belongs to the target span, in contradiction to the greedy choice. At rank four the zero images and a single proportionality group are exactly the rows in the corresponding rank-five extension, so one may skip extensions containing fewer than the target count. This argument does not require any particular chamber of the ordering to be feasible. Extra chambers only enlarge the census. The calculations modulo 769 are exact for these ranks. Every candidate row has Euclidean norm at most \(\sqrt 6\), and differences of two point labels used to test affine dependences have norm \(\sqrt 2\). Nonzero minors of the stacked rows have size at most seven (all rows have sum zero), and absolute determinant at most \((\sqrt 6)^7<769\) by Hadamard. Ranks of the quotient incidence tests and collinearity tests therefore agree with those in characteristic zero.

**COMPUTED — ternary case.** With a threshold of 12 the numbers of processed nodes at ranks 0 through 5 were \(1,70,1799,20640,94190,50\); only two rank-four nodes and the 50 retained rank-five nodes reached threshold. Affine classification of the rank-five formal point sets gives just three classes. Orient an affine triangle of maximum area, map it to \((0,0),(1,0),(0,1)\), sort the eight rational vectors, and minimize over such choices. This is a complete affine invariant. In the order of first appearance the class counts are 28, 18, 4. Representatives, given as zero-based lexicographic indices among the 70 quadruples, are
\[
 (0,5,12,14,16),\quad(0,5,12,14,18),\quad(2,7,14,23,27).
\]
The first two give zero remainders for positive products in the exact three-moment test of Part I, §7 (last factors correspond respectively to pairs \((0,3)\) and \((1,2)\)). For the third, the rows of its nullspace matrix are
\[
 (2,1,-2),(2,0,-1),(1,1,-1),(0,2,-1),(2,-1,0),e_1,e_2,e_3.
\]
The height parameters can be normalized to \((a,b,1)\). The two quadratic interpolation consistency equations from §7, up to a nonzero factor, are
\[
 (a-b)^2(b-1)(2a+b-1)=0,\quad
 (b-1)(2a+b-1)(2a^2-2ab-2a+b^2+1)=0 .
\]
Here \(2a-1>0\), \(b>0\). If \(b\ne1\) the first equation gives \(a=b\) and the last factor of the second becomes \((a-1)^2\), a contradiction. When \(b=1\), rows 3, 6, 7 are collinear at the same height one, impossible by the fixed-height parabola. All three classes are excluded.

**PROVED using the COMPUTED census.** For both rank-four survivors the formal points, in the nullspace basis, are exactly \((a,b,c,1-a-b-c)\) for \(a,b,c\in\{0,1\}\), differing only in order. Thus the moment points are an affine image of the vertices of a three-cube. Its six four-subsets formed by taking any two opposite pairs and its two parity tetrahedra have the same cardinality and centroid. Their eight indicator rows give eight polynomials with equal first three moments, hence a critical cluster. Conversely the Part I unique eight-row code gives precisely twelve different four-supports of its pair differences, each occurring twice as an edge, and four other edges of weight eight. The rational cube example in Part I realizes twelve up-to-sign primitive divisors. Thus the complete count is 12, not 28.

**COMPUTED — reproduction.** Compile with `c++ -O3 -std=c++17 ii_min_flat.cpp -o ii_min_flat`. Run `./ii_min_flat 12 > ii_min_cases.txt`, then `python3 -B ii_min_verify.py`. Python uses SymPy as in Part I. [Flat enumeration](ii_min_flat.cpp) · [Affine classes](ii_min_class.py) · [Plane interpolation](ii_min_plane.py) · [Checks](ii_min_verify.py)

**PROVED using a second COMPUTED census — allowing repeated factors.** The maximum number of *all integral kernel vectors of weight four*, including midpoints, is again 12 up to sign, with the same equality condition. On the 126-row run the numbers of processed nodes at ranks 0 through 5 are 1, 126, 5,887, 124,292, 1,086,867, 3,375. Ten rank-four nodes and the 3,375 rank-five nodes reach threshold. After affine classification there are 59 planar possibilities. Computing their actual number of midpoint and four-index relations without presuming that the displayed ordering is feasible gives
\[
\begin{array}{c|rrrrrrr}
\#\text{ relations}&12&13&14&15&17&18&19\\
\#\text{ shapes}&30&15&9&1&2&1&1 .
\end{array}
\]
The positive-product test of §7 excludes 56 of them. The three exceptions are excluded as follows, with row indices starting at zero and the same normalized height parameters as above.

* For \((0,8,19,23,33)\), the two excess interpolation equations have factors \((a-1)(a-b)^2(2a+b-1)\) and \((a-1)(a+b-2)^2(2a+b-1)\), with \(2a+b-1=y_4+b>0\). They force \(a=1\), and the midpoint triple \(4,5,7\) then has constant height.
* For \((6,13,20,24,27)\), an interpolation equation gives \((a-1)^2(b-1)(a+b-1)=0\) with \(2(a+b-1)=y_0+1>0\). At \(a=1\) the midpoint triple \(3,5,7\) has constant height, and at \(b=1\) the triple \(4,6,7\) does.
* For \((3,14,21,29,35)\), the first equation gives \((a-b)(a-5b+2)(a-2b+1)^2=0\), with \(a-5b+2=-y_0-1<0\). At \(a=b\) the second gives \((b-1)^3(2b-1)=0\). At \(a=2b-1\) it gives \((b-1)^3(3b-1)=0\). Since \(2b-1=y_2>0\), both force \(b=1\), at which the midpoint triple \(2,6,7\) has constant height.

Each constant-height triple contradicts the parabola argument. The ten rank-four records consist of the two cube records above and eight containing seven distinct points of the form \(C,C\pm U,C\pm V,C\pm(U+V)\). The latter configuration cannot occur on the moment surface. Translating the real part of every root and dilating all roots by a positive constant act linearly on moment triples and can put \(C=(1,0,-1)\). If \(C\pm(u,v,w)\) are on the surface, substitution into \(3V^2=YW+Y^4\), subtraction and addition give \(w=-3u-4u^3\) and \(v^2=u^2(1-u^2)\). For displacements \(U,V,U+V\), with height components \(u,u',u+u'\), the cubic identities force \(uu'(u+u')=0\). The displacement with zero height then vanishes in all three coordinates, a contradiction. Thus only the cube records survive for the extended count as well. Run `./ii_min_flat 12 all > ii_min_allcases.txt` followed by `python3 -B ii_min_all_verify.py`. [Full weight-four verification](ii_min_all_verify.py)

**PROVED using a COMPUTED census — the first gap.** An admissible configuration with eleven ternary directions already has twelve. Repeating the 70-row enumeration at threshold eleven produces five rank-four records, each of whose formal quotients carries twelve parallelograms. At rank five, only the classes represented by row-index sets `(0,5,12,17,30)` and `(1,10,17,24,29)` have exactly eleven parallelograms. The positive-product certificate from §7 rules out each. Compute the extra records with `./ii_min_flat 11 > ii_min_gap.txt` and verify them with `python3 -B ii_min_verify.py --gap`.

**PROVED using a COMPUTED instance — a different question about shortest vectors.** Take eight roots of the form \(x+i\), with \(x\in\{1,5,8,10,12,16,22,23\}\). Equal moment sums are precisely equal-cardinality multisets with the same sum and sum of squares. Among the multisets of cardinality \(k\), the number of unordered pairs of distinct multisets with matching signatures is \(0,0,0,0,0,0,13\) for \(k=1,\ldots,7\). The exhaustive script groups multisets by their two integer moments. The first collisions have disjoint supports because a common element would give a smaller collision. Thus this kernel has minimum 14 and thirteen shortest directions up to sign. The sharp twelve count applies to the requested degree-four shell, or to actual shortest vectors conditional on the lattice attaining weight four. [Multiset certificate](ii_min_longer.py)
