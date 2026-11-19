## Two exact saturated exclusions in degree one

**PROVED with COMPUTED polynomial identities.** No admissible linear-block solution has just one or two distinct imaginary parts among its ten roots. The two-value statement allows arbitrary unequal heights, including heights of the same sign. The first four moments already suffice.

*Proof.* Write the two nonzero heights as $A,B$. The first moment forces the edges at either height to form a regular spanning graph. For a nonconstant profile on five vertices the degree must be two by parity, so the colours are complementary pentagons. The one-colour case is excluded directly: equal imaginary moments of orders $2,3,4$ make the first three power sums in the abscissas of the two residual triples of each pair equal. Newton's identities give identical triples, contradicting distinctness.

Write $S$ for a cyclic shift. The vertex incidence maps in the two colours are $1+S$ and $1+S^2$, with cyclic edge labels. Since $(1+S^2)^{-1}(1+S)=1-S^2+S^4$, changing the order of the cyclic labels gives $T=1-S+S^2$. For each $i$ let $t_i$ be the first colour's abscissa and $y_i$ the other colour's abscissa. Set $Z=Tt$, $\alpha=A/B$. The second moment gives $y_i=d-\alpha Z_i$ with $d$ common. The five $t_i$ are distinct, as are the five $y_i$ and hence the $Z_i$.

Put $D_i=(t_i-t_{i+1})(t_{i+2}-t_{i+1})$ and $F_i=t_{i+1}D_i$. The componentwise identities for powers are
$$
 T(t^2)=Z^2-2D,\qquad T(t^3)=Z^3-3DZ-3F .
$$
From the third moment $y^2+\alpha T(t^2)$ is constant. Consequently $D_i=u_0+u_1Z_i+u_2 Z_i^2$ with $u_1=-d$, $u_2=(1+\alpha)/2$. The fourth moment makes $y^3-B^2 y+\alpha T(t^3)-\alpha A^2 Z$ constant. Substitution of the quadratic gives
$$
 F_i=v_0+v_1Z_i-2u_1u_2Z_i^2+\tfrac13u_2(1-4u_2)Z_i^3 .
$$

Normalize $(t_0,t_1,t_2,t_3,t_4)=(0,1,x,y,z)$ by a real affine change. Interpolate $u$ on the first three $Z_i$ and $v$ on the first four $Z_i$ by Cramer's rule. Require the two remaining quadratic evaluations, the one remaining cubic evaluation, and the two specified upper cubic coefficients, giving five polynomials in $x,y,z$. All interpolation denominators are nonzero by distinctness. The script keeps the squarefree part of each equation and strips only linear factors dividing a known nonzero $t_i-t_j$ or $Z_i-Z_j$. Their degrees in characteristic zero are $5,5,9,9,10$.

Take the four resultants eliminating $z$ from the first equation and each remaining equation. For each irreducible factor $f$ of the first resulting bivariate polynomial take the gcd of its nonzero $y$-resultants against the other three. Their product gives a necessary univariate equation $P(x)=0$ of degree 32 after the same safe deletions. Reversing the cyclic order preserves the conditions because $T$ reverses up to a cyclic shift and $t_{i+1}$ is the middle of the triple. Normalize the triple $(t_2,t_1,t_0)$ to $(0,1,x')$, where $x'=x/(x-1)$. The exact gcd of $P(x)$ and $(x-1)^{32}P(x/(x-1))$ is $x-2$. Every cyclic triple must therefore be an arithmetic progression, impossible since five equal increments around a closed cycle sum to zero. The same calculation over the algebraic closures at $p=17,23,103,107$ gives the same contradiction with projected degrees $31,32,34,32$. Thus the saturated slice is empty, of degree zero, in these characteristics too. [Exact two-height certificate](p3_two_height.py)

**PROVED with COMPUTED polynomial identities.** No admissible linear-block solution has a vertical reflection compatible with the marked edges, even after relabelling. Here a vertical reflection through $\operatorname{Re}r=c$ means $r_{\sigma e}=2c-\bar r_e$.

*Proof.* Translate $c$ to zero. Since the roots are distinct, $\sigma^2=1$. The normalized sixth-moment labels $h_i$ are distinct and transform by $h_{\sigma i}=-h_i$ after subtracting their common centre. Consequently at most one mark is fixed. The only possibility is $\sigma=(12)(34)$, with $h=(1,-1,k,-k,0)$ after scaling labels. In particular $k\ne0,\pm1$.

For clarity here is the abscissa reduction. Let $x_{ij}=\operatorname{Re}r_{ij}$ and $s_i=\sum_{e\ni i}x_e$. The real degree-two coefficient of the triangle equation
$$
 \sum_{\rm cyc}(h_j-h_i)\beta_{ij}\bar\beta_{kl}\bar\beta_{km}=0
$$
gives a coboundary
$$
 2(h_j-h_i)x_{ij}+h_i s_j-h_j s_i=d_j-d_i.
$$
Indeed the terms have zero sum on every triangle after substituting $x_{kl}+x_{km}=s_k-x_{ik}-x_{jk}$. Set $(Df)_i=\sum_{j\ne i}(f_j-f_i)/(h_j-h_i)$. Summing the preceding relation divided by $h_j-h_i$ over $j$ gives $Dd=hDs-2s$. A polynomial $s(h)$ of degree at most four interpolates the five $s_i$ and is odd by symmetry. Hence $s=s_1h+s_3h^3$. Here $D(h)=4$, $D(h^2)=3h$, $D(h^3)=\sum_i h_i^2+2h^2$. Thus $d$ is proportional to $h^2$ up to a constant (the kernel of $D$ consists of constants, by triangularity with leading coefficient $5-j$ on $h^j$). The coboundary formula now gives
$$
 x_{ij}=(h_i+h_j)(1-vh_ih_j)
$$
after rescaling all real and imaginary parts by the original $x_{15}\ne0$. If $x_{15}=0$, then $r_{15}=r_{25}$.

Write $X=x_{13}$, $Y=x_{14}$, and $y_{ij}=\operatorname{Im}r_{ij}$. Put $y_{13}=ka$, $y_{14}=kb$ as the two free parameters. The first two moments and reflection give
$$
 C=y_{15}=-k(Xa+Yb),\quad D=y_{35}=-Xa+Yb,\quad
 y_{12}=C+2D-k(a+b),\quad y_{34}=2C+D-k(a+b).
$$
Reflection determines all other heights. For each edge write $y_e=p_e a+q_e b$ using these formulas. Put $z=(1,a^2,ab,b^2)^t$. The real coefficients in degrees one and zero of every triangle are linear forms in $z$. Use the script's four rows of these forms, denoted by $R$, and $w$ the four signed maximal minors of its first three rows. Only factors dividing $k(k^2-1)$ are removed in row and vector normalization. Necessarily $\det R=0$ and either $w=0$ or $w$ is a nonzero multiple of $z$.

The Veronese condition $w_1w_3-w_2^2=0$ follows in both cases. Additional homogeneous equations come from the imaginary moments $3,4,5$ for marks one and three relative to mark five. Multiply each by $a$ before homogenizing in $z$. For example at power three use the signed incident-edge sum of
$$
 (p_e z_1+q_e z_2)\,[3x_e^2 z_0-(p_e^2 z_1+2p_eq_ez_2+q_e^2z_3)] .
$$
All six polynomials must vanish at $w$, including when $w=0$. For each irreducible factor of $\det R$ depending on $v$, take the gcd of its nonzero $v$-resultants against the seven homogeneous substitutions. Keep all remaining factors depending on $k$ alone except the three known exclusions. This gives a necessary $P(k)=0$ of degree 42 over $\mathbb Q$. Swapping the two transposed vertex pairs gives $P(1/k)=0$. The exact reciprocal gcd is $k^4-3k^2+1$. Sign changes and that swap permit $k=(1+\sqrt5)/2$.

In the quadratic field at this last value, the gcd of $\det R$ and the seven substitutions in $v$ is $v-(5-\sqrt5)/4$. There the exact $w$ is
$$
 (0,\ (5+\sqrt5)/16,\ -(5+\sqrt5)/16,\ (5+\sqrt5)/16).
$$
This is nonzero and cannot be proportional to $z$ because $z_0=1$. Thus the saturated slice is empty. The same algebraic-closure result holds in characteristics $101,109,139,149$ by the identical calculation with both sign and reciprocal symmetries. In those four fields the golden ratio already splits. These statements are geometrically empty slice results, not dimension claims about the full system. [Exact reflection certificate](p3_reflection.py)

**PROVED — freeness of relabelling on the entire split nondegenerate locus in characteristic zero.** Over an algebraic closure write the two independent lists of roots as $u_e,v_e$. Any permutation of the labels fixing their class modulo simultaneous affine changes of source and an optional swap of these two lists is trivial. Hence, for instance, any isolated unlabeled solution contributes 120 distinct labelled solutions modulo source and global swap.

*Proof.* Let $\sigma$ have order $d$ and $\lambda$ be the linear part of the source change. The multiplier on the distinct centred sixth labels is $\mu=\lambda^6$ without a swap and $\mu=-\lambda^6$ with one. Its order equals $d$. Without a swap, $\lambda^d=1$ because the $d$-fold change fixes every root of a list. The possible nonidentity permutation orders of five marks are $2,3,4,5,6$, so only $d=5$ remains. With a swap and odd $d$, one has $\lambda^{2d}=1$ but then $\mu^d=-1$, impossible. With a swap and even $d$, $\lambda^d=1$. The only possibility is $d=2$, $\lambda=1$ or $-1$. For $\lambda=1$ the translation is zero in characteristic zero, and some edge fixed by $\sigma$ would have $u_e=v_e$. For $\lambda=-1$ this is the reflection slice just excluded. The preceding computations on that slice make sense over the split model with $x_e=(u_e+v_e)/2$ and $y_e=(u_e-v_e)/(2i)$.

For the remaining ordinary five-cycle absorb the translation in the source and choose a primitive fifth root $\zeta$. In cyclic vertex coordinates the two edge orbits are the pentagon and its complement. Their two lists have the form $u_{i,i+j}=A_j\zeta^i$, $v_{i,i+j}=B_j\zeta^i$, $j=1,2$ (one may choose either primitive root by the author's cyclic order). At mark $i$, independence of moments for $k=1,2,3,4$ implies
$$
 G_k(A_1^k-B_1^k)+(A_2^k-B_2^k)=0,\qquad
 G_k=(1+\zeta^{-k})/(1+\zeta^{-2k}).
$$
Let $m_j=(A_j+B_j)/2$, $\delta_j=(A_j-B_j)/2$. Normalize $\delta_1=1$, as disjointness gives $\delta_1\ne0$. Then $\delta_2=-G_1$ and $m_2=c m_1$, $c=G_2/G_1$. The third equation gives $m_1^2=(G_1^3-G_3)/(3(G_3-G_1c^2))=5-8\zeta^2-8\zeta^3$, of norm one. The fourth equation divided by $8m_1$ now has left side $G_4(m_1^2+1)-G_1c(c^2m_1^2+G_1^2)=-24\zeta^3-32\zeta^2-24\zeta$, a nonzero element of $\mathbb Q(\zeta)$. [Cyclotomic reductions](p3_fivecycle.py)
