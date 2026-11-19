### A monochromatic four-cycle on an arc shorter than the critical scale

**PROVED.** There are genuine, one-unit Gaussian clusters of four points of arbitrarily large gcd-primitive norm \(N\) and angular width \(O(N^{-10/39})\) in which all four edges of a \(K_{2,2}\) have primitive cofactor **four**. Their four separation divisors are distinct. This rules out a monochromatic four-cycle exclusion even after shortening by a power of the parent norm.

*Proof.* For positive \(n,t\) write \(u(n,t)=(n+it)/(n-it)\). Choose a split prime \(M>61\) and \(j\) with \(j^2=-1\pmod M\). Put \(v=10M,\ u=20M,\ p=u^2+1,\ q=v^2+1\). These last two integers are coprime since \(p-4q=-3\) and \(3\nmid q\). Neither is divisible by \(2,5,M\). Choose \(x\in[L,2L)\), \(L=10pqM\), by the following simultaneous congruences.
\[
 x\equiv1\pmod2,\quad x\equiv1\pmod5,\quad x\equiv2+2u\pmod p,
 \quad x\equiv-2+2v\pmod q,\quad x\equiv-4/3+2j\pmod M.
\]
In particular \(x\asymp M^5\), with absolute constants. Set \(s=x^2,\ d=s^2+64,\ K_\pm=(x\pm2)^2+4\), so \(K_+K_-=d\). The first two marks on the rational unit circle are \(1\) and \(u(s,8)\).

For each of \(r=5,5p/q\) let \(k=rK_+\), \(l=K_-/r\). They are integers because \(5p\mid K_-\) and \(q\mid K_+\). They are both one modulo four by the same formulas. Put
\[
 a=(k+s)/2,\qquad b=(l+s)/2,\qquad c=(k+l+2s)/4 .
\]
Then \(a,b\) are positive odd integers, and expansion using \(kl=s^2+64\) gives
\[
 a^2+16=kc,\qquad b^2+16=lc,\qquad a+b=2c,\qquad ab-16=sc .
\]
Consequently \((a+4i)(b+4i)=c(s+8i)\) and \(u(a,4)u(b,4)=u(s,8)\). Add the mark \(u(a,4)\) for each \(r\). Its ratios to the first two marks are \(u(a,4)\) and \(u(b,4)^{-1}\), respectively. The Gaussian representatives are primitive of odd norm because \(a,b\) are odd and the imaginary coordinate is four. Their canonical Gaussian primary units are real signs, hence all four marks are in one unit class with primary prime generators. The two \(a\)'s increase strictly with \(r\), the two \(b\)'s decrease, and both \(a\)'s exceed both \(b\)'s for large \(x\), giving distinct separation divisors. The added marks are both between the original endpoints, so the angular width is \(2\arctan(8/s)\le16/s\).

To recover the size of the integer parent, clear the least common Gaussian denominator of these four unit-circle numbers. This gives an integer cluster of common, gcd-primitive norm \(N\). At each split prime the parent exponent is the range of the signed valuations of the four marks at one of the two places. For one added mark, the three pair distances in these integer valuations are \(v_P(d),v_P(kc),v_P(lc)\). For three integers their range is half the sum of their pair distances, here \(v_P(d)+v_P(c)\). Thus adding the two marks enlarges the width of the original valuation interval by at most \(v_P(c_1)+v_P(c_2)\), where \(c_1,c_2\) correspond to the two ratios \(r\).

There is a shared factor saved at \(M\). Both values of \(r\) equal five modulo \(M\), and
\[
 4rc=((r+1)x+2(r-1))^2+4(r+1)^2.
\]
The last congruence on \(x\) gives \(M\mid c_1,c_2\). Also \(M\nmid d\). Indeed at the same residue the norms from \(\mathbb Q(i)\) of the two values of \(K_\pm\) at \(-4/3+2i\) are \(16\cdot37/81\) and \(400\cdot61/81\), respectively. The two \(a\)'s agree modulo \(M\), so the factors \(a+4i\) at the two marks are divisible by the *same* Gaussian prime over \(M\); their conjugates are coprime to it. Hence the valuation intervals extend in the same direction there, saving at least one power of \(M\). Since \(5\le r<20\), this proves
\[
 x^4\le d\mid N,\qquad N\le d c_1 c_2/M\ll x^8/M\ll x^{39/5}.
\]
The lower bound ensures that the parents grow, and the upper bound gives the claimed angular exponent.

**COMPUTED.** For \(M=1000000009\), the certificate constructs the four Gaussian integers on a 398-digit primitive parent, finds exactly \(d c_1c_2/N=5000000045\), and verifies that all four reduced cofactors equal four. It certifies \(16N^{1/4}/s<1/9\) by integer arithmetic. [Integer certificate](iii_tri_cycle.py)

**PROVED — the obstruction to adding further neighbors in this family.** Fix distinct positive rational parameters \(r_1,\ldots,r_n\), with \(k,l,a,b,c\) integral and of the required parity along the specialization and \(a,b\) coprime to four. The preceding formulas put \(K_\pm,c_1,\ldots,c_n\), up to fixed nonzero scalars, in the list of quadratics \((x+h)^2+4\) with shifts \(h=\pm2,\ 2(r_i-1)/(r_i+1)\). The shifts are distinct, so these quadratics are pairwise coprime over \(\mathbb Q[x]\). Their integer values have pairwise gcd bounded in terms of the parameters after clearing fixed denominators, by their nonzero resultants. Since \(d\) and every \(c_i\) divide the normalized parent, \(N\gg_{r_1,\ldots,r_n}x^{4+2n}\). The endpoints have angular separation asymptotic to \(16/x^2\). With \(n\ge3\) the normalized arc constant grows at least as fast as \(x^{(n-2)/2}\). This obstruction concerns fixed ratios; in the four-cycle construction one ratio varies with \(M\), which makes the CRT gain effective.

**PROVED — a quantitative version for varying ratios.** If \(r_i\in[2,20]\) are distinct reduced fractions of height at most \(H\), write \(F_i(x)=((n_i+d_i)x+2(n_i-d_i))^2+4(n_i+d_i)^2=4n_i d_i c_i(x)\) for \(r_i=n_i/d_i\). All coefficients are \(O(H^2)\). The two base quadratics have bounded coefficients. Any pair of the resulting \(n+2\) polynomials has nonzero integral resultant \(O(H^8)\) by its four-by-four Sylvester determinant. At an allowed specialization the gcd of any two of \(K_+,K_-,c_1,\ldots,c_n\) divides the corresponding resultant, by the integral Bézout identity. Every listed value is \(\gg x^2\) uniformly. The product divided by their lcm is bounded by the product of their pairwise gcds, seen prime by prime. Thus for fixed \(n\ge3\) and \(H=x^{o(1)}\) one still has \(N\ge x^{4+2n-o(1)}\); the same extension fails at the critical scale. In particular three neighbors in this pencil would require rational parameter heights of power size in \(x\). This leaves general monochromatic \(K_{2,3}\) and longer cycles open.
