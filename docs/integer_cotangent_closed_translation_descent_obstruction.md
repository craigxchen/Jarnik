# A short actual clique with no radius-decreasing closed translation

The exact reciprocal-stretch grid for a triangular map is conditional on
having a map between two integral circle tuples.  An integer-cotangent
translation offers a concrete candidate: $X_i\mapsto X_i+t$ with the
same $L$.  The following four-point source lies on a short normalized
arc, but **every nonzero translation that preserves all pair cotangents
as integers increases its least squared radius**.  Thus the closed
translation operation has no universal one-step descent property, even
when its primitive matrix has determinant one.  This finite example
does not rule out a construction restricted to sufficiently many rows.

The source data are

\[
 L=1,\quad X=(157,182,447),\quad
 (Q_{12},Q_{13},Q_{23})=(-1143,-242,-307). \tag{1}
\]

All pair cotangents are integers.  The exact all-edge norm lcm is
$N=212{,}298{,}125$, and the exact normalized arc constant is

\[
 C_*=2N^{1/4}\arctan(1/157)\in(1.53,1.54). \tag{2}
\]

For a translation by $t\in\mathbb Z$, the update

\[
 Q'_{ij}=Q_{ij}+t(X_i+X_j+t)/(X_i-X_j) \tag{3}
\]

shows that closure is exactly the three divisibilities

\[
 25\mid t(t+339),\qquad
 290\mid t(t+604),\qquad
 265\mid t(t+629). \tag{4}
\]

Their ordinary prime-power decomposition is explicit:

\[
 t\equiv0\pmod2,\qquad
 t\bmod25\in\{0,11\},\qquad
 t\bmod29\in\{0,5\},\qquad
 t\bmod53\in\{0,7\}. \tag{5}
\]

At five, the alternatives modulo 25 already satisfy the two other
conditions modulo five.  Since the four moduli in (5) are coprime,
the eight closed classes modulo
$D=\operatorname{lcm}(25,290,265)=76{,}850$ are exactly

\[
 0,\ 17{,}550,\ 29{,}150,\ 31{,}436,\ 43{,}036,\
 60{,}586,\ 65{,}250,\ 72{,}186. \tag{6}
\]

The nearest nonzero class to zero is
$72{,}186-D=-4{,}664$.  Therefore any nonzero closed
translation satisfies $|t|\ge4{,}664$, including all translates of
these classes by arbitrary multiples of $D$.

One target pair already defeats radius decrease.  Its cotangent is

\[
 Q'_{12}=-\frac{(157+t)(182+t)+1}{25}. \tag{7}
\]

For $|t|\ge4{,}664$, both factors in the numerator have the same
sign and absolute values at least $4{,}507$ and $4{,}482$.  Since
closure makes (7) an integer,

\[
 |Q'_{12}|\ge\frac{4{,}507\cdot4{,}482+1}{25}
 =808{,}015. \tag{8}
\]

With $L=1$, the reduced Gaussian edge norm of an integer cotangent
$q$ is $(q^2+1)/\epsilon$, where $\epsilon\in\{1,2\}$.
Hence its norm is at least $q^2/2$.  The target all-edge least
squared radius consequently obeys

\[
 N_t\ge n(Q'_{12},1)\ge(808{,}015)^2/2
 >212{,}298{,}125=N \qquad(t\ne0\text{ closed}). \tag{9}
\]

The translation is represented in this source chart by the primitive
triangular matrix $\left(\begin{smallmatrix}1&t\\0&1\end{smallmatrix}\right)$,
so its determinant is one.  The obstruction is map existence under
integrality and radius decrease; the reciprocal-stretch grid theorem
does not assert that such a map exists.

The exact congruence classes, all-edge radius, normalized arc, and
pair-edge bound are checked in
[the companion checker](check_integer_cotangent_closed_translation_descent_obstruction.py).
The same source supports independent prime-cut orientation flips by
closed translations, but those flips incur the height cost above; see
[the orientation audit](integer_cotangent_closed_translation_prime_orientation_flip.md).
