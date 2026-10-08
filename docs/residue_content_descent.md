# Primitive residue content and the cost of removing it

The uniform endpoint bound is not proved here. This note establishes
an exact content identity for every row subset, and shows by an
explicit endpoint family that removing even bounded common residue
content can increase the least possible circle radius by an unbounded
factor. The results are prose proofs.

Throughout, points are taken in one common Gaussian-unit class, as in
`balanced_bonus_phase_audit.md`. This permits primitive pair numerators
without additional quarter-turn units.

## 1. Primitive half-angle rows and their exact minors

Fix an anchor z_0. For every other point choose the conjugate-coprime
Gaussian numerator h_i such that

    z_i/z_0=h_i/conjugate(h_i),     h_i=u_i+i y_i,

and put h_0=1. These numerators are obtained directly from signed
Gaussian valuations. In particular,

    |y_i|=t_0i,

where t_ij is the positive primitive chord residue. Signs of individual
h_i do not affect the assertions below.

For i!=j put n_ij=N(gcd_G(h_i,h_j)). Then

    h_i conjugate(h_j)/n_ij

is conjugate-coprime and represents z_i/z_j. Consequently

    |det((u_i,y_i),(u_j,y_j))|=n_ij t_ij.                 (1)

To verify primitiveness, factor out g=gcd_G(h_i,h_j).
The remaining rows a,b are coprime to one another, each to its own
conjugate, and a conjugate(b) is coprime to conjugate(a) b.
The canceled factor is g conjugate(g)=N(g), exactly as in (1).

Thus the imaginary parts are already bounded by the total residue
product T. The remaining large factors in the minors are the
pair-dependent integers n_ij.

## 2. Gcd of a star equals gcd of every residue

For any three distinct indices a,i,j,

    gcd(t_ai,t_aj) divides t_ij.                          (2)

Choose a as anchor and write its primitive numerators as
h_i=u_i+i y_i and h_j=u_j+i y_j. If a prime power q^r divides
both y_i and y_j, coordinate primitivity makes both u_i,u_j
q-adic units. Hence q divides neither N(h_i) nor N(h_j), so
n_ij is a q-adic unit. The determinant in (1) is divisible by
q^r, and therefore t_ij is divisible by q^r. This proves (2),
including at q=2.

It follows that, for every subset I of at least two points and every
anchor a in I,

    gcd_(i<j in I) t_ij = gcd_(j in I, j!=a) t_aj.         (3)

One divisibility is immediate; the other follows from (2). In
particular, the answer does not depend on the anchor.

Equivalently, at every rational prime q, the residue valuations obey

    v_q(t_ij) >= min(v_q(t_ia),v_q(t_ja)).                 (4)

Thus divisibility by q^r defines an equivalence relation on the rows
after the diagonal is included. This holds for every prime, not only
for primes outside the conductor.

## 3. Exact Smith-index identity for every row subset

More generally, take any conjugate-coprime Gaussian rows h_i whose
quotients h_i/conjugate(h_i) are distinct. For a subset I of size
at least two, define

    G_I=gcd_G(h_i : i in I),
    delta_I=gcd_(i<j in I) |det(h_i,h_j)|,
    g_I=gcd_(i<j in I) t_ij.

Then

    delta_I=N(G_I) g_I.                                 (5)

Divide all rows in I by G_I. This divides every real two-by-two
determinant by N(G_I), preserves conjugate-coprimality of each row,
and leaves every primitive pair residue unchanged. It therefore
suffices to prove (5) when the common Gaussian gcd is a unit.

In that case delta_I is coprime to every N(h_i). To see this,
suppose an odd rational prime q divides delta_I and some N(h_i).
Choose a Gaussian prime pi above q dividing h_i. Since the common
Gaussian gcd is a unit, some h_j is not divisible by pi. Also
conjugate(h_i) is not divisible by pi. But the identity

    2i det(h_i,h_j)
       =conjugate(h_i) h_j-h_i conjugate(h_j)

would imply pi divides conjugate(h_i) h_j, a contradiction.
The prime above 2 and inert primes cannot divide N(h_i), because
h_i is conjugate-coprime. Thus delta_I is coprime to every n_ij
in (1), and (1) implies delta_I divides every t_ij. The reverse
divisibility follows from (1), proving (5).

Each normalized row h_i/G_I has coprime integer coordinates.
The integer row lattice they span consequently has Smith invariants

    (1,g_I).                                            (6)

Its index is entirely controlled by primitive residues. This holds
simultaneously for every subset I.

For the anchored collection with h_0=1, the entire row lattice is
particularly explicit:

    span_Z{h_i}=Z(1,0)+Z(0,g),
    g=gcd_i |y_i|=gcd_(i<j) t_ij.                       (7)

This gives a precise content reduction, but says nothing by itself
about the heights of the necessary changes of basis or their effect
on the circle's least possible radius.

## 4. Pair-dependent content remains after the global reduction

Using rows 1 and h_b as a rational basis gives

    h_i=(det(h_i,h_b)/y_b)·1+(y_i/y_b) h_b.

The second coefficient has denominator and numerator controlled by
T. The first still contains n_ib t_ib/y_b. Removing the global
Smith index does not remove these pair-dependent integers.

Rational rescaling of individual rows, or a projective transformation
of all rows, does not remove their cross-ratio contribution. The
exact cross-ratio factorization and its forced growth in a uniform
profile are already treated in `uniformity_audit_parallel.md`:
they should not be inferred to have bounded height from (5).

Also, dividing a common Gaussian row factor rotates the half-angle
rows. It need not preserve their small imaginary coordinates.
Conversely, dividing their common imaginary-coordinate gcd is
usually a nonconformal real linear map, whose radius cost must
be computed rather than inferred from its small coefficients.

## 5. A bounded-content normalization can destroy endpoint scale

Here is an explicit instance of that last problem. For positive
integers t and j=0,1,2,3, put

    H_j=2(t+j)+i,
    z_j=H_j product_(l!=j) conjugate(H_l).

Each H_j is conjugate-coprime. For distinct indices their Gaussian
gcd divides 2(i-j), and has odd split-prime norm. Since |i-j|<=3,
this gcd is a unit. It follows that all primitive residues are

    t_ij=2|i-j|,
    T=768,
    gcd_(i<j) t_ij=2.                                   (8)

The original circle points have common Gaussian gcd a unit:
at any Gaussian prime, at most one H_j is divisible by it, and
the corresponding point with the opposite orientation has zero
valuation. Thus no common rational Gaussian multiplier can lower
their radius. Their least possible radius is

    R=product_(j=0..3)|H_j| ~ 16 t^4.

Their angular span is

    Delta=2[arctan(1/(2t))-arctan(1/(2(t+3)))]
         ~ 3/t^2,

so Delta sqrt(R) tends to 12. This is a genuine endpoint family
with fixed residue product.

Anchoring at z_0 gives h_0=1 and

    h_j=H_j conjugate(H_0)
       =(4t^2+4jt+1)-2ji,       j=1,2,3.

The natural map removing the common imaginary content is the
fixed height-two projective map

    (X,Y) -> (X,Y/2).

It produces h'_0=1 and

    h'_j=(4t^2+4jt+1)-ji=B+jL,
    B=4t^2+1,       L=4t-i.

The Gaussian columns B,L are coprime: a common divisor divides

    4B-L(L+2i)=3,

and the inert Gaussian prime 3 cannot divide L, whose imaginary
part is -1. Therefore gcd_G(h'_i,h'_j) divides i-j.

The integer coordinates of each h'_j are coprime. For j=1,3,
both are odd, and exactly one factor 1+i must be removed to
obtain conjugate-coprime d_j. For j=2 no removal is needed.
For the coordinate gcd at j=3, the first coordinate reduces to
t^2+1 modulo 3, which is never zero.
The three resulting d_1,d_2,d_3 are pairwise Gaussian-coprime.
For the pair 1,3 the only possible common factor before removal
was 1+i; the other pairs already had unit gcd.

With the anchor d_0=1, the exact least-radius formula from
`fresh_algebraic_parallel.md` therefore gives

    R'_min=|d_1 d_2 d_3|
          =|h'_1 h'_2 h'_3|/2
          ~32 t^6.                                     (9)

The new angular span is

    Delta'=2 arctan(3/(4t^2+12t+1)) ~3/(2t^2).

Consequently,

    R'_min/R ~2t^2,
    Delta' sqrt(R'_min) ~6 sqrt(2) t.                   (10)

The map halves the angular span asymptotically, but its least
possible circle radius grows enough to destroy endpoint scale.
Its height, the removed content, and the original T are all fixed.

This is a four-point example, not a counterexample to uniformity.
Its purpose is precise: even for endpoint configurations, removing
T-controlled content by a T-controlled projective map has no
automatic radius cost bounded by a power of T.

The subset-content identity was checked on 2,600 subsets of randomly
chosen conjugate-coprime Gaussian rows. The displayed family was
checked at t=1,10,100,1000, including its exact residue product and
pairwise coprimality after normalization. These checks supplement
the proofs above.

## 6. What remains

The new unconditional statement is the subset hierarchy (3)--(6).
It removes the common integer content from the rank-two
reconstruction exactly. The remaining pair-dependent Gaussian
contents are not bounded by it.

No estimate R<=const·T^A has been proved for a fixed large
near-uniform tuple. Such an estimate would require a new use of
those remaining contents together with the short arc; neither
Smith normalization nor a bounded-height projective map provides
it automatically.

## 7. What the full residue hierarchy contributes to Pluecker equations

The ultrametric statement gives an exact way to account for the entire
hierarchy of residual contents. For each rational prime q and integer
r>=1, let P_(q,r) be the partition of the k rows defined by

    i equivalent to j  iff  q^r divides t_ij,

with diagonal pairs included. Then

    log T = sum_(q,r) sum_(C in P_(q,r))
                binom(|C|,2) log q.                       (11)

More generally, for 2<=h<=k,

    sum_(|I|=h) log g_I
       = sum_(q,r) sum_(C in P_(q,r))
                binom(|C|,h) log q.                       (12)

A subset I contributes at level r precisely when all its rows belong
to one partition class. Since

    binom(c,h)/binom(c,2) <= binom(k,h)/binom(k,2),

equation (12) implies the sharp bound

    product_(|I|=h) g_I <= T^[binom(k,h)/binom(k,2)].        (13)

Thus if log T=o(W) for a fixed k, all these residual Smith contents
together still have total logarithmic size o(W). The large factors
N(G_I) in (5) are not included in this conclusion.

This calculation has not produced a new positive residue exponent
when combined with the positive Pluecker equations. The hierarchy
restricts valuations of the t_ij; the equations also constrain their
unit values modulo conductor primes and the pair-dependent Gaussian
contents. Those simultaneous global constraints remain uncontrolled.

There is a useful caution about prime powers. The elementary bound

    sum_(q divides conductor and T) log q <= log T

controls radical conductor weight only. It does not control the
full weighted sum e_q log q. For example, with pi=2+i and H=2+5i,
the four points

    pi^e H, pi^e conjugate(H),
    conjugate(pi)^e H, conjugate(pi)^e conjugate(H)

have common norm 5^e·29 and primitive common Gaussian gcd one.
The two within-allocation pair residues are 5, while the four
cross-allocation primitive residues are not divisible by 5.
Thus v_5(T)=2 independently of e, although the conductor's
5-weight is e log 5. These points are not asserted to form endpoint
arcs; they show why the valuation-support comparison itself cannot
discard prime exponents.

No estimate R<=const·T^A follows from (11)--(13), and no integral
near-uniform endpoint configuration with subpower T has been
constructed. A formal model that declares pair contents without
realizing them as actual Gaussian gcds would not settle this issue.
