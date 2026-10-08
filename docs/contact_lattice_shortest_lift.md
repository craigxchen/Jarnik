# The shortest oriented contact-lattice lift

The oriented contact lattice always has a determinant-one integral
frame, as proved in
[the contact-lattice rigidity note](oriented_contact_lattice_rigidity.md).
The issue is quantitative: the five-row construction needs a frame of
height `exp(4w+o(w))`, while an unrestricted CRT lift can be much larger.
This note reduces that shortest-lift question, for a fixed first column,
to the least representative of one ordinary integer CRT class.  It does
not prove that the representative is large.

## 1. The ten contacts reduce to four congruences

Use the notation and the exact grouping in Section 6 of
[the five-point lift note](five_point_cm_norm_lift.md).  Thus

```text
v_1=(1,0),       v_2=(r_2,b),       v_3=(r_3,d),
v_4=(r_4,s_4),

J_1=G_14G_15G_23,
J_2=G_24G_25G_13,
J_3=G_34G_35G_12.
```

The repeated source direction within each of the first three groups
combines by the Chinese remainder theorem.  The full ten-contact
lattice is exactly

```text
L={(z,w) in Z[i]^2:
       z=0                         mod J_1,
       r_2z+bw=0                  mod J_2,
       r_3z+dw=0                  mod J_3,
       r_4z+s_4w=0                mod G_45}.             (1)
```

The four Gaussian moduli have disjoint support.  Their product is the
determinant ideal

```text
Delta=J_1J_2J_3G_45,
log|Delta|=10w+o(w),       log Norm(Delta)=20w+o(w).     (2)
```

This makes the modular geometry explicit.  At a split prime power the
matrix `U=(z w) in SL_2(Z)` sends the one source direction represented
by the relevant `v_i` to the selected isotropic line determined by the
Gaussian orientation.  Different contacts occur at disjoint rational
prime powers, so no local prime sees two of the four rational source
directions.  This explains both the unconditional CRT existence theorem
and why the shared cross-ratios do not create an additional local
compatibility condition.

## 2. Fixing the first column leaves one integer shear

Fix a Gaussian integer `z=x+iy` with `gcd(x,y)=1` satisfying `J_1|z`, and
assume that it occurs as the first column of some vector of `L` with
`q(z,w)=Im(bar(z)w)=1`.  Choose one Bezout partner `w_0` with
`q(z,w_0)=1`.  Every integral partner is uniquely

```text
w=w_0+t z,       t in Z.                                (3)
```

Indeed, two partners have real determinant zero with the primitive
integer vector `(x,y)`, so their difference is an integer multiple of
that vector.

Retain the static full-profile hypotheses

```text
gcd(b,Norm(J_2))=gcd(d,Norm(J_3))=1.
```

They follow from pairwise Gaussian coprimality of the original three
anchors, but are now conditions on the fixed lattice data.  For every
determinant-one candidate with `J_2|(r_2z+bw)`, the element `z` is a
unit modulo `J_2`: if `pi|z,J_2`, then
`bar(z)w-zbar(w)=2i` makes `w` a unit modulo `pi`, while `b` is a unit,
contradicting `pi|(r_2z+bw)`.  The same argument applies at `J_3`.
Consequently every candidate satisfies

```text
gcd_G(bz,J_2)=1,       gcd_G(dz,J_3)=1.                 (4)
```

In general, put

```text
C_4=gcd_G(G_45,s_4z).
```

The fourth congruence is soluble only if
`C_4|(r_4z+s_4w_0)`, and, when it is soluble, its integer period is
`Norm(G_45/C_4)`.  This is the exact nonunit version of the reduction.
For a full-profile first column the overlap is negligible. First, every
soluble candidate satisfies the exact identity

```text
C_4=gcd_G(G_45,z^2)       up to a Gaussian unit.        (5)
```

Indeed, let `pi^e` exactly divide `G_45` and put `a=v_pi(z)`.  The
conjugate-primitive row `Q_4` has `v_pi(Q_4)>=e` and
`v_pi(bar(Q_4))=0`.  If `a<e`, the identity

```text
2i s_4=bar(z)Q_4-z bar(Q_4)
```

gives `v_pi(s_4)=a`; if `a>=e`, both gcd valuations are already `e`.
Thus

```text
min(e,v_pi(s_4z))=min(e,2a),
```

which proves (5).  The prescribed core of `z=Q_1` is disjoint from
`G_45`, so `gcd_G(G_45,z)` divides the residual correction `K'_1`.
Consequently

```text
log Norm(C_4)<=2log Norm(K'_1)=o(w).                    (6)
```

For an odd conjugate-coprime Gaussian integer `J`, the natural map

```text
Z/Norm(J)Z -> Z[i]/(J)                                  (7)
```

is an isomorphism.  Its kernel is `(J) intersect Z=(Norm(J))`, and the
two finite rings have the same order.  Substitution of (3) in (1)
therefore gives three ordinary integer residue classes:

```text
t == -r_2 b^(-1)-w_0z^(-1)       mod J_2,
t == -r_3 d^(-1)-w_0z^(-1)       mod J_3,
t == -((r_4z+s_4w_0)/C_4)((s_4z)/C_4)^(-1)
                                      mod G_45/C_4.      (8)
```

Each right side has an integer representative by the ring isomorphism
(7). The admissible-partner hypothesis ensures the nonunit fourth
congruence is soluble. Write the corresponding integer
classes as

```text
t==tau_2 mod Norm(J_2),
t==tau_3 mod Norm(J_3),
t==tau_4 mod Norm(G_45/C_4).                            (9)
```

The rational moduli are pairwise coprime, so (9) is one class

```text
t==tau mod M,
M=Norm(J_2)Norm(J_3)Norm(G_45/C_4),
log M=14w+o(w).                                        (10)
```

Thus the entire remaining determinant-one lifting problem for this
`z` is one-dimensional.

## 3. An exact formula for the shortest partner

Replace `w_0` by an integral shear so that

```text
e=Re(bar(z)w_0),       -Norm(z)/2 <= e <= Norm(z)/2.   (11)
```

Put `A=Norm(z)` and `theta=e/A`.  The determinant and scalar product
identity gives, exactly,

```text
|w_0+t z|^2=((e+tA)^2+1)/A
           =A(t+theta)^2+1/A.                         (12)
```

Define the shifted least-representative size

```text
rho(z)=min_(k in Z)|tau+kM+theta|.                     (13)
```

Then the least possible second-column modulus among all frames in `L`
with first column `z` is

```text
min |w|^2=A rho(z)^2+1/A.                             (14)
```

This is sharper than a lattice-volume estimate.  It includes the
actual congruence class, its prime-power depths, and the archimedean
choice of reduced Bezout partner.

In the five-row budget `|z|=exp(4w+o(w))`.  A required frame of height
`max(|z|,|w|)<=exp(4w+o(w))` can therefore exist only if

```text
rho(z)<=exp(o(w)),                                    (15)
```

even though its ambient modulus is `M=exp(14w+o(w))`.  Conversely, a
lower bound `rho(z)>=exp(delta w)` for every admissible
`z=J_1R_1`, `|R_1|=exp(w+o(w))`, would raise the frame height by the
same exponential factor.

The determinant-ideal rigidity theorem complements (14).  If one such
frame has height `exp(4w+o(w))`, every non-unit-equivalent second frame
has height at least `exp(6w-o(w))`, since the determinant of two lattice
vectors is divisible by `Delta`.  Equivalently, below the
`|Delta|^(1/2)=exp(5w+o(w))` independence threshold, (14) can select at
most one frame up to a Gaussian unit.

## 4. What the fixed fourth boundary value contributes

The normalized moduli coordinate of the fourth row is

```text
a=d det(v_2,v_4)/(s_4 Delta_23),
Delta_23=det(v_2,v_3).
```

Equivalently,

```text
bd r_4=s_4(d r_2-a Delta_23).                         (16)
```

At the `G_45` contact, the fixed elliptic moduli section meets the
boundary `a=b_0` at a fixed rational value `c_E`.  Away from its fixed
model primes, the full rational contact depth therefore replaces `a`
by `c_E` in (16) modulo `Norm(G_45)`.  If `bd` is not a unit, retain the
cleared equation (16); removing its exact gcd with `G_45` costs at most
`O(log T)=o(w)` and leaves the same leading modulus in (10).

This makes `tau_4` in (8) an explicit rational expression in the two
small anchor vectors and the fixed value `c_E`, together with the
common modular inverse `-w_0z^(-1)`.  It does not by itself bound the
least CRT representative.  At each prime only this one target direction
is imposed, and the surjectivity construction can place it in an
arbitrary `SL_2` residue class.

## 5. The remaining quantitative statement

The usual count of integral `SL_2` matrices of height at most `H` is of
order `H^2` up to logarithmic factors.  At `p^e`, prescribing the image
of one projective source direction has coset index
`p^e+p^(e-1)`.  Thus the exact global contact-coset index is

```text
I=Norm(Delta) product_(p|Norm(Delta))(1+1/p),
log I=20w+o(w).                                        (17)
```

For the last estimate, if the distinct prime divisors are ordered as
`p_1<...<p_r`, then `p_j>=j+1`; hence the logarithm of the extra product
is at most `sum_j 1/p_j=O(log(r+1))=O(log log(3Norm(Delta)))`.

The Gaussian ambient-lattice index is `Norm(Delta)`; it should not be
confused with this `SL_2` coset index.  The count suggests a typical
lift scale near `exp(10w)`, far above the required `exp(4w)`.  It is not a
deterministic lower bound: a congruence class may contain the identity,
and (13) describes precisely the analogous exceptional event here.

A rigorous continuation would have to prove that the classes `tau`
arising from the oriented CM allocations and the five shared row
factors cannot satisfy (15), except in a controlled exceptional family.
Neither the fixed values `infinity,0,1,c_E`, the abstract Hermitian
normal form, nor CRT surjectivity proves this.  The exact new target is
the least-representative problem (13) as `z=J_1R_1` ranges over its
height-`w` residual factor.

The exact finite checks in
[check_contact_lattice_shortest_lift.py](check_contact_lattice_shortest_lift.py)
verify a full shear period, a nonclean correction-overlap period and
identity (5), and the shortest-partner formula (14).

## 6. Changing the median anchor triple

Choosing another three-row median does not give a second short vector
in the same lattice.  Let `A,A'` be two anchor triples.  Before primitive
coordinate normalization, let their Gaussian rotation divisors be
`g_A,g_(A')`, so that

```text
Q_i^A=P_i bar(g_A)/c_i^A,
Q_i^(A')=P_i bar(g_(A'))/c_i^(A'),                     (18)
```

with positive integer coordinate contents `c_i`.  Up to the individual
positive rational row scales, the two direction configurations differ
by multiplication with the primitive Gaussian representative of the
projective class

```text
rho_(A,A') ~ bar(g_(A'))g_A.                           (19)
```

Let `C(rho)` denote its two-by-two real multiplication matrix.  If

```text
Q_i^A=U_Av_i^A,       Q_i^(A')=U_(A')v_i^(A'),
```

then the exact common projective change of coefficient models is

```text
T_(A,A')=U_(A')^(-1) C(rho_(A,A')) U_A,                (20)
v_i^(A') is a positive rational multiple of T_(A,A')v_i^A.
```

The displayed matrix is integral.  Before removing any scalar content,
it satisfies

```text
det T_(A,A')=Norm(rho_(A,A')).                         (21)
```

After primitive projective normalization, the relative rotation in
(19) keeps precisely the cut blocks on which the two Boolean median
functions differ.  For five rows their number is

```text
8       if |A intersect A'|=2,
12      if |A intersect A'|=1.                         (22)
```

Thus its logarithmic projective height is respectively `4w+o(w)` or
`6w+o(w)`.  If both median frames have height `exp(4w+o(w))`, the raw
matrix in (20) has the safe upper bounds `exp(12w+o(w))` and
`exp(14w+o(w))`.  Its primitive height is at least `exp(4w-o(w))`:
choose a row in `A\A'`.  Its primitive coefficient vector has
subexponential height in the `A` frame and height `exp(4w+o(w))` in the
`A'` frame, so an integral primitive representative of (20) must have
at least that height.

Most importantly, (20) is not a determinant-one identification.  It
satisfies

```text
U_(A') T_(A,A')=C(rho_(A,A'))U_A,                      (23)
```

whose determinant is `Norm(rho_(A,A'))`, generally exponentially
large. Applying `C(rho)^(-1)` to the whole right side of (23) restores
exactly the old integral frame `U_A`, and produces no second frame.
Applying it to `U_(A')` alone generally loses integrality and still
leaves the source model changed by `T`. The contact coefficient vectors and their selected
Gaussian orientations have changed at the same time.  Therefore two
comparable short median frames belong to different integral lattice
models; determinant-ideal uniqueness in one `L` cannot compare them.

This exact transformation also explains why reanchoring has not
improved (15): it exchanges the exceptional short CRT representative
for a projective model change of height at least `exp(4w-o(w))`.  No
second-frame contradiction follows without a new integral comparison
that controls this determinant and scalar content.
