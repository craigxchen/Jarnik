# Why short-circle phase updates do not inherit the Gaussian moat information cost

The comparison is with [*Bounded-Step Walks on Gaussian Primes*, Sections 5, 7, and 8](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/Bounded-Step-Walks-on-Gaussian-Primes-September-26-2026/paper.pdf). That proof updates residues additively along a walk with a fixed Euclidean step bound, then uses avoidance of zero modulo selected Gaussian factors to charge information to short increment words. Here all points have one fixed norm and lie on an arc of length `C sqrt(R)`. The statements below identify what survives when the update is replaced by a phase, oriented area, or anchor index. They do not improve the general endpoint point-count bound. All entropy logarithms in this note have base two.

## 1. Exact one-scalar update and its entropy cost

Let `S` be `M` distinct Gaussian integers of modulus `R` in an arc of angular width less than `pi/2`. For `Y,Z in S`, put `D(Y,Z)=Im(bar(Y) Z)`. Since `Re(bar(Y) Z)>0` and `|bar(Y) Z|=R^2`, there is an exact update

```text
Z = [sqrt(R^4-D(Y,Z)^2)+i D(Y,Z)] Y/R^2.         (1)
```

Thus the full oriented area is a Markov update letter with no side bit. For each fixed `Y`, the map `Z -> D(Y,Z)` is injective on `S`: its value is `R^2 sin(theta_Z-theta_Y)`, and sine is strictly increasing over the available angle differences. If `Y,Z` are independent uniform points of `S`, then

```text
H(D|Y)=log M,       H(D)>=log M.                  (2)
```

For a uniform ordered distinct pair the same argument gives `H(D|Y)=log(M-1)` and `H(D)>=log(M-1)`. Equation (2) is a direct entropy obstruction to treating the full area as a cheap displacement while the destination is sampled broadly.

More generally, let `F(Y)` be any deterministic state, perhaps a quotient of all prime residues, and suppose independent uniform `Y,Z` admit an exact update `F(Z)=U(F(Y),T,B)`, where `B` has at most `2^b` values. Independence gives

```text
H(T)+b >= H(F(Z)|F(Y)) = H(F(Z)).                (3)
```

Thus a state with the near-maximal residue entropy required by the source cannot be updated from an independent new anchor using a low-entropy residual and `O(1)` extra bits. Conditional sampling can correlate `Z` with `Y`; (3) then becomes the exact lower bound `H(T)+b>=H(F(Z)|F(Y))` and need not be large.

The small endpoint determinant residual is `t_ij=D_ij/G_ij`, with a pair-dependent divisor `G_ij`. Equation (1) uses the full `D_ij`. If a proposed update transmits `t_ij` and enough information to recover `G_ij`, its combined data must still meet the entropy bound (2) for independent uniform destinations. Sending only `t_ij` is not the update (1).

## 2. A cheap abstract anchor walk loses the zero test

Order any finite cluster `S={z_0,...,z_{M-1}}` and let `J` be uniform modulo `M`. At each step choose an independent fair bit `B_t in {+1,-1}` and set `J_{t+1}=J_t+B_t mod M`. This is an exact one-bit Markov update with a uniform stationary anchor, even when every `z_j` avoids zero modulo both factors of every split prime outside the common norm. Yet the entire bit word is independent of `J_0`, so

```text
I(F(z_J0); B_0,...,B_{ell-1})=0                  (4)
```

for every residue function `F`. A cheap index update and zero avoidance alone therefore do not yield the positive information rate proved in the source. Its passing-list test uses the stronger identity `F_pi(z+h)=F_pi(z)+h mod pi`: the same increment word tests every candidate starting residue by translation. The index walk has no such additive update. Its physical steps can also be as large as the arc diameter `C sqrt(R)`, so the source's fixed-step collision estimate does not apply.

## 3. Multiplicative phase updates make avoidance tautological

For a transition `Y -> Z`, the ratio `u=Z/Y` lies in `Q(i)` and has norm one. Fix a split prime `p=pi bar(pi)`.

* If `p` does not divide the common squared radius `R^2`, both `Y` and `Z` are units modulo each factor. The ratio `u` is a unit there, and the exact local update is `x -> u x`. Under any word of such updates, every candidate `x != 0` remains nonzero. A passing-list test based on avoiding zero retains all `p-1` nonzero candidates; under a uniform prior it removes only the entropy `log(p/(p-1))`, and under the actual unit-supported prior it removes none.
* If `p` divides `R^2`, every point has zero residue modulo at least one of the two factors. The source's set `A(P)`, which excludes zero modulo **both** factors, contains no point of this norm fiber. For squarefree `p` in the primitive norm, a quotient retaining only the prime-orientation choice has at most two values, rather than a nearly uniform `p`-element residue coordinate.

Thus the exact multiplicative Markov law either preserves zero avoidance automatically or lies outside the source sieve. It cannot reproduce the source's `Theta(log p)` candidate-list reduction at a selected prime.

This assertion concerns the passing test, not all information in the word.
If `X` is the actual starting residue, `W` the transition word and `E`
the passing event, then `P(E)=1` and exactly
`H(X|W,E)=H(X|W)`. It does **not** follow that `I(X;W)=0`.
The entropy loss `log(p/(p-1))` instead refers to an artificial uniform
residue prior with the word fixed or independent of that prior.

## 4. The primitive norm pins the global phase gauge

After removing the common Gaussian gcd of a cluster, its squared norm is split-supported. Write its split-prime allocations as

```text
z_i = unit_i prod_p pi_p^(a_ip) bar(pi_p)^(e_p-a_ip),
min_i a_ip=0,       max_i a_ip=e_p.              (5)
```

Let `u in Q(i)` have norm one, and suppose `u z_i` is a Gaussian integer of the **same norm** for every point of this primitive cluster. At `pi_p`, multiplication shifts every allocation by a fixed integer `c_p`. Integrality of the row attaining `a_ip=0` gives `c_p>=0`; integrality at the conjugate factor of the row attaining `a_ip=e_p` gives `c_p<=0`. Hence `c_p=0` for every split prime, including primes outside the norm. Norm one also gives zero valuations at inert and ramified primes. Therefore `u` is a Gaussian unit.

Conjugation and units give at most eight global symmetries. They can reduce anchor entropy by at most `log 8`. Reanchoring is different: it changes the rational chart of the same configuration and need not act as an integral rotation on every point. In fact the relative allocations from anchor `i` already recover its absolute allocation from (5):

```text
a_ip = -min_j (a_jp-a_ip).                       (6)
```

So quotienting the full relative prime-orientation profile by a common additive shift does not create an unpinned, entropy-saving state once the actual primitive cluster is retained. A quotient that forgets these primewise minima also forgets the named orientation data needed for residue tests.

These four calculations leave a possible but different research route: a correlated walk among actual circle points with an update alphabet much smaller than its residue entropy **and** a new informative test compatible with its update law. The Gaussian moat zero test, used with either an abstract index walk or multiplicative phase transitions, does not provide that test.
