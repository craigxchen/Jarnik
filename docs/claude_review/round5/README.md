# Round 5/6: fake circles at `log R = O(M log M)` (snapshot, refereeing in progress)

**Status of this snapshot (2026-10-09).** `construct.md` has been refereed (`referee_construct.md`:
survives, with the corrections below). `lower.md` and `round6_sharp/sharp.md` are drafts whose
adversarial referees have not finished; their statements below are the authors' claims. I
independently re-derived Lemma S of `lower.md` (the Paley `l^1`-stability lemma) and found no error.
Nothing here proves the uniform theorem or a general growth improvement.

## The question

Let `W_min(M, C)` be the least `W = log N` over *angle models* (fake circles): genuine primes, a
genuine exponent profile and real angles, satisfying every gap inequality and residue strengthening
that is valid for actual circles. Actual clusters satisfy `3 M log M - O(M) <= W` (`two_thirds`).
If `W_min/(M log M) -> infinity`, actual circles would satisfy `M = o(log R/loglog R)`.

## Results

* **Framing correction (both write-ups, refereed).** The literal strengthening of `outside.md` A.3
  (factor `Q_M = lcm(1..2M)` for every character, pairs included) is not valid for actual circles. It
  fails already on `x^2 + y^2 = 5`, and by itself it forces `W >= (8 - o(1)) M^2`. The valid class uses
  each character's *actual* collision modulus.
* **`construct.md` (refereed: survives).** Paley-I matrices have a *skew pairing*, and the fully
  flipped profile with `b` copies has exact character margins (`>= b - 4` on pairs, `>= 2n + b - 8` on
  other `+-1` characters). This is item 392 of the notes restated via rank-one rectangles. With a
  residue design from auxiliary Gaussian primes it gives angle models with point-level residue data at
  `W <= (27 + o(1)) M log M + 9 M log(1/C)` (`Q_0 = 2M`; `9(2 + max(K,1))` for `Q_0 = M^K`). So these
  inputs cannot prove `M <= (2/27 - eps) log R/loglog R`.
  Referee corrections:
  * **(major, scope)** The barrier covers *point-level* residue collisions only. Actual circles also
    obey a prime-level parity law: point residue ratios lie in square-cosets fixed by the Legendre
    symbols `(p_j/q)`. The referee proved it and checked it on 211,136 actual instances. The construct's
    `M = 12` fake violates it.
  * Axiom (R) must be restricted to `v(c) != 0`. Counterexample on `N = 65`: `1+8i, 8+i, 4+7i, 7+4i`.
  * The `(S|)` data must be in the class for the `3 M log M` lower bound.
  * `q = 2` should be included.
  * Item 394 lacks the diagonal gaps.
* **`lower.md` (referee pending).** Lemma S: for prime `q = 3 mod 4`, `q >= 43`, every zero-sum integer
  `c` other than a pair, a signed column or a half-sum has `||H^T c||_1 >= (17/16) M`. Fakes at
  `(48 + o(1)) M log M`, and `(16 + 32A + o(1)) M log M` with residue data up to `M^A`. Its §8 isolates
  the prime-level parity obstruction (Prop. 8.1).
* **`round6_sharp/sharp.md` (referee pending).** Residue data that factor through prime residues
  with the genuine norm constraint. The forced information is exactly one parity bit per prime and
  modulus, `ell_j = [(p_j/q) = -1] mod 2` (Prop. 2.1). A design with classes inside parity classes
  (max pair demand `O(log M)`) gives:
  * **P3-char** (characters strengthened): `W <= (b max(A,2) + o(1)) M log M`, `b >= 33`, moduli
    `<= 3 M^A`. So the prime-level parity law does not give a growth improvement when it is used
    through characters.
  * **P3-all** (every monomial strengthened): only `W = O(M log^2 M / loglog M)`. **Open:** is
    `W_min = O(M log M)` here? A negative answer would be a growth improvement for actual circles.
  * Moduli up to `N^A`: not proved.
  * Appendix A, Sylvester multi-block plus Green-Sanders: a sketch only.

`lead_notes/` contains the lead's notes and checks that prompted round 6. `checks/` and
`referee_construct_checks/` hold the scripts; each runs from its directory.
