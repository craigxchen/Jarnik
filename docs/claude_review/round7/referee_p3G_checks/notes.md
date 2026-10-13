# referee p3G scratch notes

Checked by hand (OK): Lemma 1.2; Lemma 2.3 (1)-(3) incl. sigma formulas; Lemma 2.4; Cor 2.6;
Thm 3.1; Thm 3.2; Prop 3.3 (rank bound); Cor 3.4; Lemma 5.1; Lemma 5.2 (closed form); Thm 5.3 seeds
(both cases), propagation (sigma_a v_a + sigma_b v_b = e_a - e_b in <E_0>), union bound numerics;
Prop 6.1 steps 1-5 (Haar argument uses outside.md A.1 connected coset; fine).

Candidate issues:
1. Sec 4.2 gloss "covers every b-copy Paley profile with b(M-1) <= r built on these primes":
   a member with b(M-1) < r columns has a DIFFERENT (smaller) prime set; closure of a subset is
   not the restriction of the full closure. Only b(M-1) = r is covered.
2. Headline "Route G does not produce fakes at O(M log M)" / "on prime sets of polynomial size the
   class contains no fakes at all": proved only for 6 explicit random sets + model (non-adaptive)
   + heuristic near-relation count. Adaptive selection = open point 8.4(b).
3. Sec 0 "A fake needs about M/2 columns of size ... (proved with c=1/2-eps in the model)" omits
   (a) near-relation hypothesis O(M/log M) (heuristic, Remark 3.5), (b) non-adaptivity.
4. Sec 7 "P3-G-char does not collapse" is not proved; only "this mechanism unavailable".
5. Heuristic 5.3 item 3 "independent set" count: avoiding certificates is weaker than independent
   set (only circuits matter) -> test adaptively.

## Run 1 (run_cases.py, rcore.py, rclosure.py): Theorem 4.1 re-run independently
All 6 collapse rows reproduced: same prime sets (pmin/pmax and tight counts S/T/P identical),
k<=2 closure complete, every certificate re-checked by exact checker (orientation via Im(A^2(-i)^s),
sin^2 via Im(A(1-i)^s)^2/(2^s|A|^2), lattice membership by sparse echelon over Z), e_k in lattice
for all k. Partial rows 636/643 and 1279/1283 also reproduced. Outputs cases_A.jsonl, cases_B.jsonl.
Cert counts differ from author's (89 vs 86 etc.) -- different firing order, irrelevant.

## Run 2: subset gloss (subset_test.py -> subset_out.jsonl)
Verified M=164,b=8,K=4 set (r=1304, complete). Subsets of 163 primes (b=1 Paley at M=164): 35-54
columns FREE in the k<=2 closure; subsets of 326: 1-9 free. So the gloss in 4.2 ("covers every
b-copy Paley profile with b(M-1) <= r built on these primes") is false for b(M-1) < r.

## Run 3: completion sum (completion_check.py) reproduces 15190.1 and 4051.9 exactly.

## Run 4: adaptive greedy (adaptive.py, eval_set.py -> adaptive_eval.jsonl)
M=164, r=163, X=328, greedy avoiding: same-v seeds (singles/doubles/pairs), n=(1,2) double-single,
and every closed walk (pairs+singles, double cover) of arc weight < pi (ties rejected too).
Got 163 primes at K=3.75 (252 cands), 3.5 (404), 3.25 (2360), 3.0 (79625 cands).
Evaluated: k<=2 closure EMPTY (163 free) for all; with tight triples added (18111 at K=3) still empty.
Random baselines same code: K=3 and K=3.5 fully pinned (collapse). Pair degree adaptive K=3: 46.6
(random 56.9). => adaptivity, not density, decides; the 4.4 threshold evidence is random-only.
Next: Farkas/LP certificate y with sigma_i <v_i,y> > 0 for all tight weight<=3 relations
=> NO Lemma-2.4 certificate of any k, n relative to {0} => restricted full closure {0}, exactly.

## Run 5: Farkas LP useless (ties & quantisation). KEY OBSERVATION (referee, proved):
For any zero combination sum n_i sigma_i v_i = 0: sum n_i sigma_i alpha_i = -(pi/4) sum n_i sigma_i s_i,
and sigma_i alpha_i = d_i mod pi, so sum n_i d_i in (pi/4)Z and the arc sum sum n_i (pi-2d_i) is in (pi/2)Z.
=> relative to Pi_0={0}, a certificate has arc EXACTLY pi/2; non-certificates have arc >= pi.
Checked numerically (lemma52_quant_out.txt): all same-v opposite-orientation combos give arc/(pi/2)
in {1,2,3} to 1e-13. Lemma 5.2 closed form exact for m<=7 (0 mismatches).
Author code on referee adaptive sets: 0 certificates, 163 free (author_code_on_adaptive.txt).

Minor: density "within 4%": 355/370 is -4.05%. p3data cite "D_max/H 3.2-3.7 at X=2M": p3data says
2.2-3.2 at X=2M, 2.9-3.7 at X=4M. Prop 6.1 constant: psi(X)~X so (66+o(1)) M X actually.

## Run 6: more adaptive
M=644 K=4.0: 643 from 2093 cands, k<=2 closure empty (pair deg 26.4; random same params collapses, Thm 4.1 row).
M=324 K=3.5: 323 from 9961, empty (pair deg 37.7).
Tie-accepting greedy (adaptive_ties.py; ties are not certificates by Lemma R.1):
M=164 K=2.75: 163 from 13456, pair deg 82.2 (random 87.7 -> collapse), 44279 triples, k<=2+w3 empty.
M=164 K=2.5: only 96 of 163 before the pool (~12.5k primes) ran out.
brute_k3.py (all level-0 certs with <=3 terms, n in {1,2}, incl doubles): random K=3.5: 43; adaptive K=3.5: 0; K=3.0: 0.

## CORRECTION (important)
eval_set "k2w3_free=163" was an ARTIFACT: rclosure.closure skips classes not representable in the
signed-component quotient (triples at level 0 have 3 roots). brute_k3 with w3 finds level-0
certificates in adaptive sets (K=3.5: 4; K=3.0: 33, incl. same-v triple pairs s=2,3).
=> Do NOT claim triples are evaded. New gclosure.py: arbitrary lattice (canonical coset reps via
echelon basis), k<=2 + k=3 (two terms share a column), checked by mpmath arc + lattice membership.
gclosure on random 164/4.0 (w<=2): 34 free (author k<=2: 41). adaptive 164/3.5 (w<=2): 163 free.

## Run 7: general closure with weight 3 (gclosure_w3*.jsonl, prune_*.txt)
adaptive (w<=2-greedy) M=164: K=3.75 -> 162 free; K=3.5, 3.25, 3.0 -> 0 free (triples seed, pairs propagate).
random M=164 b=1 w<=3: K=3.75 -> 0 free; K=4.0 -> 4 free; K=4.5 -> 160 free.
prune K=3.5: delete 4 columns -> 159 columns, no level-0 cert (<=3 terms, w<=3), gclosure empty.

## Run 8
prune K=3.0 (w<=3): 15 deleted -> 148 free.  M=1284 K=4.0 (ties accepted) 1283 from 3558; author k<=2 closure empty; random 1284/4.0 collapse.
brute_k3 w<=2: conservative sets (644/4, 324/3.5, 164/3.0, 3.5) -> 0; TIE-ACCEPTING sets 164/2.75 -> 4, 1284/4.0 -> 3
(certificates combining doubles 2e_j with pairs, k=3, not covered by the closed-walk check).
gclosure w<=2 on 164/2.75 ties set -> 0 free (4 k3 seeds + dense pairs). => caveat rows; pruning w2 running.

## Run 9 (final)
prune2 (w<=2) on tie-accepting sets: 164/2.75 -> 161 free (2 deleted); 1284/4.0 -> 1280 free (3 deleted).
Author code on 644/4.0, 324/3.5, 164/2.75 selected sets: 0 certs.
Aborted runs (killed, no output, files deleted): conservative 164/2.75, 324/3.0, 644/3.5, 1284/4.0.
Report: ../referee_p3G.md. Verdict: math correct; Thm 4.1 reproduced; headline not established;
issues A (scope), B (selection), C (subset gloss) major; D, E, F minor; Lemma R.1 added.
