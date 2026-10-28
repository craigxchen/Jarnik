# Judge notes, round 19
- Plan J1: for all 140 search records with gain > -40 (D=4,5,6; all gain -32), compute independently
  k(f) (anchor rule of referee m1), N(u) by subset enumeration, n_1 = #{T: nu_T(B) > mu_T} exactly,
  and test Thm 3.3 bound gain <= Lambda_k - N - n_1 (expect equality, k=4).
- J1 done: 140 records, all k=4, Thm3.3 ok; equality 106 ((32,16):34, (48,0):72 special); (32,0): 34 (28 lam=-1, 6 lam=1).
- J1b done: n'=48 on all 140, Cor 2.5 equality on all 140.
- (dropped) Plan J2: Lambda_k.beta = 32 - 2^{8-k} recomputation from K formula (tiny).
