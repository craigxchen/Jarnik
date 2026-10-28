# notes for width.md
- (W) at M=5 is exactly Conjecture C_bal of threshold/construct.md §2.8 (pi_bal=3/5, balanced cuts |P| in {2,3}).
- M=4 analogue: reg <= (1/4) sum_{|P|=2} w_P  = (1/2) sum over 3 cuts 2|2; data t_max(m)=3m for k=4.
- slicing by x_1 with Thm-B style induction gives 15m/4 on mR_5 if slice bounds are mean-width (fails).
## KEY RESULT (found): recursive peeling proves t_max(m) <= floor(10m/3) for all m
- T_n(a,b) := {y in Z^n : sum y = a-b, sum y_+ <= a} = lattice pts of a*Delta - b*Delta. S_m = T_5(m,m).
- slice y_1=c: c>=0 -> T_{n-1}(a-c,b); c<0 -> T_{n-1}(a,b-|c|). T_2(a,b): a+b+1 collinear pts, reg a+b.
- Peeling lemma: reg(S) <= max_i(reg(S_(i)) + i - 1) any order (proof: P = p + (phi-c_1) q induction).
- F_n = a+b+H_n; H_2=0, H_3=min(a,b), H_4=min(floor((a+b)/2), 2min), H_5=min(min+floor(max/3), 3min).
- sufficient recursion (R''): for 0<=p<=a, 0<=q<=b: min(H_{n-1}(a-p,b)+q, H_{n-1}(a,b-q)+p) <= H_n(a,b), F_{n-1} monotone.
- n=5 proof: a<=b: (i) Y<=2a+p<=3a; (ii) (1/3)X+(2/3)Y with X<=2(a-p)+q, Y<=(a+b-q)/2+p gives a+b/3.
- n=4: average gives (a+b)/2; X<=2b, Y<=2a.
- DP check: closed forms = DP exactly on [0,90]^2 (rho_verify.py/out). rho_6(m,m)=floor(10m/3) too (m<=30); rho_7 bigger.
- (W) holds for all T_5(a,b): F_5 <= 5(a+b)/3 = (1/6) sum widths.
- width.md §1, §2 written (Thm 2.5 with Lemmas 2.1-2.4). Next: §3 (W) for T_5(a,b) family + general status; §4 comps;
  §5 consequences; §0 summary; §6 audit.
- pairbox_peel.py: 60 random pair-box polytopes: iterated coordinate peeling <= W always (max ratio 1.0), reg==iterpeel 31/40.
- treg checks: T_4(a,b) a<=b<=7 reg=F_4 (full at F, defic at F-1); T_5 a<=b<=4 full at F all; (2,3) full already at 7.
- DONE: width.md complete (§0-§6). gp_test 150 certified (W) instances, 0 violations, 5 equalities.
