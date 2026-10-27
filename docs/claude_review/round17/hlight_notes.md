# hlight notes (working)

## Formula
rows inf + F_q. Patterns P_a = a - D, P'_a = a + D, D = {0} cup QR (sdc_sym.py: P_a={a} cup {x: chi(a-x)=1}).
lambda = (f on F_q, lambda_inf = -F), F = sum f.
lambda(P_a) = (h(a)+c(a))/2, lambda(P'_a) = (h(a)-c(a))/2, h = F + f, c = f*chi, (f*chi)(a)=sum_t f(a-t)chi(t).
N(lambda) := sum_S |lambda(S)| = sum_a max(|h(a)|,|c(a)|);  sigma = N - q.
chi*chi = 1 - q delta_0 (q=3 mod 4). c*chi = F - q f.
Sum c^2 = q||f||^2 - F^2. Sum h = (q+1)F.
Parity c(a) = h(a) mod 2.
Reductions: |F|>=2 => N >= (q+1)|F| => sigma>=q+2.
 q|f(b)| <= |F| + ||c||_1 - |c(b)|  => ||f||_inf>=2, |F|<=1 => sigma >= q-1.
So only ternary f, F in {0,+-1} matter.
Dual pair: f = (g_{b1}-g_{b2})/2 + corr, g_b(x)=chi(x-b): c spiky (+-q/2), N ~ 1.5 q, sigma ~ r/4 as 4-term.
Hadamard view (exact.md Thm 4.1 iii): sigma = 1 + D/(2n) for ternary, D = sum |u_a|(n-|u_a|)+|v_a|(n-|v_a|).

## KEY PROOF (found)
Defect identity: N = sum_S|t_S|, t_S = lambda(S), n=||lambda||_1, |t_S|<=n/2,
 sum_S t_S^2 = (q+1)||lambda||_2^2/2  (from formula: sum h^2 = ||f||^2+(q+2)F^2, sum c^2 = q||f||^2-F^2).
 ternary: sigma = 1 + (2/n) sum_S |t_S|(n/2-|t_S|). nontrivial entry (t notin {0,+-n/2}) gives >= (1-2/n).
Polynomial zero bound: c(a) = sum f(t) chi(a-t) == P(a) := sum f(t)(a-t)^m mod q, m=(q-1)/2.
 P == 0 iff moments M_j = sum f(t) t^j = 0 for j<=m; Vandermonde => impossible if 1<=|supp f|<=m+1 and f not == 0 mod q.
 deg P <= m - j0, j0 = first nonvanishing moment. F=0 => j0>=1 => #zeros(c) <= (q-3)/2. Pairs attain equality.
F=0 ternary, n>=4 even, n<=(q+1)/2: off-T values even; N0<=(q-3)/2; N_n <= (B-4M)/n^2, B = qn - n;
 M >= ((q+3)/2 - n - (q-1)/n)/(1-4/n^2); sigma >= n-1 + (2(n-2)/n) M  => g(n) = (q(n-2)+4n-n^2)/(n+2).
 g(4)=q/3, g(6)=q/2-3/2, g(8)=0.6q-3.2. Also sigma >= n-1 always (F=0 ternary).
|F|=1 ternary: n>=6: sigma >= 1+(1-2/n)(q-n/2+1) >= (q+1)/2; n=4: sigma >= (q-1)/2.
Non-ternary f or |F|>=2: sigma >= q-1.
Dual pairs sigma = (q-1)/2 exactly (computed q<=3251) => c_0 <= (q-1)/(4q).
Weil for n=4 (F=0): sigma >= q/2 - O(sqrt q).
