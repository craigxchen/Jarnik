/* pass-2 referee: exhaustive 4-term F=0 characters f = d_0 + d_1 - d_c - d_e (AGL(1,q) normal form, 2<=c<e<=q-1).
   sigma = sum_a max(|f(a)|,|(f*chi)(a)|) - q (formula re-derived by referee and checked against the pattern definition).
   Also S_T = sum_a chi(a(a-1)(a-c)(a-e)). Reports: min sigma, #normal forms attaining it, an argmin, max |S_T|,
   #violations of 2 sigma + 22 - q >= -3 sqrt q (Lemma 2.7, integer test), and of Hasse |S_T| <= 2 sqrt q + 1. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
int main(int argc, char **argv) {
  int q = atoi(argv[1]);
  signed char *ch2 = malloc(2*q);
  for (int u = 0; u < q; u++) ch2[u] = -1;
  for (long x = 1; x < q; x++) ch2[(x*x) % q] = 1;
  ch2[0] = 0;
  for (int u = 0; u < q; u++) ch2[u+q] = ch2[u];      /* ch2[v] = chi(v mod q) for 0<=v<2q */
  int *A = malloc(sizeof(int)*q);
  for (int a = 0; a < q; a++) A[a] = ch2[a] + ch2[(a - 1 + q)];
  long gmin = 1L<<40, cnt = 0, viol = 0, hviol = 0; int bc = -1, be = -1; long maxST = 0;
  #pragma omp parallel
  {
    long lmin = 1L<<40, lcnt = 0, lviol = 0, lhviol = 0, lmaxST = 0; int lbc = -1, lbe = -1;
    #pragma omp for schedule(dynamic, 4)
    for (int c = 2; c < q; c++) {
      for (int e = c + 1; e < q; e++) {
        long N = 0, ST = 0;
        for (int a = 0; a < q; a++) {
          int x = ch2[a - c + q], y = ch2[a - e + q];
          int cv = A[a] - x - y;
          int ac = cv < 0 ? -cv : cv;
          int fa = (a == 0 || a == 1 || a == c || a == e) ? 1 : 0;
          N += (ac > fa ? ac : fa);
          ST += ch2[a] * ch2[a - 1 + q] * x * y;
        }
        long sig = N - q;
        long ast = ST < 0 ? -ST : ST;
        if (ast > lmaxST) lmaxST = ast;
        long t = 2*sig + 22 - q;                 /* sigma >= q/2-11-1.5sqrt q  <=>  t >= -3 sqrt q  <=>  t >= 0 or t^2 <= 9q */
        if (t < 0 && t*t > 9L*q) lviol++;
        long h = ast - 1;                        /* Hasse: |S_T| <= 2 sqrt q + 1  <=>  h <= 0 or h^2 <= 4q */
        if (h > 0 && h*h > 4L*q) lhviol++;
        if (sig < lmin) { lmin = sig; lcnt = 1; lbc = c; lbe = e; }
        else if (sig == lmin) lcnt++;
      }
    }
    #pragma omp critical
    {
      if (lmin < gmin) { gmin = lmin; cnt = lcnt; bc = lbc; be = lbe; }
      else if (lmin == gmin) cnt += lcnt;
      viol += lviol; hviol += lhviol; if (lmaxST > maxST) maxST = lmaxST;
    }
  }
  printf("q=%d min_sigma=%ld count_nf=%ld argmin=(c=%d,e=%d) sigma/r=%.5f (q-3)/2=%.1f weilLB=%.3f viol_L27=%ld maxST=%ld 2sqrtq+1=%.2f 3sqrtq=%.2f hasse_viol=%ld\n",
         q, gmin, cnt, bc, be, gmin/(2.0*q), (q-3)/2.0, q/2.0 - 11 - 1.5*sqrt(q), viol, maxST, 2*sqrt(q)+1, 3*sqrt(q), hviol);
  return 0;
}
