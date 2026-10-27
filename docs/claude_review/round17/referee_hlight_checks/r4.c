/* referee: exact min of N = sum_a |(f*chi)(a)| over F=0 4-term f = d_0 + d_1 - d_c - d_e (AGL(1,q) normal form),
   sigma = N - q; also count/verify sigma >= q/2 - 11 - 1.5 sqrt q (checked as 2*sigma+22-q >= 3 sqrt q, squared, integers). */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
  long q = atol(argv[1]);
  signed char *ch = malloc(3*q);
  for (long x = 0; x < q; x++) {
    long v = 1, b = x, e = (q-1)/2; if (x == 0) { ch[x] = 0; continue; }
    while (e) { if (e & 1) v = v*b % q; b = b*b % q; e >>= 1; }
    ch[x] = (v == 1) ? 1 : -1;
  }
  for (long x = 0; x < q; x++) { ch[x+q] = ch[x]; ch[x+2*q] = ch[x]; }
  /* base[a] = chi(a) + chi(a-1) */
  int *base = malloc(sizeof(int)*q);
  for (long a = 0; a < q; a++) base[a] = ch[a] + ch[(a - 1 + q) % q];
  long best = 1L<<60, bc = -1, be = -1, nviol = 0, cnt = 0;
  #pragma omp parallel for schedule(dynamic, 4) reduction(+:nviol,cnt)
  for (long c = 2; c < q; c++) {
    for (long e = c+1; e < q; e++) {
      long N = 0;
      const signed char *pc = ch + q - c, *pe = ch + q - e;   /* chi(a-c) = ch[a + q - c] */
      for (long a = 0; a < q; a++) { int v = base[a] - pc[a] - pe[a]; N += v < 0 ? -v : v; }
      long s = N - q; cnt++;
      long lhs = 2*s + 22 - q;  /* need lhs >= -3 sqrt q */
      if (lhs < 0 && lhs*lhs > 9*q) nviol++;
      #pragma omp critical
      { if (s < best || (s == best && (c < bc || (c == bc && e < be)))) { best = s; bc = c; be = e; } }
    }
  }
  printf("q=%ld 4term_min_sigma=%ld at f=d0+d1-d%ld-d%ld ; (q-3)/2=%ld ; ratio=%.4f ; #f=%ld ; violations_of_Weil_bound=%ld\n",
         q, best, bc, be, (q-3)/2, best/(2.0*q), cnt, nviol);
  return 0;
}
