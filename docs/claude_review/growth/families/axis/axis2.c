// Integral points for a depth set D = {0, d1, d2, ...}: n = a^2 + y0^2 = (a-d)^2 + y_d^2.
// Parametrise (y_{d1} - y0)(y_{d1} + y0) = 2 d1 a - d1^2 by u = y_{d1}-y0, v = y_{d1}+y0,
// with y0 <= K sqrt(a) (bounded normalized arc constant).  Report every a <= A whose
// tuple satisfies all depths, and the count of partial matches (all but the last depth).
// usage: axis2 A K d1 d2 ... dm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll; typedef __int128 i128;
static int issq(i128 v, ll *r) {
    if (v < 0) return 0;
    long double s = sqrtl((long double)v); ll t = (ll)s;
    for (ll u = t - 2; u <= t + 2; u++) if (u >= 0 && (i128)u*u == v) { if (r) *r = u; return 1; }
    return 0;
}
int main(int argc, char **argv) {
    double A = atof(argv[1]); double K = atof(argv[2]);
    int m = argc - 3; ll D[16];
    for (int i = 0; i < m; i++) D[i] = atoll(argv[3 + i]);
    ll d1 = D[0];
    ll U = (ll)sqrt(2.0 * d1 * A) + 2;
    ll partial = 0, full = 0, tested = 0;
    for (ll u = 1; u <= U; u++) {
        if (u % 20000 == 0) { fprintf(stderr, "progress u=%lld (complete for a <= %.3e)\n", u, (double)u*u/(2.0*d1)); }
        // v >= u, same parity, (u v + d1^2) divisible by 2 d1, a <= A, y0=(v-u)/2 <= K sqrt(a)
        for (ll v = u; ; v += 2) {
            i128 num = (i128)u * v + d1 * d1;
            double a_f = (double)num / (2.0 * d1);
            if (a_f > A) break;
            ll y0 = (v - u) / 2;
            if ((double)y0 > K * sqrt(a_f) + 1) { if (v > 4 * u + 100) break; else continue; }
            if (num % (2 * d1)) continue;
            ll a = (ll)(num / (2 * d1));
            if (a <= D[m - 1]) continue;
            tested++;
            int ok = 1; int j;
            ll ys[16];
            for (j = 1; j < m; j++) {
                ll d = D[j];
                if (!issq((i128)y0 * y0 + (i128)2 * d * a - (i128)d * d, &ys[j])) { ok = 0; break; }
            }
            if (j >= m - 1) partial++;
            if (ok) {
                full++;
                double n = (double)a * a + (double)y0 * y0;
                double ymax = (double)ys[m - 1];
                printf("FULL a=%lld y0=%lld y_d1=%lld", a, y0, (u + v) / 2);
                for (int t = 1; t < m; t++) printf(" y%lld=%lld", D[t], ys[t]);
                printf(" n=%.6e C=%.6f\n", n, 2.0 * ymax / pow(n, 0.25));
                fflush(stdout);
            }
        }
    }
    fprintf(stderr, "A=%.3g K=%.2f tested=%lld partial(all but last)=%lld full=%lld\n", A, K, tested, partial, full);
    return 0;
}
