// Faster integral-point search for depth sets D={0,d1,...,dm}: n=a^2+y0^2=(a-d)^2+y_d^2.
// (y_{d1}-y0)(y_{d1}+y0) = 2 d1 a - d1^2 ; u=y_{d1}-y0, v=y_{d1}+y0 (same parity), y0<=K sqrt(a).
// Quadratic-residue filters mod 64, 63, 65, 11 before an exact integer square root.
// usage: axis3 A K d1 d2 ... ; progress to stderr.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef long long ll; typedef unsigned long long u64; typedef __int128 i128;
static uint8_t q64[64], q63[63], q65[65], q11[11];
static int issq(u64 v, ll *r) {
    if (!q64[v & 63] || !q63[v % 63] || !q65[v % 65] || !q11[v % 11]) return 0;
    u64 t = (u64)sqrtl((long double)v);
    while (t * t > v) t--;
    while ((t + 1) * (t + 1) <= v) t++;
    if (t * t == v) { *r = (ll)t; return 1; }
    return 0;
}
int main(int argc, char **argv) {
    double A = atof(argv[1]); double K = atof(argv[2]);
    int m = argc - 3; ll D[16];
    for (int i = 0; i < m; i++) D[i] = atoll(argv[3 + i]);
    for (int i = 0; i < 64; i++) q64[(i * i) % 64] = 1;
    for (int i = 0; i < 63; i++) q63[(i * i) % 63] = 1;
    for (int i = 0; i < 65; i++) q65[(i * i) % 65] = 1;
    for (int i = 0; i < 11; i++) q11[(i * i) % 11] = 1;
    ll d1 = D[0];
    ll U = (ll)sqrt(2.0 * d1 * A) + 2;
    ll full = 0; double tested = 0;
    for (ll u = 1; u <= U; u++) {
        if (u % 50000 == 0) { fprintf(stderr, "progress u=%lld complete for a <= %.3e tested=%.3e\n", u, (double)u * u / (2.0 * d1), tested); fflush(stderr); }
        // y0 <= K sqrt(a):  (v-u)^2/4 <= K^2 (u v + d1^2)/(2 d1)  ->  quadratic in v
        // v^2 - (2u + c u) v + u^2 - c d1^2 <= 0 with c = 2K^2/d1
        double c = 2.0 * K * K / d1;
        double B = 2.0 * u + c * u, Cc = (double)u * u - c * d1 * d1;
        double disc = B * B - 4 * Cc;
        double vmaxK = (B + sqrt(disc > 0 ? disc : 0)) / 2.0 + 1;
        double vmaxA = (2.0 * d1 * A - (double)d1 * d1) / (double)u;
        double vmax = vmaxK < vmaxA ? vmaxK : vmaxA;
        if (vmax < u) break;
        for (ll v = u; v <= (ll)vmax; v += 2) {
            i128 num = (i128)u * v + d1 * d1;
            if (num % (2 * d1)) continue;
            ll a = (ll)(num / (2 * d1));
            if (a <= D[m - 1]) continue;
            ll y0 = (v - u) / 2;
            tested++;
            ll ys[16]; int ok = 1;
            for (int j = 1; j < m; j++) {
                ll d = D[j];
                i128 val = (i128)y0 * y0 + (i128)2 * d * a - (i128)d * d;
                if (val < 0 || !issq((u64)val, &ys[j])) { ok = 0; break; }
            }
            if (ok) {
                full++;
                double n = (double)a * a + (double)y0 * y0;
                printf("FULL a=%lld y0=%lld y%lld=%lld", a, y0, d1, (u + v) / 2);
                for (int t = 1; t < m; t++) printf(" y%lld=%lld", D[t], ys[t]);
                printf(" n=%.6e C~%.6f\n", n, 2.0 * ys[m - 1] / pow(n, 0.25));
                fflush(stdout);
            }
        }
    }
    fprintf(stderr, "DONE A=%.3g K=%.2f tested=%.3e full=%lld\n", A, K, tested, full);
    return 0;
}
