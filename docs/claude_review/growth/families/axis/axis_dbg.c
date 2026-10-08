// Search for circles with lattice points at m consecutive abscissae a, a-1, ..., a-m+1:
//   y_j^2 = y_0^2 + 2 j a - j^2  (j = 0..m-1), equivalently n = (a-j)^2 + y_j^2.
// Enumerate y1 <= Y; write 2*y1^2 - 2 = y0^2 + y2^2 (m=3 condition) with y0 < y1 < y2,
// a = (y1^2 - y0^2 + 1)/2, then test y3^2 = 3y1^2 - 2y0^2 - 6, y4^2 = y0^2 + 8a - 16, and the
// backward extension y_{-1}^2 = y0^2 - 2a - 1.  Representations of M = 2(y1-1)(y1+1) as sums of
// two squares are enumerated from the factorisation (smallest-prime-factor sieve).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef unsigned long long u64; typedef __int128 i128; typedef long long ll;
static uint32_t *spf;
static int issq(i128 v, ll *r) {
    if (v < 0) return 0;
    long double s = sqrtl((long double)v);
    ll t = (ll)s;
    for (ll u = t - 2; u <= t + 2; u++) if (u >= 0 && (i128)u * u == v) { if (r) *r = u; return 1; }
    return 0;
}
// gaussian multiply
typedef struct { i128 x, y; } G;
static G gm(G a, G b) { G c = { a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x }; return c; }
static G prim[64]; static int pe[64]; static int np;
static G reps[1<<14]; static int nr;
static void gen(int i, G cur) {
    if (i == np) { reps[nr++] = cur; return; }
    // prime pi with exponent e: choose t copies of pi, e-t of conj
    G p = prim[i], pc = { p.x, -p.y };
    for (int t = 0; t <= pe[i]; t++) {
        G w = cur;
        for (int s = 0; s < t; s++) w = gm(w, p);
        for (int s = 0; s < pe[i]-t; s++) w = gm(w, pc);
        gen(i+1, w);
        if (nr >= (1<<14) - 64) return;
    }
}
static ll two_sq_a, two_sq_b;
static void two_sq(ll p) { // p prime 1 mod 4, find a^2+b^2 = p by brute modular sqrt
    // find x with x^2 = -1 mod p
    ll x = 0;
    for (ll c = 2;; c++) {
        // pow c^((p-1)/4)
        unsigned __int128 r = 1, b = c % p; ll e = (p-1)/4;
        while (e) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; }
        if ((r * r) % p == (unsigned __int128)(p - 1)) { x = (ll)r; break; }
    }
    ll a = p, b = x, lim = (ll)sqrtl((long double)p);
    while (b > lim) { ll t = a % b; a = b; b = t; }
    ll c2 = p - b*b; ll c = (ll)sqrtl((long double)c2);
    while (c*c > c2) c--; while ((c+1)*(c+1) <= c2) c++;
    two_sq_a = b; two_sq_b = c;
}
int main(int argc, char **argv) {
    ll Y = atoll(argv[1]);
    spf = calloc(Y + 3, sizeof(uint32_t));
    for (ll i = 2; i <= Y + 2; i++) if (!spf[i]) for (ll j = i; j <= Y + 2; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    ll cnt3 = 0, cnt4 = 0, cnt5 = 0;
    for (ll y1 = 3; y1 <= Y; y1++) {
        // M = 2 (y1-1)(y1+1); collect prime factorization
        ll f[64]; int fe[64]; int nf = 0;
        ll parts[2] = { y1 - 1, y1 + 1 };
        // merge factorizations
        ll e2 = 1; // factor 2 from the leading 2
        int ok = 1;
        np = 0;
        ll tmp[128]; int nt = 0;
        for (int k = 0; k < 2; k++) { ll v = parts[k]; while (v > 1) { ll p = spf[v]; tmp[nt++] = p; v /= p; } }
        // count exponents
        for (int i = 0; i < nt; i++) {
            int found = 0;
            for (int j = 0; j < nf; j++) if (f[j] == tmp[i]) { fe[j]++; found = 1; break; }
            if (!found) { f[nf] = tmp[i]; fe[nf] = 1; nf++; }
        }
        // include the leading factor 2
        { int found = 0; for (int j = 0; j < nf; j++) if (f[j] == 2) { fe[j]++; found = 1; } if (!found) { f[nf] = 2; fe[nf] = 1; nf++; } }
        G cur = { 1, 0 };
        for (int j = 0; j < nf; j++) {
            ll p = f[j];
            if (p == 2) { G t = { 1, 1 }; for (int s = 0; s < fe[j]; s++) cur = gm(cur, t); }
            else if (p % 4 == 3) { if (fe[j] % 2) { ok = 0; break; } ll pp = 1; for (int s = 0; s < fe[j]/2; s++) pp *= p; G t = { pp, 0 }; cur = gm(cur, t); }
            else { two_sq(p); G t = { two_sq_a, two_sq_b }; prim[np] = t; pe[np] = fe[j]; np++; }
        }
        if (!ok) continue;
        nr = 0; gen(0, cur);
        for (int r = 0; r < nr; r++) {
            for (int u = 0; u < 4; u++) {
                G z = reps[r];
                for (int s = 0; s < u; s++) { G t = { -z.y, z.x }; z = t; }
                ll y0 = (ll)z.x, y2 = (ll)z.y;
                if (y0 < 0 || y2 <= 0) continue;
                if (!(y0 < y1 && y1 < y2)) continue;
                i128 a2 = (i128)y1*y1 - (i128)y0*y0 + 1;
                if (a2 % 2) continue;
                i128 a = a2 / 2;
                if (a < 4) continue;
                cnt3++;
                ll y3, y4, ym;
                int has3 = issq(3*(i128)y1*y1 - 2*(i128)y0*y0 - 6, &y3);
                int hasm = issq((i128)y0*y0 - 2*a - 1, &ym);
                if (getenv("ALL") || has3 || hasm) {
                    cnt4++;
                    int has4 = has3 && issq((i128)y0*y0 + 8*a - 16, &y4);
                    if (has4) cnt5++;
                    i128 n = a*a + (i128)y0*y0;
                    double nf_ = (double)n;
                    // C of the symmetric 8-point cluster (m=4): span 2*y_last over sqrt n * n^(1/4)
                    double C8 = has3 ? 2.0*y3/pow(nf_, 0.25) : 2.0*y2/pow(nf_, 0.25);
                    printf("m4 y0=%lld y1=%lld y2=%lld %s=%lld a=%lld n=%.6e C=%.6f %s\n", y0, y1, y2, has3 ? "y3" : "y-1", has3 ? y3 : ym, (ll)a, nf_, C8, has4 ? "M5!" : "");
                }
            }
        }
    }
    fprintf(stderr, "Y=%lld m3=%lld m4=%lld m5=%lld\n", Y, cnt3, cnt4, cnt5);
    return 0;
}
