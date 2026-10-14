/* freescan.c -- exact forced/free status of many moduli for one monomial set V_K (round 14).
 * q is FREE for norm K iff some c in (Z/n)^r (with c_j = par_j mod 2 when n is even) has <v,c> != 0 mod n for
 * every v in V_K.  Exhaustive DFS; column 0 over gcd-orbit representatives (divisors d of n; units of Z/n act
 * on the constraint set {<v,c> != 0} and preserve parity because they are odd when n is even).
 * Input: r nV ; nV monomials (sign-normalised, last nonzero > 0) ; nq ; then per modulus: q n needpar par_0..par_{r-1}
 * Output per modulus: "q n FREE c_0 .. c_{r-1} nodes"  or "q n FORCED nodes"  or "q n ABORT nodes".
 * maxnodes (argv[1], 0 = unlimited).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int r, nV, n; static int (*V)[16]; static int *order, lev_start[17], lev_cnt[17];
static int par[16], needpar; static int c[16]; static int *stamp[16]; static int stampv[16];
static long long nodes, maxnodes; static int found, aborted;
static long long egcd(long long a, long long b, long long *x, long long *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    long long x1, y1; long long g = egcd(b, a % b, &x1, &y1); *x = y1; *y = x1 - (a / b) * y1; return g;
}
static int mark(int j) {
    int st = ++stampv[j]; int *S = stamp[j];
    for (int t = 0; t < lev_cnt[j]; t++) {
        int *v = V[order[lev_start[j] + t]];
        long long p = 0; for (int i = 0; i < j; i++) if (v[i]) p += (long long)v[i] * c[i];
        p %= n; if (p < 0) p += n;
        long long a = v[j], x, y; long long g = egcd(a, n, &x, &y); long long ng = n / g;
        long long inv = ((x % ng) + ng) % ng;
        long long rhs = (n - p) % n; if (rhs % g) continue;
        long long x0 = ((rhs / g) % ng) * inv % ng;
        for (long long s = 0; s < g; s++) S[x0 + s * ng] = st;
    }
    return st;
}
static void dfs(int j) {
    if (found || aborted) return;
    nodes++; if (maxnodes && nodes > maxnodes) { aborted = 1; return; }
    if (j == r) { found = 1; return; }
    int st = mark(j); int *S = stamp[j];
    if (j == 0) {
        for (int d = 1; d <= n; d++) { if (n % d) continue; int x = d % n;
            if (S[x] == st) continue; if (needpar && (x & 1) != par[0]) continue;
            c[0] = x; dfs(1); if (found || aborted) return; }
        return;
    }
    for (int x = 0; x < n; x++) {
        if (S[x] == st) continue; if (needpar && (x & 1) != par[j]) continue;
        c[j] = x; dfs(j + 1); if (found || aborted) return;
    }
}
int main(int argc, char **argv) {
    maxnodes = argc > 1 ? atoll(argv[1]) : 0;
    if (scanf("%d %d", &r, &nV) != 2) return 1;
    V = malloc(sizeof(*V) * nV); order = malloc(sizeof(int) * nV);
    memset(lev_cnt, 0, sizeof lev_cnt);
    for (int k = 0; k < nV; k++) { int last = -1; for (int j = 0; j < r; j++) { if (scanf("%d", &V[k][j]) != 1) return 1; if (V[k][j]) last = j; }
        if (V[k][last] <= 0) return 2; lev_cnt[last]++; }
    lev_start[0] = 0; for (int j = 1; j <= r; j++) lev_start[j] = lev_start[j - 1] + lev_cnt[j - 1];
    { int fill[17] = {0}; for (int k = 0; k < nV; k++) { int last = -1; for (int j = 0; j < r; j++) if (V[k][j]) last = j; order[lev_start[last] + fill[last]++] = k; } }
    int nq; if (scanf("%d", &nq) != 1) return 1;
    int cap = 0;
    for (int t = 0; t < nq; t++) {
        int q; if (scanf("%d %d %d", &q, &n, &needpar) != 3) return 1;
        for (int j = 0; j < r; j++) if (scanf("%d", &par[j]) != 1) return 1;
        if (n > cap) { for (int j = 0; j < r; j++) { free(stamp[j]); stamp[j] = calloc(n + 1, sizeof(int)); stampv[j] = 0; } cap = n; }
        nodes = 0; found = 0; aborted = 0; dfs(0);
        if (found) { printf("%d %d FREE", q, n); for (int j = 0; j < r; j++) printf(" %d", c[j]); printf(" %lld\n", nodes); }
        else printf("%d %d %s %lld\n", q, n, aborted ? "ABORT" : "FORCED", nodes);
        fflush(stdout);
    }
    return 0;
}
