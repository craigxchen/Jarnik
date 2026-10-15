// rfree.c -- referee's forced/free decider.  Input lines: "r K n par_1..par_r" (par_j = required parity of c_j
// mod 2 if n even, -1 = none).  Decides whether some c in (Z/n)^r with the parity has <v,c> != 0 mod n for all
// v with 1<=||v||_1<=K.  Complete DFS: c_1 over divisor representatives (unit orbits), c_2.. exhaustive.
// Output "FREE c1 .. cr" or "FORCED nodes" (or "ABORT").
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int r, K, n, par[16];
static long long nodes, cap = 4000000000LL;
static int c[16];
// S[k][j] = bitset of values <v',c> mod n over v' supported on coords < k with ||v'||_1 <= j
static u64 *S[17][16];
static int W;
static inline void bset(u64 *b, int x) { b[x >> 6] |= 1ULL << (x & 63); }
static inline int bget(const u64 *b, int x) { return (b[x >> 6] >> (x & 63)) & 1; }
static int found = 0, aborted = 0;

static void build_next(int k, int ck) {  // from S[k] build S[k+1] with c_k = ck
    for (int j = 0; j < K; j++) {
        u64 *dst = S[k + 1][j];
        memset(dst, 0, W * 8);
        for (int t = -j; t <= j; t++) {
            int sh = (int)(((long long)t * ck) % n); if (sh < 0) sh += n;
            int jj = j - abs(t);
            u64 *src = S[k][jj];
            for (int x = 0; x < n; x++) if (bget(src, x)) { int y = x + sh; if (y >= n) y -= n; bset(dst, y); }
        }
    }
}
static int ok_value(int k, int x) {  // c_k = x allowed given S[k]?
    for (int t = 1; t <= K; t++) {
        int val = (int)(((long long)t * x) % n);
        int need = (n - val) % n;  // <v',c> = -t x
        if (bget(S[k][K - t], need)) return 0;
    }
    return 1;
}
static void dfs(int k) {
    if (found || aborted) return;
    if (k == r) { found = 1; return; }
    if (++nodes > cap) { aborted = 1; return; }
    int start = 0;
    for (int i = 0; i < n; i++) {
        int x = (start + i) % n;
        if (par[k] >= 0 && (x & 1) != par[k]) continue;
        if (!ok_value(k, x)) continue;
        c[k] = x;
        if (k + 1 < r) build_next(k, x);
        dfs(k + 1);
        if (found || aborted) return;
    }
}
int main(int argc, char **argv) {
    if (argc > 1) cap = atoll(argv[1]);
    while (scanf("%d %d %d", &r, &K, &n) == 3) {
        for (int j = 0; j < r; j++) scanf("%d", &par[j]);
        if (n % 2) for (int j = 0; j < r; j++) par[j] = -1;
        W = (n + 63) / 64;
        for (int k = 0; k <= r; k++) for (int j = 0; j < K; j++) { S[k][j] = calloc(W, 8); }
        for (int j = 0; j < K; j++) bset(S[0][j], 0);  // empty v' gives 0
        nodes = 0; found = 0; aborted = 0;
        // c_0 over divisor reps
        for (int d = 1; d <= n && !found && !aborted; d++) {
            if (n % d) continue;
            int x = d % n;
            if (par[0] >= 0 && (x & 1) != par[0]) continue;
            if (!ok_value(0, x)) continue;
            c[0] = x;
            if (r > 1) build_next(0, x);
            dfs(1);
        }
        if (found) { printf("FREE"); for (int j = 0; j < r; j++) printf(" %d", c[j]); printf("\n"); }
        else if (aborted) printf("ABORT %lld\n", nodes);
        else printf("FORCED %lld\n", nodes);
        fflush(stdout);
        for (int k = 0; k <= r; k++) for (int j = 0; j < K; j++) free(S[k][j]);
    }
    return 0;
}
