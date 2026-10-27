/* deepsig.c -- (derived from csp.c) multi-level signatures at a small modulus q: labels mod M = m q^(a-1);
 * forbidden conditions <v,ell> = f mod D (D | M); at leaves every monomial's collision level b (largest b<=a with
 * n q^(b-1) | <v,ell>) is computed and leaves with b > bmax[v][t] are discarded (t = level-1 residue index).
 * csp.c -- exact exhaustive DFS for single-modulus label problems (round 14, rigcompute).
 *
 * Variables ell_0..ell_{r-1} in Z/M.  Constraint per monomial v (sign-normalised: last nonzero coeff > 0):
 *     <v, ell> mod M  not in Forb(v)   (a list of <= 8 residues mod M).
 * Parity: par[j] in {-1,0,1}; if >=0, ell_j = par[j] mod 2 (M even required).
 * Column 0 ranges over a supplied list of (value, weight) pairs (orbit representatives under a symmetry group
 * that preserves the constraint system; weight = orbit size).  The python driver is responsible for the
 * validity of the symmetry reduction.
 * Classification: c = ell mod n; leaf is "genuine class" iff c is in the supplied hash set.
 * Modes: 0 = count all leaves (weighted) and non-genuine leaves; also record up to maxwit non-genuine witnesses
 *        1 = stop after maxwit non-genuine witnesses (exhaustive if fewer exist)
 *        2 = stop at the first leaf (existence test; used for free vectors)
 * maxnodes: abort (reported) if exceeded (0 = unlimited).
 * Output lines:  "W e0 e1 ... "  witnesses;  "RESULT complete|aborted nodes leaves_w nongen_w leaves_raw"
 * All arithmetic exact (64-bit ints; M < 2^31; counts in unsigned __int128).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
static int M, n, r, mode, maxwit; static long long maxnodes;
static int par[16];
static int ncol0; static int *col0v; static long long *col0w;
static int nV; static int (*V)[16]; static int *nf; static int (*forb)[8]; static int (*forbD)[8]; static int (*bmax)[4]; static int qq, nn1;
static int *lev_start, *lev_cnt; static int *order;  /* monomials grouped by last index */
static int nclass; static uint64_t *hkeys; static int *hvals; static int hsize; static int (*cls)[16]; static int *cmask;
static int ell[16];
static int *stamp[16]; static int stampv[16];
static long long nodes = 0; static int aborted = 0;
static u128 leaves_w = 0, nongen_w = 0; static unsigned long long leaves_raw = 0; static int nwit = 0;
static long long curw = 1;

static void print128(u128 x) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (x == 0) { printf("0"); return; }
    while (x) { buf[--i] = '0' + (int)(x % 10); x /= 10; }
    printf("%s", buf + i);
}
static uint64_t hashvec(const int *c) {
    uint64_t h = 1469598103934665603ULL;
    for (int j = 0; j < r; j++) { h ^= (uint64_t)(c[j] + 1); h *= 1099511628211ULL; h ^= h >> 29; }
    return h | 1ULL;
}
static void hinsert(const int *c, int idx) {
    uint64_t h = hashvec(c); int p = (int)(h % (uint64_t)hsize);
    while (hkeys[p]) { p = (p + 1) % hsize; }
    hkeys[p] = h; hvals[p] = idx;
}
static int hfind(const int *c);
static int hfind(const int *c) {
    uint64_t h = hashvec(c); int p = (int)(h % (uint64_t)hsize);
    while (hkeys[p]) {
        if (hkeys[p] == h) { int k = hvals[p]; int ok = 1; for (int j = 0; j < r; j++) if (cls[k][j] != c[j]) { ok = 0; break; } if (ok) return cmask[k]; }
        p = (p + 1) % hsize;
    }
    return 0;
}
static long long egcd(long long a, long long b, long long *x, long long *y) {
    if (b == 0) { *x = 1; *y = 0; return a; }
    long long x1, y1; long long g = egcd(b, a % b, &x1, &y1); *x = y1; *y = x1 - (a / b) * y1; return g;
}
static int mark(int j) {  /* mark forbidden values for column j given ell_0..ell_{j-1}; returns stamp */
    int st = ++stampv[j]; int *S = stamp[j];
    for (int t = 0; t < lev_cnt[j]; t++) {
        int k = order[lev_start[j] + t]; int *v = V[k];
        long long p = 0;
        for (int i = 0; i < j; i++) if (v[i]) p += (long long)v[i] * ell[i];
        p %= M; if (p < 0) p += M;
        long long a = v[j];  /* > 0 */
        for (int f = 0; f < nf[k]; f++) {
            long long D = forbD[k][f];
            long long x, y; long long g = egcd(a, D, &x, &y); long long Dg = D / g;
            long long inv = ((x % Dg) + Dg) % Dg;
            long long rhs = ((long long)forb[k][f] - p) % D; if (rhs < 0) rhs += D;
            if (rhs % g) continue;
            long long x0 = ((rhs / g) % Dg) * inv % Dg;
            for (long long t = x0; t < M; t += Dg) S[t] = st;
        }
    }
    return st;
}

/* ---- mode 3: distinct collision signatures of all leaves ---- */
typedef struct { uint64_t h; int len; int *el; int orm; int none; int notL; u128 cnt; int rep[16]; int repN[16]; int nmk; int mks[32]; int mkrep[32][16]; } Sig;
static Sig *sigs; static int nsig = 0, capsig = 0; static int *sht; static int shsize = 0;
static int *coll_buf;
static void sig_init(void) { capsig = 1024; sigs = malloc(sizeof(Sig) * capsig); shsize = 1 << 23; sht = malloc(sizeof(int) * shsize); memset(sht, -1, sizeof(int) * shsize); coll_buf = malloc(sizeof(int) * (4 * nV + 4)); }
static void sig_leaf(int mk) {
    int len = 0;
    for (int k = 0; k < nV; k++) {
        long long d = 0; for (int j = 0; j < r; j++) d += (long long)V[k][j] * ell[j];
        d %= M; if (d < 0) d += M;
        if (d % nn1 == 0) {
            int t = (int)((d % (4LL * nn1)) / nn1); int b = 1; long long dd = (long long)nn1 * qq;
            while (dd <= M / 4 && d % dd == 0) { b++; dd *= qq; }
            if (b > bmax[k][t]) return;      /* inadmissible depth: leaf discarded */
            coll_buf[len++] = 64 * k + 16 * t + b; }
    }
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)len;
    for (int t = 0; t < len; t++) { h ^= (uint64_t)(coll_buf[t] + 7); h *= 1099511628211ULL; h ^= h >> 31; }
    int p = (int)(h & (shsize - 1));
    while (sht[p] >= 0) {
        Sig *s = &sigs[sht[p]];
        if (s->h == h && s->len == len && !memcmp(s->el, coll_buf, sizeof(int) * len)) {
            s->orm |= mk; { int f = 0; for (int t = 0; t < s->nmk; t++) if (s->mks[t] == mk) f = 1; if (!f && s->nmk < 32) { s->mks[s->nmk] = mk; memcpy(s->mkrep[s->nmk], ell, sizeof(int) * 16); s->nmk++; } } if (mk == 0) s->none = 1; if (!(mk & 1) && !s->notL) { s->notL = 1; memcpy(s->repN, ell, sizeof(int) * 16); } s->cnt += (u128)curw; return; }
        p = (p + 1) & (shsize - 1);
    }
    if (nsig == capsig) { capsig *= 2; sigs = realloc(sigs, sizeof(Sig) * capsig); }
    if (nsig * 2 > shsize) { fprintf(stderr, "signature table full\n"); exit(3); }
    Sig *s = &sigs[nsig]; s->h = h; s->len = len; s->el = malloc(sizeof(int) * (len + 1)); memcpy(s->el, coll_buf, sizeof(int) * len);
    s->orm = mk; s->none = (mk == 0); s->notL = !(mk & 1); s->cnt = (u128)curw; memcpy(s->rep, ell, sizeof(int) * 16); memcpy(s->repN, ell, sizeof(int) * 16); s->nmk = 1; s->mks[0] = mk; memcpy(s->mkrep[0], ell, sizeof(int) * 16); sht[p] = nsig++;
}
static void sig_dump(void) {
    for (int i = 0; i < nsig; i++) { Sig *s = &sigs[i];
        printf("S %d %d %d ", s->orm, s->none, s->notL); print128(s->cnt); printf(" %d", s->len);
        for (int j = 0; j < r; j++) printf(" %d", s->rep[j]); for (int j = 0; j < r; j++) printf(" %d", s->repN[j]);
        printf(" %d", s->nmk); for (int t = 0; t < s->nmk; t++) { printf(" %d", s->mks[t]); for (int j = 0; j < r; j++) printf(" %d", s->mkrep[t][j]); }
        for (int t = 0; t < s->len; t++) printf(" %d", s->el[t]); printf("\n"); }
}
static void leaf(void) {
    leaves_raw++; leaves_w += (u128)curw;
    int c[16]; for (int j = 0; j < r; j++) c[j] = ell[j] % n;
    int mk = hfind(c);
    if (mode == 3) { sig_leaf(mk); return; }
    if (!(mk & 1)) {
        nongen_w += (u128)curw;
        if (nwit < maxwit) { printf("W"); for (int j = 0; j < r; j++) printf(" %d", ell[j]); printf("\n"); nwit++; }
    }
}
static int stop(void) {
    if (aborted) return 1;
    if (mode == 2 && leaves_raw > 0) return 1;
    if (mode == 1 && nwit >= maxwit) return 1;
    return 0;
}
static void dfs(int j) {
    if (stop()) return;
    nodes++; if (maxnodes && nodes > maxnodes) { aborted = 1; return; }
    if (j == r) { leaf(); return; }
    if (j == 0) {
        int st = mark(0); int *S = stamp[0];
        for (int t = 0; t < ncol0; t++) {
            int x = col0v[t]; if (S[x] == st) continue;
            ell[0] = x; curw = col0w[t]; dfs(1); if (stop()) return;
        }
        return;
    }
    int st = mark(j); int *S = stamp[j];
    int step = (par[j] >= 0) ? 2 : 1; int x0 = (par[j] >= 0) ? par[j] : 0;
    for (int x = x0; x < M; x += step) {
        if (S[x] == st) continue;
        ell[j] = x; dfs(j + 1); if (stop()) return;
    }
}
int main(void) {
    if (scanf("%d %d %d %d %d %lld %d %d", &M, &n, &r, &mode, &maxwit, &maxnodes, &qq, &nn1) != 8) return 1;
    for (int j = 0; j < r; j++) scanf("%d", &par[j]);
    scanf("%d", &ncol0); col0v = malloc(sizeof(int) * ncol0); col0w = malloc(sizeof(long long) * ncol0);
    for (int t = 0; t < ncol0; t++) scanf("%d %lld", &col0v[t], &col0w[t]);
    scanf("%d", &nV); V = malloc(sizeof(*V) * nV); nf = malloc(sizeof(int) * nV); forb = malloc(sizeof(*forb) * nV); forbD = malloc(sizeof(*forbD) * nV); bmax = malloc(sizeof(*bmax) * nV);
    for (int k = 0; k < nV; k++) {
        for (int j = 0; j < r; j++) scanf("%d", &V[k][j]);
        scanf("%d", &nf[k]); for (int f = 0; f < nf[k]; f++) scanf("%d %d", &forb[k][f], &forbD[k][f]);
        for (int t = 0; t < 4; t++) scanf("%d", &bmax[k][t]);
    }
    lev_start = calloc(r + 1, sizeof(int)); lev_cnt = calloc(r + 1, sizeof(int)); order = malloc(sizeof(int) * nV);
    for (int k = 0; k < nV; k++) { int last = -1; for (int j = 0; j < r; j++) if (V[k][j]) last = j; lev_cnt[last]++; }
    for (int j = 1; j <= r; j++) lev_start[j] = lev_start[j - 1] + lev_cnt[j - 1];
    { int *fill = calloc(r, sizeof(int));
      for (int k = 0; k < nV; k++) { int last = -1; for (int j = 0; j < r; j++) if (V[k][j]) last = j;
        if (V[k][last] <= 0) { fprintf(stderr, "monomial not normalised\n"); return 2; }
        order[lev_start[last] + fill[last]++] = k; } free(fill); }
    scanf("%d", &nclass); cls = malloc(sizeof(*cls) * (nclass + 1)); cmask = malloc(sizeof(int) * (nclass + 1));
    hsize = 4 * nclass + 7; hkeys = calloc(hsize, sizeof(uint64_t)); hvals = calloc(hsize, sizeof(int));
    for (int k = 0; k < nclass; k++) { for (int j = 0; j < r; j++) scanf("%d", &cls[k][j]); scanf("%d", &cmask[k]); hinsert(cls[k], k); }
    if (mode == 3) sig_init();
    for (int j = 0; j < r; j++) { stamp[j] = calloc(M, sizeof(int)); stampv[j] = 0; }
    dfs(0);
    if (mode == 3) sig_dump();
    printf("RESULT %s %lld ", aborted ? "aborted" : "complete", nodes);
    print128(leaves_w); printf(" "); print128(nongen_w); printf(" %llu\n", leaves_raw);
    return 0;
}
