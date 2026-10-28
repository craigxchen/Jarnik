/* Independent rank check for toric.md Thm 3.4(b)(c).
   S = { x in Z^k : sum x = 0, ||x||_1 <= 2m }, coordinates x_1..x_{k-1}.
   Matrix: rows = points of S, columns = monomials of total degree <= d in k-1 variables,
   entries mod p = 2^31-1.  Prints rank.  Full rank mod p  => full rank over Q
   => no nonzero Laurent polynomial supported on S vanishes at e to order > d (Lemma 3.1).
   Rank < |S| mod p does NOT certify deficiency over Q (only consistent with it).
   usage: regrank k m d */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define P 2147483647ULL
static int k, m, d, nv;
static int npts = 0, cap = 0; static int *pts = NULL;
static void addpt(int *x) {
    if (npts == cap) { cap = cap ? 2*cap : 1024; pts = realloc(pts, sizeof(int)*cap*nv); }
    memcpy(pts + (size_t)npts*nv, x, sizeof(int)*nv); npts++;
}
static void genpts(int i, int *x, int sum, int l1) {
    if (i == nv) { int last = -sum; if (l1 + abs(last) <= 2*m) addpt(x); return; }
    for (int v = -2*m; v <= 2*m; v++) { if (l1 + abs(v) > 2*m) continue; x[i] = v; genpts(i+1, x, sum+v, l1+abs(v)); }
}
static int nmon = 0, mcap = 0; static int *mons = NULL;
static void addmon(int *e) {
    if (nmon == mcap) { mcap = mcap ? 2*mcap : 1024; mons = realloc(mons, sizeof(int)*mcap*nv); }
    memcpy(mons + (size_t)nmon*nv, e, sizeof(int)*nv); nmon++;
}
static void genmon(int i, int *e, int deg) {
    if (i == nv) { addmon(e); return; }
    for (int t = 0; t + deg <= d; t++) { e[i] = t; genmon(i+1, e, deg+t); }
}
static uint64_t pw(int64_t b, int e) {
    uint64_t r = 1, bb = (uint64_t)(((b % (int64_t)P) + (int64_t)P) % (int64_t)P);
    while (e) { if (e & 1) r = r*bb % P; bb = bb*bb % P; e >>= 1; }
    return r;
}
static uint64_t inv(uint64_t a) { return pw((int64_t)a, (int)(P-2)); }
int main(int argc, char **argv) {
    k = atoi(argv[1]); m = atoi(argv[2]); d = atoi(argv[3]); nv = k-1;
    int x[16], e[16];
    genpts(0, x, 0, 0); genmon(0, e, 0);
    int R = npts, C = nmon;
    uint32_t *A = malloc(sizeof(uint32_t)*(size_t)R*C);
    if (!A) { fprintf(stderr, "alloc fail\n"); return 1; }
    for (int r = 0; r < R; r++) {
        int *pp = pts + (size_t)r*nv;
        /* power table */
        uint64_t pt[16][64];
        for (int i = 0; i < nv; i++) for (int t = 0; t <= d; t++) pt[i][t] = pw(pp[i], t);
        for (int c = 0; c < C; c++) {
            int *ee = mons + (size_t)c*nv; uint64_t v = 1;
            for (int i = 0; i < nv; i++) v = v*pt[i][ee[i]] % P;
            A[(size_t)r*C + c] = (uint32_t)v;
        }
    }
    int rank = 0;
    for (int c = 0; c < C && rank < R; c++) {
        int piv = -1;
        for (int r = rank; r < R; r++) if (A[(size_t)r*C + c]) { piv = r; break; }
        if (piv < 0) continue;
        if (piv != rank) { for (int j = c; j < C; j++) { uint32_t t = A[(size_t)piv*C+j]; A[(size_t)piv*C+j] = A[(size_t)rank*C+j]; A[(size_t)rank*C+j] = t; } }
        uint32_t *pr = A + (size_t)rank*C;
        uint64_t iv = inv(pr[c]);
        for (int j = c; j < C; j++) pr[j] = (uint32_t)((uint64_t)pr[j]*iv % P);
        for (int r = rank+1; r < R; r++) {
            uint32_t *rr = A + (size_t)r*C;
            uint64_t f = rr[c]; if (!f) continue;
            uint64_t nf = P - f;
            for (int j = c; j < C; j++) rr[j] = (uint32_t)((rr[j] + nf*pr[j]) % P);
        }
        rank++;
    }
    printf("k=%d m=%d d=%d |S|=%d monomials=%d rank_mod_p=%d %s\n", k, m, d, R, C, rank, rank == R ? "FULL" : "deficient");
    return 0;
}
