/* Distinct level-1 signatures at one modulus (round 14).  Enumerates parity-respecting l in (Z/m)^r, l_1 over
   U_1-orbit representatives (signatures are U_1-invariant: u = 1 mod 4 preserves units), skipping vectors with a
   collision (v,s) that is individually infeasible (feas[v*4+s] == 0).  Signature = bitset over (v,s).
   Each distinct signature gets flags: bit0 = realised by a genuine-class vector, bit1 = by a non-class vector. */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXR 8
static int R, M, NN, NV, NW;
static const int *VV; static const unsigned char *FEAS; static const int64_t *CLS; static int NCLS; static int NE; static const int *CLSOFF;
static const int *NUV; static int **byidx; static int *nby; static int l[MAXR];
static uint64_t *tab; static unsigned char *tflag; static int64_t tcap, tcnt; static int overflow;
static uint64_t cur[64];
static int incls(void) {   /* returns bitmask of eps-classes containing the current vector */
    int64_t code = 0;
    for (int j = 0; j < R; j++) code = code * NN + (l[j] % NN);
    int mask = 0;
    for (int k = 0; k < NE; k++) {
        int lo = CLSOFF[k], hi = CLSOFF[k + 1] - 1;
        while (lo <= hi) { int mid = (lo + hi) / 2; if (CLS[mid] == code) { mask |= 1 << k; break; } if (CLS[mid] < code) lo = mid + 1; else hi = mid - 1; }
    }
    return mask;
}
static void insert(int flag) {
    cur[NW - 1] = (uint64_t)flag;
    uint64_t h = 1469598103934665603ULL;
    for (int w = 0; w < NW; w++) { h ^= cur[w]; h *= 1099511628211ULL; h ^= h >> 29; }
    int64_t pos = (int64_t)(h % (uint64_t)tcap);
    while (1) {
        if (tflag[pos] == 0) {
            if (tcnt * 2 >= tcap) { overflow = 1; return; }
            memcpy(tab + pos * NW, cur, NW * 8); tflag[pos] = 1; tcnt++; return;
        }
        if (memcmp(tab + pos * NW, cur, NW * 8) == 0) { return; }
        pos++; if (pos == tcap) pos = 0;
    }
}
static void leaf(void) {
    memset(cur, 0, NW * 8);
    for (int i = 0; i < NV; i++) {
        long t = 0; const int *v = VV + (size_t)i * R;
        for (int j = 0; j < R; j++) t += (long)v[j] * l[j];
        t %= M; if (t < 0) t += M;
        if (t % NN == 0) { int b = i * 4 + (int)(t / NN); cur[b >> 6] |= 1ULL << (b & 63); }
    }
    insert(incls());
}
static void dfs(int j) {
    if (overflow) return;
    if (j == R) { leaf(); return; }
    for (int x = NUV[j]; x < M; x += 2) {
        l[j] = x; int ok = 1;
        for (int k = 0; k < nby[j] && ok; k++) {
            int i = byidx[j][k]; const int *v = VV + (size_t)i * R;
            long t = 0; for (int jj = 0; jj <= j && jj < MAXR; jj++) t += (long)v[jj] * l[jj];
            t %= M; if (t < 0) t += M;
            if (t % NN == 0 && !FEAS[i * 4 + t / NN]) ok = 0;
        }
        if (ok) dfs(j + 1);
    }
}
/* returns number of distinct signatures (or -1 on overflow); signatures in outtab (cap*nw words), flags in outflag */
int64_t run(int r, int m, int nvec, const int *V, const unsigned char *feas, const int *lastidx, const int64_t *cls, int ne, const int *clsoff,
            const int *nuv, int nreps, const int *reps, int64_t cap, uint64_t *outtab, unsigned char *outflag) {
    R = r; M = m; NN = m / 4; NV = nvec; VV = V; FEAS = feas; CLS = cls; NE = ne; CLSOFF = clsoff; NUV = nuv;
    NW = (nvec * 4 + 63) / 64 + 1; tab = outtab; tflag = outflag; tcap = cap; tcnt = 0; overflow = 0;
    memset(outflag, 0, cap);
    byidx = malloc(sizeof(int*) * r); nby = calloc(r, sizeof(int));
    for (int j = 0; j < r; j++) byidx[j] = malloc(sizeof(int) * (nvec + 1));
    for (int i = 0; i < nvec; i++) { int j = lastidx[i]; byidx[j][nby[j]++] = i; }
    for (int a = 0; a < nreps && !overflow; a++) {
        l[0] = reps[a]; int ok = 1;
        for (int k = 0; k < nby[0] && ok; k++) {
            int i = byidx[0][k]; long t = ((long)V[(size_t)i * r] * l[0]) % m; if (t < 0) t += m;
            if (t % NN == 0 && !feas[i * 4 + t / NN]) ok = 0;
        }
        if (ok) dfs(1);
    }
    for (int j = 0; j < r; j++) free(byidx[j]);
    free(byidx); free(nby);
    return overflow ? -1 : tcnt;
}
