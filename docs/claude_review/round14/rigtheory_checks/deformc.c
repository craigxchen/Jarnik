/* One-modulus deformation enumerator (round 14).
   Enumerate parity-respecting l in (Z/m)^r with l_1 restricted to given orbit representatives (weights = orbit sizes
   under U_1 = {u in (Z/m)^*, u = 1 mod 4}); a vector survives iff every v (one per +-pair) colliding at level 1,
   <v,l> = s*n mod m, has allowed[v][s] = 1.  Survivors are classified:
     inclass  : c = l mod n lies in the genuine level-1 class {u L mod n : u unit mod n}   (sorted codes given)
     strict   : collision signature {(v,s)} equals the genuine one
     F2       : outside class, every colliding v is a genuine level-1 collision (<v,L> = 0 mod n)
     F1       : outside class, some colliding v is not a genuine collision.
   Counts are weighted (exact totals).  Also stores up to MAXEX example survivors outside the class. */
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXR 8
static int R, M, NN, NV;
static const int *VV;      /* NV x R */
static const unsigned char *ALW;  /* NV x 4 */
static const int *LASTIDX; /* NV */
static const unsigned char *GCOLL; /* NV: genuine level-1 collision (any unit) */
static const int *GSIG;    /* NV: genuine unit s if collides at level 1 else -1 */
static const int64_t *CLS; static int NCLS;
static const int *NUV;
static int **byidx; static int *nby;
static int l[MAXR];
static int64_t cnt[6]; /* 0 total,1 inclass,2 strict,3 F2,4 F1 */
static int64_t W;
static int *exbuf; static int nex, maxex;
static int incls(void) {
    int64_t code = 0;
    for (int j = 0; j < R; j++) code = code * NN + (l[j] % NN);
    int lo = 0, hi = NCLS - 1;
    while (lo <= hi) { int mid = (lo + hi) / 2; if (CLS[mid] == code) return 1; if (CLS[mid] < code) lo = mid + 1; else hi = mid - 1; }
    return 0;
}
static void leaf(void) {
    int strict = 1, f1 = 0;
    for (int i = 0; i < NV; i++) {
        long t = 0; const int *v = VV + (size_t)i * R;
        for (int j = 0; j < R; j++) t += (long)v[j] * l[j];
        t %= M; if (t < 0) t += M;
        int s = -1;
        if (t % NN == 0) s = (int)(t / NN);
        if (s != GSIG[i]) strict = 0;
        if (s >= 0 && !GCOLL[i]) f1 = 1;
    }
    cnt[0] += W;
    if (incls()) { cnt[1] += W; if (strict) cnt[2] += W; }
    else {
        if (f1) cnt[4] += W; else cnt[3] += W;
        if (nex < maxex) { for (int j = 0; j < R; j++) exbuf[nex * R + j] = l[j]; nex++; }
    }
}
static void dfs(int j) {
    if (j == R) { leaf(); return; }
    int start = NUV[j];
    for (int x = start; x < M; x += 2) {
        l[j] = x;
        int ok = 1;
        for (int k = 0; k < nby[j] && ok; k++) {
            int i = byidx[j][k]; const int *v = VV + (size_t)i * R;
            long t = 0; for (int jj = 0; jj <= j; jj++) t += (long)v[jj] * l[jj];
            t %= M; if (t < 0) t += M;
            if (t % NN == 0 && !ALW[i * 4 + t / NN]) ok = 0;
        }
        if (ok) dfs(j + 1);
    }
}
/* reps/weights for l_1; returns counts in out[0..4]; examples in ex (maxex*R) and number in *nexout */
void run(int r, int m, int nvec, const int *V, const unsigned char *alw, const int *lastidx,
         const unsigned char *gcoll, const int *gsig, const int64_t *cls, int ncls, const int *nuv,
         int nreps, const int *reps, const int64_t *wts, int64_t *out, int maxex_, int *ex, int *nexout) {
    R = r; M = m; NN = m / 4; NV = nvec; VV = V; ALW = alw; LASTIDX = lastidx; GCOLL = gcoll; GSIG = gsig;
    CLS = cls; NCLS = ncls; NUV = nuv; maxex = maxex_; exbuf = ex; nex = 0;
    memset(cnt, 0, sizeof cnt);
    byidx = malloc(sizeof(int*) * r); nby = calloc(r, sizeof(int));
    for (int j = 0; j < r; j++) byidx[j] = malloc(sizeof(int) * (nvec + 1));
    for (int i = 0; i < nvec; i++) { int j = lastidx[i]; byidx[j][nby[j]++] = i; }
    for (int a = 0; a < nreps; a++) {
        l[0] = reps[a]; W = wts[a];
        int ok = 1;
        for (int k = 0; k < nby[0] && ok; k++) {
            int i = byidx[0][k]; long t = ((long)V[(size_t)i * r] * l[0]) % m; if (t < 0) t += m;
            if (t % NN == 0 && !alw[i * 4 + t / NN]) ok = 0;
        }
        if (ok) dfs(1);
    }
    for (int k = 0; k < 5; k++) out[k] = cnt[k];
    *nexout = nex;
    for (int j = 0; j < r; j++) free(byidx[j]);
    free(byidx); free(nby);
}
