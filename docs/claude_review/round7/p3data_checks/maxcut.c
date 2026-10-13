/* Exact weighted max-cut by Gray-code enumeration (vertex 0 fixed on side 0).
 * best_cuts(M, nu[M*M] symmetric, K, masks_out[K], vals_out[K]):
 *   enumerates all 2^(M-1) subsets T of {1..M-1} (side 1), cut value = sum_{x in T, y notin T} nu[x][y],
 *   returns the K largest cut values (with their masks, bit x = side of vertex x) in descending order.
 *   The empty set (constant cut) is skipped.
 * Used only for LP pricing; every value it returns is recomputed in Python before use.
 * Build: gcc -O3 -march=native -shared -fPIC -o libmaxcut.so maxcut.c
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void heap_push(double *hv, uint64_t *hm, int *n, int K, double v, uint64_t m) {
    /* min-heap of size <= K keeping the K largest values */
    if (*n < K) {
        int i = (*n)++;
        hv[i] = v; hm[i] = m;
        while (i > 0) {
            int p = (i - 1) / 2;
            if (hv[p] <= hv[i]) break;
            double tv = hv[p]; hv[p] = hv[i]; hv[i] = tv;
            uint64_t tm = hm[p]; hm[p] = hm[i]; hm[i] = tm;
            i = p;
        }
    } else if (v > hv[0]) {
        hv[0] = v; hm[0] = m;
        int i = 0;
        for (;;) {
            int l = 2 * i + 1, r = l + 1, s = i;
            if (l < K && hv[l] < hv[s]) s = l;
            if (r < K && hv[r] < hv[s]) s = r;
            if (s == i) break;
            double tv = hv[s]; hv[s] = hv[i]; hv[i] = tv;
            uint64_t tm = hm[s]; hm[s] = hm[i]; hm[i] = tm;
            i = s;
        }
    }
}

int best_cuts(int M, const double *nu, int K, uint64_t *masks_out, double *vals_out) {
    double *g = (double *)malloc(sizeof(double) * M);   /* g[v] = sum_y nu[v][y] * (side(y) differs ? -1 : +1) gain of flipping v */
    double *hv = (double *)malloc(sizeof(double) * K);
    uint64_t *hm = (uint64_t *)malloc(sizeof(uint64_t) * K);
    int n = 0;
    /* start: all on side 0, cut = 0; flipping v from side s to 1-s changes cut by
       sum_{y same side} nu[v][y] - sum_{y other side} nu[v][y] */
    int side[64];
    memset(side, 0, sizeof(side));
    for (int v = 0; v < M; v++) {
        double s = 0;
        for (int y = 0; y < M; y++) if (y != v) s += nu[v * M + y];
        g[v] = s;
    }
    double cut = 0;
    uint64_t mask = 0;
    uint64_t total = (uint64_t)1 << (M - 1);
    for (uint64_t i = 1; i < total; i++) {
        int b = __builtin_ctzll(i);
        int v = b + 1;                    /* vertices 1..M-1 */
        cut += g[v];
        side[v] ^= 1;
        mask ^= ((uint64_t)1 << v);
        /* update gains: for y != v, nu[v][y] moves between same/other side */
        const double *row = nu + (size_t)v * M;
        for (int y = 0; y < M; y++) {
            if (y == v) continue;
            if (side[y] == side[v]) g[y] += 2 * row[y];   /* v joined y's side: y's flip gain increases */
            else g[y] -= 2 * row[y];
        }
        g[v] = -g[v];
        heap_push(hv, hm, &n, K, cut, mask);
    }
    /* sort descending */
    for (int a = 0; a < n; a++)
        for (int c = a + 1; c < n; c++)
            if (hv[c] > hv[a]) {
                double tv = hv[a]; hv[a] = hv[c]; hv[c] = tv;
                uint64_t tm = hm[a]; hm[a] = hm[c]; hm[c] = tm;
            }
    for (int a = 0; a < n; a++) { masks_out[a] = hm[a]; vals_out[a] = hv[a]; }
    free(g); free(hv); free(hm);
    return n;
}
