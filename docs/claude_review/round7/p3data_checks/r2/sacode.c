/* Simulated annealing for the equal-weight code problem (round 7, p3data, finishing session).
 * Find S in {+-1}^(M x r) (row-major int8, S[x*r+j]) with Hamming distances d_xy >= T_xy for all
 * x < y, all columns nonconstant.  Objective V = sum_{x<y} max(0, T_xy - d_xy).
 * Moves: flip one entry; with probability pf the entry is chosen from a violated pair (WalkSAT style).
 * Returns the best V found; S is overwritten by the best state.  Output is a candidate only: every
 * solution is re-checked in Python (exact integer distances) and then in exact prime arithmetic.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static inline int viol(int t, int d) { return t > d ? t - d : 0; }

long sa(int M, int r, const int *T, int8_t *S, long iters, double T0, double T1, double pf,
        uint64_t seed)
{
    rs = seed * 2654435761ULL + 88172645463325252ULL;
    if (rs == 0) rs = 1;
    int *d = calloc((size_t)M * M, sizeof(int));
    int *cp = calloc(r, sizeof(int));
    int8_t *best = malloc((size_t)M * r);
    long V = 0, bestV;
    for (int j = 0; j < r; j++)
        for (int x = 0; x < M; x++) cp[j] += S[x * r + j] > 0;
    for (int x = 0; x < M; x++)
        for (int y = x + 1; y < M; y++) {
            int c = 0;
            for (int j = 0; j < r; j++) c += S[x * r + j] != S[y * r + j];
            d[x * M + y] = d[y * M + x] = c;
            V += viol(T[x * M + y], c);
        }
    bestV = V;
    memcpy(best, S, (size_t)M * r);
    double lr = log(T1 / T0);
    for (long it = 0; it < iters && bestV > 0; it++) {
        double temp = T0 * exp(lr * (double)it / (double)iters);
        int x = -1, j = -1;
        if (V > 0 && urand() < pf) {
            for (int tries = 0; tries < 64; tries++) {
                int a = rnd() % M, b = rnd() % M;
                if (a == b) continue;
                if (T[a * M + b] > d[a * M + b]) {
                    for (int tj = 0; tj < 64; tj++) {
                        int c = rnd() % r;
                        if (S[a * r + c] == S[b * r + c]) { j = c; break; }
                    }
                    if (j >= 0) x = (rnd() & 1) ? a : b;
                    break;
                }
            }
        }
        if (x < 0) { x = rnd() % M; j = rnd() % r; }
        int sx = S[x * r + j];
        if (sx > 0 && cp[j] == 1) continue;          /* would make column constant (-1) */
        if (sx < 0 && cp[j] == M - 1) continue;      /* would make column constant (+1) */
        long dv = 0;
        for (int y = 0; y < M; y++) {
            if (y == x) continue;
            int dd = d[x * M + y];
            int nd = dd + (S[y * r + j] == sx ? 1 : -1);
            int t = T[x * M + y];
            dv += viol(t, nd) - viol(t, dd);
        }
        if (dv <= 0 || urand() < exp(-(double)dv / temp)) {
            for (int y = 0; y < M; y++) {
                if (y == x) continue;
                int nd = d[x * M + y] + (S[y * r + j] == sx ? 1 : -1);
                d[x * M + y] = d[y * M + x] = nd;
            }
            S[x * r + j] = -sx;
            cp[j] += sx > 0 ? -1 : 1;
            V += dv;
            if (V < bestV) { bestV = V; memcpy(best, S, (size_t)M * r); }
        }
    }
    memcpy(S, best, (size_t)M * r);
    free(d); free(cp); free(best);
    return bestV;
}
