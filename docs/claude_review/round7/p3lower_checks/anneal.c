/* Simulated annealing for prescribed-correlation sign matrices (evidence only).
 *
 * Input (stdin):  M r K iters seed restarts
 *                 then M*M integers k_xy (symmetric, zero diagonal)
 *   K >= 0 : number of fixed Sylvester Hadamard blocks (M must be a power of 2 if K > 0);
 *            they contribute exactly -K to every off-diagonal correlation.
 *   r      : number of FREE +-1 columns searched.
 * Goal: <s_x, s_y>_free - K <= -k_xy for all x != y.
 * Energy: sum_{x<y} h(G_xy - K + k_xy), h(v) = max(0, v) (+ v^2/8 to smooth).
 * Output: best violation count, best energy, and (if feasible) "FEASIBLE" and the free block,
 *         re-verified from scratch with integer arithmetic.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static unsigned long long rs;
static inline unsigned long long xr(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double ur(void){ return (xr() >> 11) * (1.0/9007199254740992.0); }

static inline double h(int v){ return v > 0 ? (double)v + 0.125*(double)v*(double)v : 0.0; }

int main(void){
  int M, r, K, restarts; long long iters; unsigned long long seed;
  if (scanf("%d %d %d %lld %llu %d", &M, &r, &K, &iters, &seed, &restarts) != 6) return 1;
  int *k = malloc(sizeof(int)*M*M);
  for (int i = 0; i < M*M; i++) if (scanf("%d", &k[i]) != 1) return 1;
  signed char *S = malloc(M*r), *best = malloc(M*r);
  int *G = malloc(sizeof(int)*M*M);
  rs = seed * 2654435761ULL + 12345ULL; if (!rs) rs = 1;
  int bestviol = 1<<30; double bestE = 1e300;
  for (int rep = 0; rep < restarts; rep++){
    for (int i = 0; i < M*r; i++) S[i] = (xr() & 1) ? 1 : -1;
    for (int x = 0; x < M; x++) for (int y = 0; y < M; y++){
      int g = 0; for (int j = 0; j < r; j++) g += S[x*r+j]*S[y*r+j]; G[x*M+y] = g; }
    double E = 0; for (int x = 0; x < M; x++) for (int y = x+1; y < M; y++) E += h(G[x*M+y] - K + k[x*M+y]);
    double T0 = 2.0, T1 = 0.02;
    for (long long it = 0; it < iters; it++){
      double T = T0 * pow(T1/T0, (double)it/(double)iters);
      int x = xr() % M, j = xr() % r;
      int sx = S[x*r+j];
      double dE = 0;
      for (int z = 0; z < M; z++){ if (z == x) continue;
        int g = G[x*M+z], d = -2*sx*S[z*r+j];
        dE += h(g + d - K + k[x*M+z]) - h(g - K + k[x*M+z]); }
      if (dE <= 0 || ur() < exp(-dE/T)){
        S[x*r+j] = -sx;
        for (int z = 0; z < M; z++){ if (z == x) continue; int d = -2*sx*S[z*r+j]; G[x*M+z] += d; G[z*M+x] += d; }
        E += dE;
        if (E < 0.5){ break; }
      }
    }
    int viol = 0; for (int x = 0; x < M; x++) for (int y = x+1; y < M; y++) if (G[x*M+y] - K + k[x*M+y] > 0) viol++;
    if (viol < bestviol || (viol == bestviol && E < bestE)){ bestviol = viol; bestE = E; memcpy(best, S, M*r); }
    if (bestviol == 0) break;
  }
  /* exact re-verification */
  int viol = 0;
  for (int x = 0; x < M; x++) for (int y = x+1; y < M; y++){
    int g = 0; for (int j = 0; j < r; j++) g += best[x*r+j]*best[y*r+j];
    if (g - K + k[x*M+y] > 0) viol++; }
  printf("%d %.3f\n", viol, bestE);
  if (viol == 0){ printf("FEASIBLE\n");
    for (int x = 0; x < M; x++){ for (int j = 0; j < r; j++) putchar(best[x*r+j] > 0 ? '+' : '-'); putchar('\n'); } }
  return 0;
}
