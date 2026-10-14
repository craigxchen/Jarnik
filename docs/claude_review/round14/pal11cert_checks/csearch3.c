/* Exhaustive enumeration of multisets {z_1..z_m} (pool indices i_1 <= ... <= i_m, i_1 a Z_5 orbit representative)
   satisfying the width counts N_v <= tau_A(v) for all test directions v (fatpoint Prop 2.6(c)).
   Pruning (all exact consequences of the counts, no heuristics):
    (a) bitset: a generator touching a saturated direction is discarded;
    (b) knapsack bounds: for each nonnegative integer weight vector Y_q on directions, any completion by r more
        generators satisfies sum w_q(z) <= sum_j Y_q[j] (tau_j - cnt_j), w_q(z) = sum_{j touched by z} Y_q[j];
        every generator can occur at most twice (it touches a coordinate direction e_x, tau = 2), so the left side is
        at least the sum of the r smallest values of w_q over the remaining candidates, each counted twice.
   usage: csearch3 pool.txt weights.txt m out.txt [nparts part]  (top-level orbit reps split round-robin) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define W 32
#define QMAX 64
int NP, ND, m, Q;
int *tau, *rep, **tl, *tn;
uint64_t (*bits)[W];
int cnt[W * 64]; uint64_t full[W];
long long *Y[QMAX], *wq[QMAX], cap[QMAX];
int nparts = 1, part = 0, repcount = 0;
int chosen[64]; long long leaves = 0, nodes[65], pruned_k = 0; FILE *fo;
int *candbuf[65]; int ncand[65];
static int bound_ok(int *c, int n, int r) {
  if (r == 0) return 1;
  for (int q = 0; q < Q; q++) {
    long long best[16]; int nb = 0; /* r smallest values, each candidate counted twice */
    for (int s = 0; s < n; s++) {
      long long v = wq[q][c[s]];
      for (int rep2 = 0; rep2 < 2; rep2++) {
        if (nb < r) { int t = nb++; while (t > 0 && best[t - 1] > v) { best[t] = best[t - 1]; t--; } best[t] = v; }
        else if (v < best[r - 1]) { int t = r - 1; while (t > 0 && best[t - 1] > v) { best[t] = best[t - 1]; t--; } best[t] = v; }
      }
    }
    if (nb < r) return 0;
    long long sum = 0; for (int t = 0; t < r; t++) sum += best[t];
    if (sum > cap[q]) return 0;
  }
  return 1;
}
void rec(int k) {
  nodes[k]++;
  if (k == m) {
    leaves++;
    for (int i = 0; i < m; i++) fprintf(fo, "%d ", chosen[i]);
    fprintf(fo, "\n"); fflush(fo);
    return;
  }
  int *c = candbuf[k], nc = ncand[k];
  for (int t = 0; t < nc; t++) {
    int i = c[t];
    if (k == 0 && !rep[i]) continue;
    if (k == 0) { repcount++; if ((repcount - 1) % nparts != part) continue; }
    uint64_t saved[W]; memcpy(saved, full, sizeof(full));
    for (int j = 0; j < tn[i]; j++) { int d = tl[i][j]; if (++cnt[d] >= tau[d]) full[d >> 6] |= 1ULL << (d & 63); }
    for (int q = 0; q < Q; q++) cap[q] -= wq[q][i];
    chosen[k] = i;
    int okcap = 1; for (int q = 0; q < Q; q++) if (cap[q] < 0) okcap = 0;
    if (okcap) {
      if (k + 1 < m) {
        int *c2 = candbuf[k + 1], n2 = 0;
        for (int s = t; s < nc; s++) {
          int i2 = c[s]; int ok = 1;
          for (int w = 0; w < W; w++) if (bits[i2][w] & full[w]) { ok = 0; break; }
          if (ok) for (int q = 0; q < Q; q++) if (wq[q][i2] > cap[q]) { ok = 0; break; }
          if (ok) c2[n2++] = i2;
        }
        ncand[k + 1] = n2;
        if (n2 > 0 && bound_ok(c2, n2, m - k - 1)) rec(k + 1); else pruned_k++;
      } else rec(k + 1);
    }
    for (int q = 0; q < Q; q++) cap[q] += wq[q][i];
    for (int j = 0; j < tn[i]; j++) cnt[tl[i][j]]--;
    memcpy(full, saved, sizeof(full));
    if (k == 0 && (t % 50 == 0)) { fprintf(stderr, "top %d/%d leaves %lld nodes1 %lld\n", t, nc, leaves, nodes[1]); fflush(stderr); }
    if (k == 0 && rep[i] && ((repcount - 1) % nparts == part)) { fprintf(stderr, "done rep item %d (pool index %d) leaves %lld\n", repcount - 1, i, leaves); fflush(stderr); }
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); FILE *fw = fopen(argv[2], "r"); m = atoi(argv[3]); fo = fopen(argv[4], "w");
  if (argc > 6) { nparts = atoi(argv[5]); part = atoi(argv[6]); }
  if (fscanf(f, "%d %d", &NP, &ND) != 2) return 1;
  if (ND > W * 64) { fprintf(stderr, "too many dirs\n"); return 1; }
  tau = malloc(ND * sizeof(int));
  for (int j = 0; j < ND; j++) if (fscanf(f, "%d", &tau[j]) != 1) return 1;
  rep = malloc(NP * sizeof(int)); tn = malloc(NP * sizeof(int)); tl = malloc(NP * sizeof(int *));
  bits = calloc(NP, sizeof(uint64_t[W]));
  for (int i = 0; i < NP; i++) {
    int zz;
    if (fscanf(f, "%d", &rep[i]) != 1) return 1;
    for (int x = 0; x < 12; x++) if (fscanf(f, "%d", &zz) != 1) return 1;
    if (fscanf(f, "%d", &tn[i]) != 1) return 1;
    tl[i] = malloc((tn[i] + 1) * sizeof(int));
    for (int j = 0; j < tn[i]; j++) { if (fscanf(f, "%d", &tl[i][j]) != 1) return 1; bits[i][tl[i][j] >> 6] |= 1ULL << (tl[i][j] & 63); }
  }
  if (fscanf(fw, "%d", &Q) != 1 || Q > QMAX) return 1;
  for (int q = 0; q < Q; q++) {
    Y[q] = malloc(ND * sizeof(long long)); wq[q] = calloc(NP, sizeof(long long)); cap[q] = 0;
    for (int j = 0; j < ND; j++) { if (fscanf(fw, "%lld", &Y[q][j]) != 1) return 1; cap[q] += Y[q][j] * tau[j]; }
    for (int i = 0; i < NP; i++) for (int j = 0; j < tn[i]; j++) wq[q][i] += Y[q][tl[i][j]];
  }
  for (int k = 0; k <= m; k++) candbuf[k] = malloc(NP * sizeof(int));
  ncand[0] = NP; for (int i = 0; i < NP; i++) candbuf[0][i] = i;
  if (bound_ok(candbuf[0], NP, m)) rec(0); else fprintf(stderr, "root bound fails\n");
  fclose(fo);
  printf("pool %d dirs %d m %d weightvecs %d : leaves %lld ; nodes per level:", NP, ND, m, Q, leaves);
  for (int k = 0; k <= m; k++) printf(" %lld", nodes[k]);
  printf("\n");
  return 0;
}
