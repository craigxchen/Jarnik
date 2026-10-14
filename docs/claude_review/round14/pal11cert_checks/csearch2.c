/* Exhaustive enumeration of multisets {z_1..z_m} (pool indices i_1 <= ... <= i_m) such that
   #{t : <z_{i_t}', v> != 0} <= tau_A(v) for every test direction v (width counts, necessary for a product
   certificate: fatpoint Prop 2.6(c)).  Symmetry: i_1 must be a Z_5 orbit representative (complete up to G).
   Bitset implementation: 'full' = set of directions whose count has reached tau; a generator is admissible iff its
   touched set misses 'full'.  Candidate lists are filtered level by level (exact; no heuristic pruning).
   usage: csearch2 pool.txt m out.txt [maxleaves_to_write] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define W 32
int NP, ND, m;
int *tau, *rep, **tl, *tn;
uint64_t (*bits)[W];
int cnt[W * 64]; uint64_t full[W];
int chosen[64]; long long leaves = 0, nodes[65]; FILE *fo; long long maxw = 100000000LL;
int *candbuf[65]; int ncand[65];
void rec(int k) {
  nodes[k]++;
  if (k == m) {
    leaves++;
    if (leaves <= maxw) { for (int i = 0; i < m; i++) fprintf(fo, "%d ", chosen[i]); fprintf(fo, "\n"); }
    return;
  }
  int *c = candbuf[k], nc = ncand[k];
  for (int t = 0; t < nc; t++) {
    int i = c[t];
    if (k == 0 && !rep[i]) continue;
    /* add i */
    uint64_t saved[W]; memcpy(saved, full, sizeof(full));
    for (int j = 0; j < tn[i]; j++) { int d = tl[i][j]; if (++cnt[d] >= tau[d]) full[d >> 6] |= 1ULL << (d & 63); }
    chosen[k] = i;
    if (k + 1 < m) {
      int *c2 = candbuf[k + 1], n2 = 0;
      for (int s = t; s < nc; s++) {
        int i2 = c[s]; int ok = 1;
        for (int w = 0; w < W; w++) if (bits[i2][w] & full[w]) { ok = 0; break; }
        if (ok) c2[n2++] = i2;
      }
      ncand[k + 1] = n2;
      if (n2 > 0) rec(k + 1);
    } else rec(k + 1);
    for (int j = 0; j < tn[i]; j++) cnt[tl[i][j]]--;
    memcpy(full, saved, sizeof(full));
    if (k == 0) { fprintf(stderr, "top %d/%d leaves %lld\n", t, nc, leaves); fflush(stderr); }
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); m = atoi(argv[2]); fo = fopen(argv[3], "w");
  if (argc > 4) maxw = atoll(argv[4]);
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
  for (int k = 0; k <= m; k++) candbuf[k] = malloc(NP * sizeof(int));
  ncand[0] = NP; for (int i = 0; i < NP; i++) candbuf[0][i] = i;
  rec(0);
  fclose(fo);
  printf("pool %d dirs %d m %d : leaves %lld ; nodes per level:", NP, ND, m, leaves);
  for (int k = 0; k <= m; k++) printf(" %lld", nodes[k]);
  printf("\n");
  return 0;
}
