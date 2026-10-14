/* Exhaustive enumeration of multisets {z_1..z_m} from a generator pool that satisfy the width counts
   #{i : <z_i',v> != 0} <= tau_A(v) for every test direction v (necessary for prod <z_i',u> in J(A)_m^perp,
   fatpoint Prop 2.6(c)).  Symmetry: only multisets whose minimum element is a Z_5 orbit representative
   (complete up to the group action; see pal11cert.md §2).
   usage: csearch pool.txt m out.txt */
#include <stdio.h>
#include <stdlib.h>
int NP, ND, m;
int *tau, *rep, **tl, *tn, (*z)[12];
int *cnt; int chosen[64]; long long leaves = 0, nodes = 0; FILE *fo;
int inc[12];
void rec(int k, int start) {
  nodes++;
  if (k == m) {
    leaves++;
    for (int i = 0; i < m; i++) fprintf(fo, "%d ", chosen[i]);
    fprintf(fo, "\n");
    return;
  }
  /* capacity: remaining incidence must allow m-k more generators of support >= 2 */
  int cap = 0; for (int x = 0; x < 12; x++) cap += (x == 0 ? 4 : 2) - inc[x];
  if (cap < 2 * (m - k)) return;
  for (int i = start; i < NP; i++) {
    if (k == 0 && !rep[i]) continue;
    if (k == 0) { fprintf(stderr, "top %d nodes %lld leaves %lld\n", i, nodes, leaves); fflush(stderr); }
    int ok = 1;
    for (int j = 0; j < tn[i]; j++) if (cnt[tl[i][j]] >= tau[tl[i][j]]) { ok = 0; break; }
    if (!ok) continue;
    for (int j = 0; j < tn[i]; j++) cnt[tl[i][j]]++;
    for (int x = 0; x < 12; x++) if (z[i][x]) inc[x]++;
    chosen[k] = i;
    rec(k + 1, i);
    for (int j = 0; j < tn[i]; j++) cnt[tl[i][j]]--;
    for (int x = 0; x < 12; x++) if (z[i][x]) inc[x]--;
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); m = atoi(argv[2]); fo = fopen(argv[3], "w");
  if (fscanf(f, "%d %d", &NP, &ND) != 2) return 1;
  tau = malloc(ND * sizeof(int)); cnt = calloc(ND, sizeof(int));
  for (int j = 0; j < ND; j++) if (fscanf(f, "%d", &tau[j]) != 1) return 1;
  rep = malloc(NP * sizeof(int)); tn = malloc(NP * sizeof(int)); tl = malloc(NP * sizeof(int *));
  z = malloc(NP * sizeof(int[12]));
  for (int i = 0; i < NP; i++) {
    if (fscanf(f, "%d", &rep[i]) != 1) return 1;
    for (int x = 0; x < 12; x++) if (fscanf(f, "%d", &z[i][x]) != 1) return 1;
    if (fscanf(f, "%d", &tn[i]) != 1) return 1;
    tl[i] = malloc((tn[i] + 1) * sizeof(int));
    for (int j = 0; j < tn[i]; j++) if (fscanf(f, "%d", &tl[i][j]) != 1) return 1;
  }
  rec(0, 0);
  fclose(fo);
  printf("pool %d dirs %d m %d : nodes %lld leaves (multisets passing all counts, min element an orbit rep) %lld\n",
         NP, ND, m, nodes, leaves);
  return 0;
}
