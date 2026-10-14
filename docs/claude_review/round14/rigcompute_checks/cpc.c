/* cpc.c -- C version of the exact joint DFS of cp_joint.py / cp_mix.py (round 14).
 * Input (stdin):
 *   npairs ; B[0..npairs-1] (doubles, budgets log cap/2) ; full (eps bitmask of all classes, 0 = no coherence test)
 *   nmod ; per modulus:  lq isC nopt ; per option: emask len e1 m1 e2 m2 ...   (pair index, multiplicity)
 * Semantics: choose one option per modulus; feasible iff for every pair e: sum_q lq * mult_q(e) <= B[e] + TOL.
 * If full != 0: additionally require AND over isC moduli of emask == 0 (non-coherence).
 * Propagation (sound): mandatory load of a modulus = lq * min over its viable options of mult(e) (0 if some viable
 * option lacks e); an option is dropped if some element overflows on top of the other moduli's mandatory load;
 * if full != 0, prune when some eps bit is present in every viable option of every isC modulus.
 * Output: "INFEASIBLE nodes" | "FOUND nodes k_1 ... k_nmod" | "ABORTED nodes".  argv[1] = nodecap, argv[2] = maxsol
 * (number of FOUND solutions to print before stopping; the driver verifies each).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define TOL 1e-9
typedef struct { int len; int *e; int *mu; int emask; } Opt;
static int npairs, nmod, full; static double *B; static double *lq; static int *isC; static int *nopt; static Opt **opt;
static long long nodes = 0, nodecap; static int aborted = 0, nsol = 0, maxsol;
static int *cnt; static int *mn; static double *load; static double **mand; /* per modulus mandatory vectors (sparse via touched list) */
static int **touched; static int *ntouched;

static int propagate(char **via, int *nv) {
    for (;;) {
        for (int e = 0; e < npairs; e++) load[e] = 0;
        for (int q = 0; q < nmod; q++) {
            /* mandatory of q */
            int nt = 0;
            for (int k = 0; k < nopt[q]; k++) if (via[q][k]) { Opt *o = &opt[q][k];
                for (int t = 0; t < o->len; t++) { int e = o->e[t];
                    if (cnt[e] == 0) { touched[q][nt++] = e; mn[e] = o->mu[t]; } else if (o->mu[t] < mn[e]) mn[e] = o->mu[t];
                    cnt[e]++; } }
            ntouched[q] = 0;
            for (int t = 0; t < nt; t++) { int e = touched[q][t];
                if (cnt[e] == nv[q]) { touched[q][ntouched[q]] = e; mand[q][ntouched[q]] = lq[q] * mn[e]; ntouched[q]++; load[e] += lq[q] * mn[e]; }
                cnt[e] = 0; }
        }
        for (int e = 0; e < npairs; e++) if (load[e] > B[e] + TOL) return 0;
        int changed = 0;
        for (int q = 0; q < nmod; q++) {
            if (nv[q] <= 1) continue;
            /* other = load - mandatory_q */
            for (int t = 0; t < ntouched[q]; t++) load[touched[q][t]] -= mand[q][t];
            for (int k = 0; k < nopt[q]; k++) if (via[q][k]) { Opt *o = &opt[q][k]; int ok = 1;
                for (int t = 0; t < o->len; t++) if (load[o->e[t]] + lq[q] * o->mu[t] > B[o->e[t]] + TOL) { ok = 0; break; }
                if (!ok) { via[q][k] = 0; nv[q]--; changed = 1; } }
            for (int t = 0; t < ntouched[q]; t++) load[touched[q][t]] += mand[q][t];
            if (nv[q] == 0) return 0;
        }
        if (full) {
            int killable = 0;
            for (int q = 0; q < nmod; q++) if (isC[q]) for (int k = 0; k < nopt[q]; k++) if (via[q][k]) killable |= (~opt[q][k].emask) & full;
            if (killable != full) return 0;
        }
        if (!changed) return 1;
    }
}
static void dfs(char **via, int *nv) {
    if (aborted || nsol >= maxsol) return;
    nodes++; if (nodes > nodecap) { aborted = 1; return; }
    if (!propagate(via, nv)) return;
    int best = -1, bestC = -1;
    for (int q = 0; q < nmod; q++) if (nv[q] > 1) {
        if (isC[q] && (bestC < 0 || nv[q] < nv[bestC])) bestC = q;
        if (best < 0 || nv[q] < nv[best]) best = q; }
    if (full && bestC >= 0) best = bestC;
    if (best < 0) {
        if (full) { int a = full; for (int q = 0; q < nmod; q++) if (isC[q]) for (int k = 0; k < nopt[q]; k++) if (via[q][k]) a &= opt[q][k].emask; if (a) return; }
        printf("FOUND %lld", nodes);
        for (int q = 0; q < nmod; q++) for (int k = 0; k < nopt[q]; k++) if (via[q][k]) printf(" %d", k);
        printf("\n"); fflush(stdout); nsol++; return;
    }
    /* copy state */
    char **v2 = malloc(sizeof(char *) * nmod); int *n2 = malloc(sizeof(int) * nmod);
    for (int q = 0; q < nmod; q++) { v2[q] = malloc(nopt[q]); }
    for (int k = 0; k < nopt[best]; k++) if (via[best][k]) {
        for (int q = 0; q < nmod; q++) { memcpy(v2[q], via[q], nopt[q]); n2[q] = nv[q]; }
        memset(v2[best], 0, nopt[best]); v2[best][k] = 1; n2[best] = 1;
        dfs(v2, n2);
        if (aborted || nsol >= maxsol) break;
    }
    for (int q = 0; q < nmod; q++) free(v2[q]); free(v2); free(n2);
}
int main(int argc, char **argv) {
    nodecap = argc > 1 ? atoll(argv[1]) : 1000000; maxsol = argc > 2 ? atoi(argv[2]) : 1;
    if (scanf("%d", &npairs) != 1) return 1;
    B = malloc(sizeof(double) * npairs); for (int e = 0; e < npairs; e++) if (scanf("%lf", &B[e]) != 1) return 1;
    if (scanf("%d %d", &full, &nmod) != 2) return 1;
    lq = malloc(sizeof(double) * nmod); isC = malloc(sizeof(int) * nmod); nopt = malloc(sizeof(int) * nmod); opt = malloc(sizeof(Opt *) * nmod);
    int maxlen = 0;
    for (int q = 0; q < nmod; q++) {
        if (scanf("%lf %d %d", &lq[q], &isC[q], &nopt[q]) != 3) return 1;
        opt[q] = malloc(sizeof(Opt) * nopt[q]);
        for (int k = 0; k < nopt[q]; k++) { Opt *o = &opt[q][k];
            if (scanf("%d %d", &o->emask, &o->len) != 2) return 1;
            o->e = malloc(sizeof(int) * (o->len + 1)); o->mu = malloc(sizeof(int) * (o->len + 1));
            for (int t = 0; t < o->len; t++) if (scanf("%d %d", &o->e[t], &o->mu[t]) != 2) return 1; }
    }
    cnt = calloc(npairs, sizeof(int)); mn = calloc(npairs, sizeof(int)); load = calloc(npairs, sizeof(double));
    mand = malloc(sizeof(double *) * nmod); touched = malloc(sizeof(int *) * nmod); ntouched = calloc(nmod, sizeof(int));
    for (int q = 0; q < nmod; q++) { mand[q] = malloc(sizeof(double) * (npairs + 1)); touched[q] = malloc(sizeof(int) * (npairs + 1)); }
    char **via = malloc(sizeof(char *) * nmod); int *nv = malloc(sizeof(int) * nmod);
    for (int q = 0; q < nmod; q++) { via[q] = malloc(nopt[q]); memset(via[q], 1, nopt[q]); nv[q] = nopt[q]; }
    dfs(via, nv);
    if (nsol == 0) printf("%s %lld\n", aborted ? "ABORTED" : "INFEASIBLE", nodes);
    else printf("DONE %lld %d %s\n", nodes, nsol, aborted ? "aborted" : "");
    return 0;
}
