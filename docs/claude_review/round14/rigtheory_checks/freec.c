/* Exact free/forced test (round 14).  q is free for norm K iff some c in (Z/n)^r with
   (n even => c_j = nu_j mod 2) has <v,c> != 0 mod n for every v with 1<=||v||_1<=K.
   DFS over columns; c_1 restricted to divisors of n (orbits of (Z/n)^* on Z/n; the unit action
   preserves parity of c_j when n is even since units are odd).  Reachable sets P[k] = {<w,c>: ||w||_1<=k, supp w in prefix}. */
#include <stdlib.h>
#include <string.h>
static int N, R, KK, EVEN; static const int *NU;
static int *sol;
static int dfs(int j, unsigned char *P) {
    if (j == R) return 1;
    unsigned char *Q = malloc((size_t)(KK + 1) * N);
    for (int c = (j == 0 ? 1 : 0); c < N; c++) {
        if (j == 0 && N % c != 0) continue;
        if (EVEN && (c % 2) != NU[j]) continue;
        int ok = 1;
        for (int t = 1; t <= KK && ok; t++) if (P[(size_t)(KK - t) * N + ((long)t * c) % N]) ok = 0;
        if (!ok) continue;
        for (int k = 0; k <= KK; k++) {
            unsigned char *dst = Q + (size_t)k * N; memcpy(dst, P + (size_t)k * N, N);
            for (int t = 1; t <= k; t++) {
                unsigned char *src = P + (size_t)(k - t) * N;
                int sh = (int)(((long)t * c) % N);
                for (int x = 0; x < N; x++) if (src[x]) { dst[(x + sh) % N] = 1; dst[(x - sh + N) % N] = 1; }
            }
        }
        sol[j] = c;
        if (dfs(j + 1, Q)) { free(Q); return 1; }
    }
    free(Q); return 0;
}
/* returns 1 if free (solution in out), 0 if forced */
int free_vector(int n, int r, int K, const int *nu, int *out) {
    N = n; R = r; KK = K; EVEN = (n % 2 == 0); NU = nu; sol = out;
    if (n == 1) return 0;
    unsigned char *P = calloc((size_t)(K + 1) * n, 1);
    for (int k = 0; k <= K; k++) P[(size_t)k * n] = 1;
    int res = dfs(0, P); free(P); return res;
}
