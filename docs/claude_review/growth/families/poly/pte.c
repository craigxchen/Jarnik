// Exhaustive search of polynomial templates z(u) = prod_m (u - beta_m^{(+/-)}),
// beta_m = x_m + i y_m in Z[i], x in [0,X], y in [1,Y], multisets of size D with min x = 0
// (real translation normalises the template; the endpoint condition is translation invariant).
//
// A choice picks, for each distinct root beta of multiplicity e, a copies of beta and e-a of
// conj(beta).  Signature coordinate r: S_r = sum_m (2a_m - e_m) Im(beta_m^r) = Im p_r(choice)
// (p_r = r-th power sum of the chosen roots; Re p_r does not depend on the choice).
// Endpoint (bounded C) for even D: S_1..S_{D/2-1} agree on the fibre; then
//    C_base = (2/D) * (max S_{D/2} - min S_{D/2})            (no gcd improvement).
// If additionally S_{D/2} agrees, or D is odd and S_1..S_{(D-1)/2} agree, C -> 0.
//
// Output: for each fibre size k, the least C_base found and an example; and the largest
// C->0 fibre.  usage: pte D X Y [kmin_print]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int D, X, Y, NC;
static int cx[256], cy[256];
static ll W[256][12];  // Im(beta^r)
static int sel[32];

#define MAXK 64
static ll bestspread[MAXK + 1];
static int bestex[MAXK + 1][32];
static ll bestcount[MAXK + 1];
static int zero_best = 0;
static int zero_ex[32];
static ll nmulti = 0;

typedef struct {
    ll key[6];
    ll val;
} Item;
static Item items[1 << 12];

static int cmpitem(const void *a, const void *b) {
    const Item *p = a, *q = b;
    for (int r = 0; r < 6; r++) {
        if (p->key[r] < q->key[r]) return -1;
        if (p->key[r] > q->key[r]) return 1;
    }
    return 0;
}

static int nkey, half;

static void process(void) {
    // distinct roots with multiplicity
    int dr[32], de[32], nd = 0;
    for (int i = 0; i < D; i++) {
        if (nd > 0 && dr[nd - 1] == sel[i]) de[nd - 1]++;
        else { dr[nd] = sel[i]; de[nd] = 1; nd++; }
    }
    // enumerate choices
    int a[32];
    for (int i = 0; i < nd; i++) a[i] = 0;
    int n = 0;
    while (1) {
        Item it;
        memset(&it, 0, sizeof it);
        for (int r = 1; r <= nkey; r++) {
            ll s = 0;
            for (int i = 0; i < nd; i++) s += (ll)(2 * a[i] - de[i]) * W[dr[i]][r];
            it.key[r - 1] = s;
        }
        if (D % 2 == 0) {
            ll s = 0;
            for (int i = 0; i < nd; i++) s += (ll)(2 * a[i] - de[i]) * W[dr[i]][half];
            it.val = s;
        }
        items[n++] = it;
        int j = 0;
        while (j < nd && a[j] == de[j]) { a[j] = 0; j++; }
        if (j == nd) break;
        a[j]++;
    }
    qsort(items, n, sizeof(Item), cmpitem);
    int i = 0;
    while (i < n) {
        int j = i + 1;
        while (j < n && cmpitem(&items[i], &items[j]) == 0) j++;
        int k = j - i;
        if (k >= 2) {
            if (D % 2 == 0) {
                ll mn = items[i].val, mx = items[i].val;
                for (int t = i; t < j; t++) {
                    if (items[t].val < mn) mn = items[t].val;
                    if (items[t].val > mx) mx = items[t].val;
                }
                ll sp = mx - mn;
                if (k > MAXK) k = MAXK;
                bestcount[k]++;
                if (sp < bestspread[k]) {
                    bestspread[k] = sp;
                    memcpy(bestex[k], sel, sizeof(int) * D);
                }
                // C->0 subfibre: equal val as well
                int t = i;
                while (t < j) {
                    int u = t + 1;
                    // items within the fibre are not sorted by val; count equal vals by O(k^2)
                    (void)u;
                    t++;
                }
                for (int p = i; p < j; p++) {
                    int c = 0;
                    for (int q = i; q < j; q++) if (items[q].val == items[p].val) c++;
                    if (c > zero_best) { zero_best = c; memcpy(zero_ex, sel, sizeof(int) * D); }
                }
            } else {
                if (k > zero_best) { zero_best = k; memcpy(zero_ex, sel, sizeof(int) * D); }
                if (k > MAXK) k = MAXK;
                bestcount[k]++;
                if (bestspread[k] > 0) { bestspread[k] = 0; memcpy(bestex[k], sel, sizeof(int) * D); }
            }
        }
        i = j;
    }
    nmulti++;
}

static void rec(int pos, int start, int haszero) {
    if (pos == D) {
        if (haszero) process();
        return;
    }
    for (int c = start; c < NC; c++) {
        sel[pos] = c;
        rec(pos + 1, c, haszero || cx[c] == 0);
    }
}

int main(int argc, char **argv) {
    D = atoi(argv[1]);
    X = atoi(argv[2]);
    Y = atoi(argv[3]);
    int kprint = argc > 4 ? atoi(argv[4]) : 2;
    NC = 0;
    for (int x = 0; x <= X; x++)
        for (int y = 1; y <= Y; y++) {
            cx[NC] = x; cy[NC] = y;
            ll re = 1, im = 0;
            for (int r = 1; r <= 11; r++) {
                ll nre = re * x - im * y, nim = re * y + im * x;
                re = nre; im = nim;
                W[NC][r] = im;
            }
            NC++;
        }
    half = D / 2;
    nkey = (D % 2 == 0) ? half - 1 : (D - 1) / 2;
    if (nkey > 6) { fprintf(stderr, "D too large\n"); return 1; }
    for (int k = 0; k <= MAXK; k++) { bestspread[k] = (ll)4e18; bestcount[k] = 0; }
    rec(0, 0, 0);
    printf("D=%d X=%d Y=%d multisets=%lld\n", D, X, Y, nmulti);
    for (int k = kprint; k <= MAXK; k++) {
        if (bestcount[k] == 0) continue;
        printf("k=%d fibres=%lld ", k, bestcount[k]);
        if (D % 2 == 0) {
            // C_base = 2*spread/D
            printf("Cbase=%.6f roots:", 2.0 * bestspread[k] / D);
        } else printf("C->0 roots:");
        for (int i = 0; i < D; i++) printf(" %d+%di", cx[bestex[k][i]], cy[bestex[k][i]]);
        printf("\n");
    }
    printf("largest C->0 fibre: %d roots:", zero_best);
    for (int i = 0; i < D; i++) printf(" %d+%di", cx[zero_ex[i]], cy[zero_ex[i]]);
    printf("\n");
    return 0;
}
