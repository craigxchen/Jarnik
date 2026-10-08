// Exact finite exclusion of all normalized odd-positive 8-row restrictions
// of the Paley12 tensor Sylvester2 matrix, after deleting the constant and
// any balanced column. Compile: c++ -O3 -std=c++17 this_file.cpp -o checker
// No floating point, external library, or assumed endpoint realization.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>
using namespace std;

using Cuts = array<uint8_t, 22>;
using Rows = array<int, 8>;
struct RawClass {
    Rows rows;
    int deleted;
    uint8_t deleted_mask;
    int code_count;
};
struct CanonicalClass {
    RawClass representative;
    int raw_count = 0;
    int code_count = 0;
};

static array<array<int, 24>, 24> tensor24() {
    array<array<int, 12>, 12> h12{};
    h12[0].fill(1);
    bool square[11]{};
    for (int i = 1; i < 11; ++i) square[i * i % 11] = true;
    for (int i = 0; i < 11; ++i) {
        h12[i + 1][0] = 1;
        for (int j = 0; j < 11; ++j)
            h12[i + 1][j + 1] = i == j ? -1 : (square[(i - j + 11) % 11] ? 1 : -1);
    }
    array<array<int, 24>, 24> h{};
    for (int i = 0; i < 12; ++i)
        for (int u = 0; u < 2; ++u)
            for (int j = 0; j < 12; ++j)
                for (int v = 0; v < 2; ++v)
                    h[2 * i + u][2 * j + v] = h12[i][j] * (u && v ? -1 : 1);
    for (int i = 0; i < 24; ++i)
        for (int j = 0; j < 24; ++j) {
            int dot = 0;
            for (int k = 0; k < 24; ++k) dot += h[i][k] * h[j][k];
            assert(dot == (i == j ? 24 : 0));
        }
    return h;
}

static uint8_t unoriented(int mask) { return min(mask, 255 ^ mask); }
static uint8_t transform(uint8_t mask, const array<int, 8>& perm) {
    int out = 0;
    for (int i = 0; i < 8; ++i) if (mask & (1 << i)) out |= 1 << perm[i];
    return unoriented(out);
}

static vector<array<int, 8>> partition_permutations() {
    vector<array<int, 8>> group;
    array<int, 4> a{{0, 1, 2, 3}};
    do {
        array<int, 4> b{{4, 5, 6, 7}};
        do {
            for (int swap_sides = 0; swap_sides < 2; ++swap_sides) {
                array<int, 8> p{};
                for (int i = 0; i < 4; ++i) p[i] = a[i] + 4 * swap_sides;
                for (int i = 0; i < 4; ++i) p[4 + i] = b[i] - 4 * swap_sides;
                group.push_back(p);
            }
        } while (next_permutation(b.begin(), b.end()));
    } while (next_permutation(a.begin(), a.end()));
    assert(group.size() == 1152);
    return group;
}

static Cuts canonical_key(const Cuts& cuts, int deleted_mask,
                          const vector<array<int, 8>>& group) {
    // Send the deleted balanced cut to {0,1,2,3}; its complement is
    // equivalent because column signs may be reversed. The residual group
    // permutes within the two four-row classes and may swap the classes.
    array<int, 8> base{};
    int lo = 0, hi = 4;
    for (int i = 0; i < 8; ++i)
        base[i] = deleted_mask & (1 << i) ? lo++ : hi++;
    assert(lo == 4 && hi == 8);
    Cuts standard{};
    for (int j = 0; j < 22; ++j) standard[j] = transform(cuts[j], base);
    Cuts best{};
    best.fill(255);
    for (const auto& p : group) {
        Cuts candidate{};
        for (int j = 0; j < 22; ++j)
            candidate[j] = transform(standard[j], p);
        sort(candidate.begin(), candidate.end());
        if (candidate < best) best = candidate;
    }
    return best;
}

struct CriticalSearch {
    array<array<int8_t, 22>, 16> c{};
    array<uint32_t, 16> support{};
    array<vector<int>, 22> incident{};
    unordered_set<uint64_t> seen;
    array<int, 23> by_depth{};
    int certificate_count = 0;
    static constexpr array<int, 11> pattern{{1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1}};

    CriticalSearch(const array<array<int, 24>, 24>& h,
                   const Rows& rows, int deleted) {
        for (int i = 0; i < 8; ++i)
            for (int k = i + 1; k < 8; ++k) {
                array<int8_t, 22> diff{};
                int used = 0, j = 0;
                for (int col = 1; col < 24; ++col) {
                    if (col == deleted) continue;
                    diff[j] = (h[rows[i]][col] - h[rows[k]][col]) / 2;
                    used += diff[j] != 0;
                    ++j;
                }
                assert(j == 22);
                if (used != 11) continue;
                assert(certificate_count < 16);
                int t = certificate_count++;
                c[t] = diff;
                for (int col = 0; col < 22; ++col)
                    if (diff[col]) {
                        support[t] |= 1u << col;
                        incident[col].push_back(t);
                    }
            }
        assert(certificate_count == 16);
    }

    bool feasible(uint32_t mask, uint16_t orientations) {
        if (mask == ((1u << 22) - 1)) return true;
        uint64_t state = (uint64_t(mask) << 16) | orientations;
        if (!seen.insert(state).second) return false;
        ++by_depth[__builtin_popcount(mask)];
        for (int col = 0; col < 22; ++col) {
            if (mask & (1u << col)) continue;
            int required = 0;
            bool conflict = false;
            for (int t : incident[col]) {
                int position = __builtin_popcount(mask & support[t]);
                if (!position) continue;
                int orientation = orientations & (1u << t) ? 1 : -1;
                int sign = orientation * pattern[position] * c[t][col];
                if (required && required != sign) { conflict = true; break; }
                required = sign;
            }
            if (conflict) continue;
            for (int sign : {1, -1}) {
                if (required && sign != required) continue;
                uint16_t next = orientations;
                for (int t : incident[col])
                    if (!(mask & support[t]) && c[t][col] * sign == 1)
                        next |= uint16_t(1u << t);
                if (feasible(mask | (1u << col), next)) return true;
            }
        }
        return false;
    }
};

int main() {
    auto h = tensor24();
    map<Cuts, RawClass> raw;
    Rows rows{{0, 1, 2, 3, 4, 5, 6, 7}};
    int row_sets = 0, codes = 0;
    while (true) {
        array<int, 24> masks{};
        for (int j = 0; j < 24; ++j)
            for (int i = 0; i < 8; ++i)
                if (h[rows[i]][j] < 0) masks[j] |= 1 << i;
        vector<int> balanced;
        bool has_odd = false;
        for (int j = 0; j < 24; ++j) {
            int n = __builtin_popcount(unsigned(masks[j]));
            if (n == 4) balanced.push_back(j);
            if (n & 1) has_odd = true;
        }
        if (balanced.size() >= 14 && has_odd) {
            assert(masks[0] == 0 && balanced.size() == 16);
            ++row_sets;
            for (int deleted : balanced) {
                Cuts key{};
                int k = 0;
                for (int j = 1; j < 24; ++j)
                    if (j != deleted) key[k++] = unoriented(masks[j]);
                assert(k == 22);
                sort(key.begin(), key.end());
                auto it = raw.find(key);
                if (it == raw.end())
                    raw.emplace(key, RawClass{rows, deleted,
                                              uint8_t(masks[deleted]), 1});
                else {
                    // The retained cut multiset determines the missing
                    // balanced partition in this enumerated family.
                    assert(unoriented(it->second.deleted_mask) ==
                           unoriented(masks[deleted]));
                    ++it->second.code_count;
                }
                ++codes;
            }
        }
        int i = 7;
        while (i >= 0 && rows[i] == 24 - 8 + i) --i;
        if (i < 0) break;
        ++rows[i];
        for (int j = i + 1; j < 8; ++j) rows[j] = rows[j - 1] + 1;
    }
    assert(row_sets == 1320 && codes == 21120 && raw.size() == 9056);
    auto group = partition_permutations();
    map<Cuts, CanonicalClass> classes;
    for (const auto& kv : raw) {
        Cuts key = canonical_key(kv.first, kv.second.deleted_mask, group);
        auto it = classes.find(key);
        if (it == classes.end())
            classes.emplace(key, CanonicalClass{kv.second, 1,
                                                kv.second.code_count});
        else {
            ++it->second.raw_count;
            it->second.code_count += kv.second.code_count;
        }
    }
    assert(classes.size() == 4);
    int total_raw = 0, total_codes = 0;
    for (const auto& kv : classes) {
        const auto& x = kv.second;
        total_raw += x.raw_count;
        total_codes += x.code_count;
        CriticalSearch search(h, x.representative.rows, x.representative.deleted);
        assert(!search.feasible(0, 0));
        int max_depth = 0;
        for (int d = 0; d <= 22; ++d) if (search.by_depth[d]) max_depth = d;
        cout << "deleted " << x.representative.deleted
             << " raw " << x.raw_count << " codes " << x.code_count
             << " states " << search.seen.size()
             << " max_depth " << max_depth << " rows";
        for (int r : x.representative.rows) cout << ' ' << r;
        cout << '\n';
    }
    assert(total_raw == 9056 && total_codes == 21120);
    cout << "PASS: 1320 row sets; 21120 deletions; 9056 raw keys; "
            "4 classes; all fail the critical 11-root pattern.\n";
}
