// Enumerate ordered relation flats. Modulo 769 all relevant minors have
// the same rank as over Q by the Hadamard bound in Part II.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

constexpr int prime = 769;
using Row = std::array<int, 8>;
std::array<int, prime> inverse;
std::vector<Row> equations;
std::array<long long, 6> nodes{}, heavy{};

struct Basis {
    std::array<Row, 8> rows{};
    int rank = 0;

    Row reduce(Row v) const {
        for (int j = 0; j < 8; ++j) {
            if (!v[j] || !rows[j][j]) continue;
            int complement = prime-v[j];
            for (int k = j; k < 8; ++k)
                v[k] = (v[k] + complement*rows[j][k]) % prime;
        }
        return v;
    }

    void append(Row v) {
        v = reduce(v);
        auto pivot = std::find_if(v.begin(), v.end(), [](int x) { return x != 0; });
        assert(pivot != v.end());
        int j = static_cast<int>(pivot-v.begin()), scale = inverse[*pivot];
        for (int k = j; k < 8; ++k)
            rows[j][k] = v[k]*scale % prime;
        ++rank;
    }
};

bool allowed(const Basis& basis) {
    std::array<Row,8> points;
    for (int i = 0; i < 8; ++i) {
        Row label{};
        label[i] = 1;
        points[i] = basis.reduce(label);
        for (int j = 0; j < i; ++j)
            if (points[i] == points[j]) return false;
    }
    for (int i = 0; i < 8; ++i) {
        for (int j = i+1; j < 8; ++j) {
            Row direction;
            int pivot = 0;
            for (int a = 0; a < 8; ++a) {
                direction[a] = (points[i][a]+prime-points[j][a]) % prime;
                if (direction[a]) pivot=a;
            }
            int extra = 0;
            for (int t = j+1; t < 8; ++t) {
                bool collinear = true;
                for (int a = 0; a < 8; ++a) {
                    int minor = (points[t][a]-points[j][a])*direction[pivot]
                              - (points[t][pivot]-points[j][pivot])*direction[a];
                    if (minor % prime) { collinear=false; break; }
                }
                if (collinear && ++extra >= 2) return false;
            }
        }
    }
    return true;
}

void enumerate(const Basis& basis, int last, std::vector<int>& ids, int target) {
    ++nodes[basis.rank];
    int contained = 0;
    std::vector<Row> quotient_directions;
    std::vector<std::vector<int>> groups;
    for (size_t i = 0; i < equations.size(); ++i) {
        Row v = basis.reduce(equations[i]);
        auto pivot = std::find_if(v.begin(), v.end(), [](int x) { return x != 0; });
        if (pivot == v.end()) { ++contained; continue; }
        int scale = inverse[*pivot];
        for (auto it = pivot; it != v.end(); ++it) *it = *it*scale % prime;
        auto group = std::find(quotient_directions.begin(), quotient_directions.end(), v);
        size_t index = group-quotient_directions.begin();
        if (group == quotient_directions.end()) {
            quotient_directions.push_back(v);
            groups.emplace_back();
        }
        groups[index].push_back(static_cast<int>(i));
    }
    if (contained >= target) {
        ++heavy[basis.rank];
        std::cout << basis.rank << ' ' << contained;
        for (int i: ids) std::cout << ' ' << i;
        std::cout << '\n';
    }
    if (basis.rank == 5) return;

    for (const auto& group: groups) {
        int next = group.front();
        if (next <= last) continue;
        if (basis.rank == 4 && contained+static_cast<int>(group.size()) < target) continue;
        Basis child=basis;
        child.append(equations[next]);
        if (child.rank >= 2 && !allowed(child)) continue;
        ids.push_back(next);
        enumerate(child, next, ids, target);
        ids.pop_back();
    }
}

int main(int argc, char** argv) {
    int target = argc>1 ? std::atoi(argv[1]) : 12;
    inverse[1]=1;
    for (int i=2; i<prime; ++i) inverse[i] = prime-(prime/i)*inverse[prime%i]%prime;
    for (int i=0; i<8; ++i) for (int j=i+1; j<8; ++j)
        for (int k=j+1; k<8; ++k) for (int l=k+1; l<8; ++l) {
            Row a{};
            a[i] = a[l] = 1;
            a[j] = a[k] = prime-1;
            equations.push_back(a);
        }
    if (argc>2 && std::string(argv[2])=="all")
        for (int i=0; i<8; ++i) for (int j=i+1; j<8; ++j) for (int k=j+1; k<8; ++k) {
            Row a{};
            a[i] = a[k] = 1;
            a[j] = prime-2;
            equations.push_back(a);
        }
    std::vector<int> ids;
    enumerate(Basis{}, -1, ids, target);
    for (int i=0; i<=5; ++i)
        std::cerr << i << " nodes " << nodes[i] << " heavy " << heavy[i] << '\n';
}
