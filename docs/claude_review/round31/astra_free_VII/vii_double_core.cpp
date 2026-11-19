// Exhaust the degree-e first Walsh core with eight marks in F2[t].
// c++ -O3 -std=c++17 vii_double_core.cpp -o vii_double_core
// ./vii_double_core [largest even e, default 16, maximum 16]
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using BinPoly=std::uint64_t;
using Shell=std::array<unsigned,4>;
using PairCode=std::uint32_t;

static int degree(BinPoly f) {return f ? 63-__builtin_clzll(f) : -1;}
static BinPoly remainder(BinPoly f,BinPoly g) {
    int d=degree(g);
    while(degree(f)>=d) f^=g<<(degree(f)-d);
    return f;
}
static BinPoly quotient(BinPoly f,BinPoly g) {
    BinPoly result=0;
    int d=degree(g);
    while(degree(f)>=d) {
        int k=degree(f)-d;
        result^=BinPoly(1)<<k;
        f^=g<<k;
    }
    return result;
}
static BinPoly multiply(BinPoly f,BinPoly g) {
    BinPoly result=0;
    while(g) {
        result^=f<<__builtin_ctzll(g);
        g&=g-1;
    }
    return result;
}
static std::vector<BinPoly> irreducibles(int e) {
    std::vector<BinPoly> small,primes;
    auto irreducible=[&](BinPoly f,int deg){
        for(auto p:small) {
            if(2*degree(p)>deg) break;
            if(!remainder(f,p)) return false;
        }
        return true;
    };
    for(int d=1;d<=e/2;++d)
        for(BinPoly f=BinPoly(1)<<d;f<BinPoly(2)<<d;++f)
            if(irreducible(f,d)) small.push_back(f);
    for(BinPoly f=BinPoly(1)<<e;f<BinPoly(2)<<e;++f)
        if(irreducible(f,e)) primes.push_back(f);
    return primes;
}
static BinPoly compress_even(BinPoly x) {
    x &= 0x5555555555555555ULL;
    x=(x|(x>>1))&0x3333333333333333ULL;
    x=(x|(x>>2))&0x0f0f0f0f0f0f0f0fULL;
    x=(x|(x>>4))&0x00ff00ff00ff00ffULL;
    x=(x|(x>>8))&0x0000ffff0000ffffULL;
    return (x|(x>>16))&0xffffffffULL;
}
static PairCode pair_code(unsigned a,unsigned b){return (a<<16)|b;}
static BinPoly shell_code(const Shell &a) {
    BinPoly result=0; for(auto i:a) result=(result<<16)|i;return result;
}
// Symmetric difference, with the two-shell test |A intersection B| == 2.
static bool add_shell(const Shell&a,const Shell&b,Shell&out) {
    unsigned merged[8];
    auto end=std::set_symmetric_difference(a.begin(),a.end(),b.begin(),b.end(),merged);
    if(end!=merged+4) return false;
    std::copy_n(merged,4,out.begin());return true;
}
struct PairSignature {unsigned signature;PairCode pair;};
static void run(int e) {
    const auto primes=irreducibles(e);
    const int b=(e/2-1)/2, L=2*e-b-1;
    std::vector<PairSignature> pairs;
    pairs.reserve(primes.size()*(primes.size()-1)/2);
    for(unsigned i=0;i<primes.size();++i)for(unsigned j=i+1;j<primes.size();++j) {
        BinPoly f=multiply(primes[i],primes[j]);
        BinPoly even=compress_even(f),odd=compress_even(f>>1);
        assert(degree(even)==e);
        // f=even^2+t odd^2. Two monic pair products give derivative
        // square root of degree <=b iff the Laurent quotients odd/even
        // agree through coefficient t^{-L} (including that coefficient).
        unsigned signature=quotient(odd<<L,even);
        pairs.push_back({signature,pair_code(i,j)});
    }
    std::sort(pairs.begin(),pairs.end(),[](auto x,auto y){
        return x.signature<y.signature;
    });
    std::vector<Shell> shells;
    std::unordered_map<BinPoly,unsigned> ids;
    std::uint64_t collisions=0;
    for(std::size_t begin=0;begin<pairs.size();) {
        std::size_t end=begin+1;
        while(end<pairs.size() && pairs[end].signature==pairs[begin].signature) ++end;
        for(std::size_t i=begin;i<end;++i)for(std::size_t j=i+1;j<end;++j) {
            auto p=pairs[i].pair,q=pairs[j].pair;
            Shell s={p>>16,p&65535,q>>16,q&65535};std::sort(s.begin(),s.end());
            if(std::adjacent_find(s.begin(),s.end())!=s.end()) continue;
            ++collisions;
            auto key=shell_code(s);
            if(ids.emplace(key,shells.size()).second) shells.push_back(s);
        }
        begin=end;
    }
    // Every shell is encountered via its three pairings.
    assert(collisions==3*shells.size());
    if(e==10) {
        // Independent check: every unordered quartet, with an ordinary
        // polynomial product and literal formal derivative.
        std::size_t direct_count=0;
        for(unsigned i=0;i<primes.size();++i)
        for(unsigned j=i+1;j<primes.size();++j) {
            const BinPoly pair=multiply(primes[i],primes[j]);
            for(unsigned k=j+1;k<primes.size();++k) {
                const BinPoly triple=multiply(pair,primes[k]);
                for(unsigned l=k+1;l<primes.size();++l) {
                    const BinPoly H=multiply(triple,primes[l]);
                    const BinPoly derivative=(H>>1)&0x5555555555555555ULL;
                    bool qualifies=degree(derivative)<=e/2-1;
                    assert(qualifies==(ids.count(shell_code({i,j,k,l}))!=0));
                    if(qualifies) ++direct_count;
                }
            }
        }
        assert(direct_count==shells.size());
    }
    std::unordered_map<PairCode,std::vector<unsigned>> with_pair;
    for(unsigned i=0;i<shells.size();++i) {
        const auto &a=shells[i];
        for(int j=0;j<4;++j)for(int k=j+1;k<4;++k)
            with_pair[pair_code(a[j],a[k])].push_back(i);
    }
    std::vector<std::array<unsigned,3>> lines;
    std::vector<std::vector<unsigned>> neighbours(shells.size());
    for(const auto &entry:with_pair) {
        const auto &v=entry.second;
        for(std::size_t i=0;i<v.size();++i)for(std::size_t j=i+1;j<v.size();++j) {
            unsigned a=v[i],bidx=v[j];Shell c;
            if(!add_shell(shells[a],shells[bidx],c)) continue;
            auto it=ids.find(shell_code(c));
            if(it==ids.end()) continue;
            unsigned ci=it->second;
            neighbours[a].push_back(bidx);neighbours[bidx].push_back(a);
            if(a<bidx && bidx<ci) lines.push_back({a,bidx,ci});
        }
    }
    std::uint64_t bases=0;
    auto in_shells=[&](unsigned x,unsigned y) {
        Shell s;return add_shell(shells[x],shells[y],s) && ids.count(shell_code(s))!=0;
    };
    for(auto line:lines) for(auto u:neighbours[line[0]]) {
        if(u==line[1] || u==line[2]) continue;
        if(in_shells(u,line[1]) && in_shells(u,line[2])) ++bases;
    }
    // In a full 3-space each of 7 lines has 4 possible outside words.
    assert(bases%28==0);
    std::printf("e=%d irreducibles=%zu pairs=%zu four_prime_shells=%zu lines=%zu full_simplex_cores=%llu\n",
                e,primes.size(),pairs.size(),shells.size(),lines.size(),
                (unsigned long long)(bases/28));
    std::fflush(stdout);
}
int main(int argc,char**argv) {
    int max_e=argc>1?std::atoi(argv[1]):16;
    assert(10<=max_e && max_e<=16 && max_e%2==0);
    for(int e=10;e<=max_e;e+=2)run(e);
}
