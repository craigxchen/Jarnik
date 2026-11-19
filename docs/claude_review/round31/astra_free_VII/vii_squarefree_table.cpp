// c++ -O3 -std=c++17 vii_squarefree_table.cpp -o /path/in/cwd/vii_squarefree_table
// ./vii_squarefree_table [maximum degree, default 16]
// F9 coefficient a+ib is a+3b (0<=a,b<3); F9[t] polynomials use base 9.
//
// Distinct points of one cluster have the same leading coefficient. Divide all
// by that coefficient to normalize the points and their parent to monic.
// A squarefree parent splits into conjugate pairs of monic irreducibles of
// F9[t]. Each point chooses one prime from every pair. Among any five points,
// each pair of primes has uniform orientation on at least one of the ten
// triples. Count these contributions weighted by deg(pi); some triple shares
// a common factor G of degree >= ceil(D/10). On such a triple write
// z_i=G U_i prod_{j!=i} bar(U_j): each nonconstant orientation column has
// exactly one exceptional mark i, whose prime goes in U_i. The norms of the
// four blocks are squarefree and disjoint.
//
// Let K=ceil(D/2). The ratio of two points in one cluster is 1+O(t^{-K}).
// The condition on two of the above marks is that U_i bar(U_j) is real up
// to depth K relative to its leading monomial. It is equivalent to equality
// of the slopes Im(U_i)/Re(U_i) and Im(U_j)/Re(U_j) truncated mod t^{-K}.
// Unless the marks coincide one has deg U_i + deg U_j >= K. Thus
// max_i deg U_i <= D-ceil(D/10)-K and deg G <= D-ceil(3K/2).
//
// The loops enumerate all such triples and G; the subset phase count counts
// every independent choice of conjugates on the resulting parent producing
// the triple's prefix. For D=1..5 a direct parent census is smaller. When
// the five-point covering criterion does not yet give a four-point witness,
// we also enumerate triples with unrestricted common block. This covers all
// sets of size at least three. All three modes output a witness via explicit
// polynomial multiplication, independent of the truncated-angle group law.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

using Word = std::uint64_t;
using Poly = std::vector<int>;
static int sum9[9][9], prod9[9][9], neg9[9], conj9[9];
static Word power9[19];
static int power3[17];

static Poly unpack(Word w) {
    Poly p;
    while (w) { p.push_back(w%9); w/=9; }
    return p;
}
static Word pack(const Poly &p) {
    Word w=0;
    for (auto it=p.rbegin(); it!=p.rend(); ++it) w=9*w+*it;
    return w;
}
static Word conjugate(Word f) {
    Poly p=unpack(f);
    for (int &c:p) c=conj9[c];
    return pack(p);
}
static Word multiply(Word f, Word g) {
    Poly p=unpack(f), q=unpack(g), u(p.size()+q.size()-1,0);
    for (unsigned i=0;i<p.size();++i)
        for (unsigned j=0;j<q.size();++j)
            u[i+j]=sum9[u[i+j]][prod9[p[i]][q[j]]];
    return pack(u);
}
static Word exact_quotient(const Poly &f, const Poly &g) {
    Poly r=f, q(r.size()-g.size()+1,0);
    for(int i=q.size()-1;i>=0;--i) {
        q[i]=r[i+g.size()-1];
        if(q[i])
            for(unsigned j=0;j<g.size();++j)
                r[i+j]=sum9[r[i+j]][neg9[prod9[q[i]][g[j]]]];
    }
    for(unsigned j=0;j+1<g.size();++j) if(r[j]) return 0;
    return pack(q);
}
static void tables() {
    for(int x=0;x<9;++x) {
        neg9[x]=(3-x%3)%3+3*((3-x/3)%3);
        conj9[x]=x%3+3*((3-x/3)%3);
        for(int y=0;y<9;++y) {
            sum9[x][y]=(x%3+y%3)%3+3*((x/3+y/3)%3);
            prod9[x][y]=(9+(x%3)*(y%3)-(x/3)*(y/3))%3+
                        3*(((x%3)*(y/3)+(x/3)*(y%3))%3);
        }
    }
    power3[0]=1; for(int i=1;i<17;++i) power3[i]=3*power3[i-1];
    power9[0]=1; for(int i=1;i<19;++i) power9[i]=9*power9[i-1];
}
struct Block {
    Word poly;
    int degree;
    std::vector<Word> factors;
    std::vector<Word> base;  // min(pi,bar pi) for each split prime
    int angle;
    std::vector<int> subsets;
};
struct Factors {
    std::map<Word, std::vector<Word>> fact;
    std::vector<Block> valid[8];
    Factors(int max_degree) {
        struct Irred { Word code; Poly coeff; };
        std::vector<std::vector<Irred>> irred(max_degree+1);
        fact[1]={};
        valid[0].push_back({1,0,{}, {},0,{}});
        for(int d=1;d<=max_degree;++d) for(Word n=power9[d]; n<2*power9[d];++n) {
            Poly p=unpack(n);
            Word divisor=0, quotient=0;
            for(int e=1;e<=d/2 && !divisor;++e) for(const auto &v:irred[e]) {
                quotient=exact_quotient(p,v.coeff);
                if(quotient) {divisor=v.code;break;}
            }
            std::vector<Word> primes;
            if(divisor) { primes=fact.at(quotient); primes.push_back(divisor); }
            else { irred[d].push_back({n,p}); primes.push_back(n); }
            fact.emplace(n,primes);
            bool ok=true; std::vector<Word> keys;
            for(Word f:primes) {
                Word fb=conjugate(f), key=std::min(f,fb);
                if(f==fb || std::find(keys.begin(),keys.end(),key)!=keys.end())
                    ok=false;
                keys.push_back(key);
            }
            if(ok) valid[d].push_back({n,d,primes,keys,0,{}});
        }
    }
};
static bool disjoint(const Block &a,const Block &b) {
    for(Word i:a.base) for(Word j:b.base) if(i==j) return false;
    return true;
}
// The angle is B/A modulo x^K for the reversed monic f=A+iB.
// The group law induced by multiplication is (u+v)/(1-uv).
struct Angles {
    int K, size;
    std::vector<std::vector<short>> add;
    std::vector<int> neg;
    Angles(int k):K(k),size(power3[k-1]),add(size,std::vector<short>(size)),neg(size) {
        std::vector<Poly> digit(size,Poly(K,0));
        for(int u=0;u<size;++u){
            int w=u, n=0;
            for(int j=1;j<K;++j){digit[u][j]=w%3; w/=3;n+=((3-digit[u][j])%3)*power3[j-1];}
            neg[u]=n;
        }
        for(int u=0;u<size;++u) for(int v=0;v<=u;++v) {
            int product[16]={}, result[16]={}, code=0;
            for(int j=1;j<K;++j) for(int h=1;h<j;++h)
                product[j]+=digit[u][h]*digit[v][j-h];
            for(int j=1;j<K;++j) {
                int b=digit[u][j]+digit[v][j];
                for(int h=1;h<j;++h) b+=product[h]*result[j-h];
                result[j]=b%3; code+=power3[j-1]*result[j];
            }
            add[u][v]=add[v][u]=code;
        }
    }
    int angle(Word n) const {
        Poly p=unpack(n);
        int d=p.size()-1, digit[16]={}, code=0;
        for(int j=1;j<K;++j) {
            int b=j<=d?p[d-j]/3:0;
            for(int h=1;h<=std::min(d,j);++h) b-=p[d-h]%3*digit[j-h];
            digit[j]=(b%3+3)%3;code+=digit[j]*power3[j-1];
        }
        return code;
    }
    void prepare(Block &b) const {
        b.angle=angle(b.poly);b.subsets={0};
        for(Word f:b.factors) {
            int a=angle(f), old=b.subsets.size();
            for(int j=0;j<old;++j) b.subsets.push_back(add[b.subsets[j]][a]);
        }
        assert(b.subsets.back()==b.angle);
    }
};
struct Witness {
    int best=0; Word top=0; std::vector<Word> points;
    // Exact audit of each improved lower bound using ordinary F9 products.
    void record(const std::vector<Word> &primes, int D, int count) {
        std::map<Word,std::vector<Word>> by_top;
        for(unsigned s=0;s<(1U<<primes.size());++s) {
            Word f=1;
            for(unsigned i=0;i<primes.size();++i)
                f=multiply(f,(s>>i&1)?primes[i]:conjugate(primes[i]));
            assert(f/power9[D]==1);
            by_top[f/power9[D/2+1]].push_back(f);
        }
        for(auto &entry:by_top) if((int)entry.second.size()>best) {
            best=entry.second.size(); top=entry.first; points=entry.second;
        }
        assert(best>=count);
    }
    void print(int D) {
        std::printf("D=%d M=%d base9_marks=",D,best);
        for(auto z:points) std::printf("%llu,",(unsigned long long)z);
    }
};
// Small direct census of all squarefree parents by their prime supports.
static void small_parent(int D,int max_degree_seen,Factors &f,Witness &w,
                         int remain,std::vector<Word> &primes,Word last) {
    if(!remain) {w.record(primes,D,0);return;}
    for(int e=max_degree_seen;e<=remain;++e) for(const auto &b:f.valid[e]) {
        if(b.factors.size()!=1 || b.poly>conjugate(b.poly)) continue;
        if(e==max_degree_seen && b.poly<=last) continue;
        primes.push_back(b.poly);
        small_parent(D,e,f,w,remain-e,primes,b.poly);
        primes.pop_back();
    }
}
struct Search {
    int D, K;
    Factors &f;
    Angles group;
    Witness w;
    std::uint64_t triples=0, completers=0;
    std::vector<int> hist, changed;
    Search(int d, Factors &fs):D(d),K((d+1)/2),f(fs),group(K),hist(group.size,0) {}
    void run(int common_floor, int max_degree) {
        std::vector<std::vector<std::vector<Block*>>> bucket(
            max_degree+1,std::vector<std::vector<Block*>>(group.size));
        for(int a=0;a<=max_degree;++a) for(auto &v:f.valid[a]) {
            group.prepare(v);bucket[a][v.angle].push_back(&v);
        }
        for(int a=0;a<=max_degree;++a) for(int b=a;b<=max_degree;++b)
        for(int c=b;c<=max_degree;++c) {
            int g=D-a-b-c;
            if(a+b<K || g<common_floor || g>max_degree) continue;
            for(auto &u:f.valid[a]) for(auto v:bucket[b][u.angle]) {
                if((a==b && u.poly>=v->poly) || !disjoint(u,*v)) continue;
                for(auto t:bucket[c][u.angle]) {
                    if((b==c && v->poly>=t->poly) || !disjoint(*v,*t) || !disjoint(u,*t))
                        continue;
                    ++triples;
                    for(int x:u.subsets) for(int y:v->subsets) {
                        int xy=group.add[x][y];
                        for(int z:t->subsets) {
                            int h=group.add[xy][z];
                            if(hist[h]++==0) changed.push_back(h);
                        }
                    }
                    for(auto &G:f.valid[g]) {
                        if(!disjoint(u,G)||!disjoint(*v,G)||!disjoint(*t,G)) continue;
                        ++completers;
                        int count=0;
                        // Positive primes of the common block are taken by the
                        // starting triple. Switching some subtracts their phases.
                        for(int y:G.subsets) count+=hist[group.add[u.angle][y]];
                        if(count>w.best) {
                            std::vector<Word> all=u.factors;
                            for(auto block:{v,t,&G})
                                all.insert(all.end(),block->factors.begin(),block->factors.end());
                            w.record(all,D,count);
                        }
                    }
                    for(int h:changed) hist[h]=0;
                    changed.clear();
                }
            }
        }
    }
};
int main(int argc,char**argv) {
    int limit=argc>1?std::atoi(argv[1]):16;
    assert(limit>=1&&limit<=16);tables();
    const int factor_limit=std::min(limit,limit==16 ? 6 : 5);
    Factors factors(factor_limit);
    for(int D=1;D<=limit;++D) {
        if(D<=5){
            Witness w;std::vector<Word> block;
            small_parent(D,1,factors,w,D,block,0);
            w.print(D); std::puts(" method=all_parents");
        } else {
            Search search(D,factors);
            int gmin=(D+9)/10, max_degree=D-gmin-search.K;
            // From A+B>=K, C<=D-gmin-K. Also G<=D-ceil(3K/2).
            max_degree=std::max(max_degree,D-(3*search.K+1)/2);
            assert(max_degree<=factor_limit);
            search.run(gmin,max_degree);
            if(search.w.best<4) {
                // Any 4-cluster has a triple; allow all triples here.
                int mx=D-search.K;
                if(mx>factor_limit) {std::printf("D=%d requires an additional four-point lower witness\n",D);return 2;}
                search.run(0,mx);
            }
            assert(search.w.best>=3);
            search.w.print(D);
            std::printf(" triples=%llu completers=%llu\n",(unsigned long long)search.triples,
                        (unsigned long long)search.completers);
        }
        std::fflush(stdout);
    }
}
