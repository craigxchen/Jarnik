#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
struct Rec { int s,c,m; bool operator<(Rec const&o) const { return s<o.s || (s==o.s && c<o.c); } };
struct Hit { int M=0, fiber=0; vector<int>a, rows; int sum=0,cube=0; };
uint64_t choose_count=0, tuple_count=0, fiber_count=0, source_count=0;
Hit globalMax;
vector<Hit> hits;

string maskset(const vector<int>& rows){ string z="["; for(size_t k=0;k<rows.size();k++){if(k)z+=",";z+=to_string(rows[k]);} return z+"]"; }
string fracstr(uint64_t n,uint64_t d){ uint64_t g=gcd(n,d); n/=g; d/=g; return to_string(n)+"/"+to_string(d); }
void inspect(int M,const vector<int>&a,const vector<Rec>& rec,int lo,int hi){
    int f=hi-lo;
    vector<int> rows; for(int q=lo;q<hi;q++) rows.push_back(rec[q].m);
    if(f>globalMax.fiber){ globalMax={M,f,a,rows,rec[lo].s,rec[lo].c}; }
    if(f<5) return;
    fiber_count++;
    vector<int> fifth;
    for(int mask:rows){ int t=0; for(int j=0;j<10;j++) if(mask>>j&1){ int x=a[j]; int x2=x*x; int x4=x2*x2; t += x4*x; } fifth.push_back(t); }
    vector<int> ord(f); iota(ord.begin(),ord.end(),0);
    sort(ord.begin(),ord.end(),[&](int i,int j){return fifth[i]<fifth[j];});
    bool strict=true; for(int i=1;i<f;i++) if(fifth[ord[i]]==fifth[ord[i-1]]) strict=false;
    int source=0; vector<int> sourceCols;
    if(strict){
      for(int j=0;j<10;j++){ int b0=(rows[ord[0]]>>j)&1,b4=(rows[ord[f-1]]>>j)&1; bool mid=true; for(int k=1;k<f-1;k++) if(((rows[ord[k]]>>j)&1)==b0) mid=false; if(b0==b4 && mid){source++;sourceCols.push_back(j+1);} }
    }
    if(source){source_count++; if(hits.size()<20){ Hit h;h.M=M;h.fiber=f;h.a=a;h.rows=rows;h.sum=rec[lo].s;h.cube=rec[lo].c;hits.push_back(h);} }
}
void reccomb(int M,int at,int need,vector<int>&a){
 if(need==0){
   tuple_count++;
   static vector<int> sm(1024),cu(1024),fifth(1024); sm[0]=cu[0]=fifth[0]=0;
   for(int m=1;m<1024;m++){int b=__builtin_ctz(m), p=m&(m-1);int x=a[b],x2=x*x;sm[m]=sm[p]+x;cu[m]=cu[p]+x*x2;fifth[m]=fifth[p]+x*x2*x2*x;}
   vector<Rec> rec;rec.reserve(1024);for(int m=0;m<1024;m++)rec.push_back({sm[m],cu[m],m});sort(rec.begin(),rec.end());
   int lo=0;while(lo<1024){int hi=lo+1;while(hi<1024&&rec[hi].s==rec[lo].s&&rec[hi].c==rec[lo].c)hi++;inspect(M,a,rec,lo,hi);lo=hi;}
   if(tuple_count%10000==0) cerr<<"M="<<M<<" tuples="<<tuple_count<<" max="<<globalMax.fiber<<"\n";
   return;
 }
 for(int x=at;x<=M-need+1;x++){a.push_back(x);reccomb(M,x+1,need-1,a);a.pop_back();}
}
int main(int argc,char**argv){
 vector<int> Ms={10,12,14,16,18};
 for(int M:Ms){tuple_count=0;fiber_count=source_count=0;hits.clear();globalMax=Hit{};vector<int>a;reccomb(M,1,10,a);cout<<"M="<<M<<" tuples="<<tuple_count<<" maxFiber="<<globalMax.fiber<<" maxA=";for(int x:globalMax.a)cout<<x<<",";cout<<" maxRows="<<maskset(globalMax.rows)<<" maxKey="<<globalMax.sum<<","<<globalMax.cube<<" groups>=5="<<fiber_count<<" sourceHits="<<source_count<<"\n"; for(auto&h:hits){uint64_t prod=1,r5=0;for(int x:h.a){prod*=x;}for(int m:h.rows){int t=0;for(int j=0;j<10;j++)if(m>>j&1){uint64_t x=h.a[j];r5+=x*x*x*x*x;} }
   vector<int> vals;for(int m:h.rows){uint64_t t=0;for(int j=0;j<10;j++)if(m>>j&1){uint64_t x=h.a[j];t+=x*x*x*x*x;}vals.push_back((int)t);} auto mm=minmax_element(vals.begin(),vals.end()); uint64_t rg=*mm.second-*mm.first;cout<<" hit a=";for(int x:h.a)cout<<x<<",";cout<<" rows="<<maskset(h.rows)<<" range5="<<rg<<" C2="<<fracstr(4*rg*rg,25*prod)<<"\n"; }
 }
}
