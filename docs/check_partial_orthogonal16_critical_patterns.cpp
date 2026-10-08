// Exhaust every ordered 8-by-14 signed-column sequence compatible with the
// critical/successor root orders, one common spectral ordering, and exact
// pair distances of an 8-by-16 orthogonal-row sign matrix punctured at a
// constant column and a balanced column. Compile with -std=c++17.
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
struct Edge{int i,j,type,orient,target;};
struct Key{uint64_t lo,hi;bool operator==(Key const& z)const{return lo==z.lo&&hi==z.hi;}};
struct Hash{size_t operator()(Key const& z)const{return size_t(z.lo*0x9e3779b97f4a7c15ULL ^ z.hi);}};
int pat(int type,int k){
 static int C[7]={1,1,-1,-1,1,1,-1};
 static int A[8]={1,-1,-1,1,1,-1,-1,1};
 static int B[8]={1,1,1,-1,-1,1,1,-1};
 return type==0?C[k]:type==1?A[k]:B[k];
}
struct Solver{
 vector<Edge> edges;array<unordered_set<Key,Hash>,15> seen;
 array<int,14> path{};uint64_t nodes=0;
 Solver(int low){
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){
   int li=i<low?0:(i<low+4?1:2), lj=j<low?0:(j<low+4?1:2);
   int type=li==lj?1:((li==0&&lj==2)?2:0);
   edges.push_back({i,j,type,type==1?-1:1,type==0?7:8});
  }
 }
 int val(Key k,int e){return e<16?(k.lo>>(4*e))&15:(k.hi>>(4*(e-16)))&15;}
 void inc(Key &k,int e){if(e<16)k.lo+=uint64_t(1)<<(4*e);else k.hi+=uint64_t(1)<<(4*(e-16));}
 bool dfs(int depth,Key key){
  if(depth==14){for(int e=0;e<28;e++)if(val(key,e)!=edges[e].target)return false;return true;}
  if(seen[depth].find(key)!=seen[depth].end())return false;
  seen[depth].insert(key);nodes++;
  // mask=0 covers a retained constant column. mask=255 is its global
  // negative and has the same pairwise effect, so one no-op is enough.
  // Other repeated columns are allowed.
  for(int mask=0;mask<255;mask++){
   Key next=key;bool okay=true;
   for(int e=0;e<28;e++){
    auto x=edges[e];
    int a=(mask>>x.i&1)?1:-1,b=(mask>>x.j&1)?1:-1;
    int c=(a-b)/2,k=val(key,e);
    if(c){
     if(k==x.target||c!=x.orient*pat(x.type,k)){okay=false;break;}
     inc(next,e);k++;
    }
    if(k+13-depth<x.target){okay=false;break;}
   }
   if(!okay)continue;
   path[depth]=mask;
   if(dfs(depth+1,next))return true;
  }
  return false;
 }
};
int main(){
 const uint64_t expected[5]={59684,87000,77636,87000,59684};
 for(int low=0;low<=4;low++){
  Solver s(low);bool found=s.dfs(0,{0,0});
  assert(!found && s.nodes==expected[low]);
  cout<<"low B rows "<<low<<": "<<s.nodes<<" exact dead states\n";
  if(found){for(int x:s.path)cout<<x<<' ';cout<<"\n";}
 }
 cout<<"PASS: no punctured orthogonal 8-by-16 pattern assignment\n";
}
