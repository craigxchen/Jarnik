#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using Clock=chrono::steady_clock;
struct State{
 array<int8_t,8> r{};array<uint8_t,8> order{};
 bool operator==(State const& z)const{return r==z.r&&order==z.order;}
};
struct Hash{
 size_t operator()(State const& s)const{
  uint64_t h=1469598103934665603ULL;
  for(auto x:s.r){h^=uint8_t(x+2);h*=1099511628211ULL;}
  for(auto x:s.order){h^=x;h*=1099511628211ULL;}
  return size_t(h);
 }
};
struct Pair {uint8_t i,j,kind;int8_t ks,orientation,distance;}; // 0=C,1=A,2=B
struct Transition {uint8_t mask;State out;uint32_t inc;};
struct Edge {uint32_t dest;uint32_t inc;uint8_t mask;};

struct Machine{
 array<int8_t,8> delta{{0,0,0,0,-1,-1,1,1}};
 array<uint8_t,8> kappa{{6,7,0,1,2,3,4,5}};
 array<Pair,28> pair{};
 State initial{},terminal{};
 Machine(){
  int rank[8]{};for(int k=0;k<8;k++)rank[kappa[k]]=k;
  int e=0;
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){
   int kind=(i<4)!=(j<4)?0:(i<4||((i<6)==(j<6))?1:2);
   int ks=rank[i]>rank[j]?1:-1;
   pair[e++]={uint8_t(i),uint8_t(j),uint8_t(kind),int8_t(ks),int8_t(kind==1?ks:-ks),int8_t(kind==0?1:2)};
  }
  assert(e==28);
  array<int8_t,28> comp{};
  for(int k=0;k<28;k++){
   auto p=pair[k];int exponent=(delta[p.i]+delta[p.j]-p.distance)/2;
   comp[k]=int8_t(p.ks*((exponent&1)?-1:1));
  }
  assert(order_from_comparisons(comp,initial.order));
  terminal.r=delta;terminal.order=kappa;
  assert(valid(initial)&&valid(terminal));
 }
 static bool order_from_comparisons(array<int8_t,28> const& c,array<uint8_t,8>& order){
  int below[8]{};int e=0;
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,e++)below[c[e]==1?i:j]++;
  bool used[8]{};
  for(int i=0;i<8;i++){
   if(below[i]<0||below[i]>7||used[below[i]])return false;
   used[below[i]]=true;order[below[i]]=uint8_t(i);
  }
  int rank[8]{};for(int i=0;i<8;i++)rank[order[i]]=i;
  e=0;for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,e++)
   if(c[e]!=(rank[i]>rank[j]?1:-1))return false;
  return true;
 }
 static int mod4(int x){x%=4;return x<0?x+4:x;}
 int phase(State const& s,int e)const{
  auto p=pair[e];int rank[8]{};for(int k=0;k<8;k++)rank[s.order[k]]=k;
  int comp=rank[p.i]>rank[p.j]?1:-1;
  int mismatch=comp!=p.ks;
  return mod4(p.distance+s.r[p.i]+s.r[p.j]-delta[p.i]-delta[p.j]+2*mismatch);
 }
 bool valid(State const& s)const{
  if(s.r[0])return false;
  bool used[8]{};for(int k=0;k<8;k++){
   if(s.order[k]>7||used[s.order[k]])return false;
   used[s.order[k]]=true;
  }
  int amin=2,amax=-2;
  for(int a=0;a<4;a++){
   if(s.r[a]<-1||s.r[a]>1)return false;
   amin=min(amin,int(s.r[a]));amax=max(amax,int(s.r[a]));
  }
  if(amax-amin>1)return false;
  for(int b=4;b<8;b++){
   if(s.r[b]!=0&&s.r[b]!=delta[b]&&s.r[b]!=2*delta[b])return false;
  }
  int prev=-1;
  for(int k=0;k<8;k++){
   int i=s.order[k],w=((delta[i]-s.r[i])&1)?-1:1;
   if(w<prev)return false;prev=w;
  }
  for(int e=0;e<28;e++){
   auto p=pair[e];int h=p.orientation*(s.r[p.i]-s.r[p.j]),k=phase(s,e);
   if(p.kind==0){int want[4]={0,1,2,1};if(h!=want[k])return false;}
   else if(p.kind==1){int want[4]={0,1,0,-1};if(h!=want[k])return false;}
   else{int want[4]={2,1,2,3};if(!((h==0&&k==0)||h==want[k]))return false;}
  }
  return true;
 }
 vector<Transition> transitions(State const& s)const{
  int zero=0;for(int k=0;k<8;k++){
   int i=s.order[k];if((delta[i]-s.r[i])&1)zero++;
  }
  vector<uint8_t> masks{0,255};
  for(int left=0;left<=zero;left++)for(int right=zero;right<=8;right++)
   if(left<right){
    int mask=0;for(int k=left;k<right;k++)mask|=1<<s.order[k];
    if(mask!=255)masks.push_back(uint8_t(mask));
   }
  sort(masks.begin(),masks.end());
  masks.erase(unique(masks.begin(),masks.end()),masks.end());
  vector<Transition> out;
  for(auto mask:masks){
   State t=s;
   if(mask!=0&&mask!=255){
    int left=8,right=0;
    for(int k=0;k<8;k++)if(mask&(1<<s.order[k])){left=min(left,k);right=max(right,k+1);}
    reverse(t.order.begin()+left,t.order.begin()+right);
    if(mask&1)reverse(t.order.begin(),t.order.end());
    int b0=mask&1;
    for(int i=0;i<8;i++)t.r[i]+=((mask>>i)&1)-b0;
    if(!valid(t))continue;
   }
   uint32_t inc=0;
   for(int e=0;e<28;e++){
    auto p=pair[e];if(((mask>>p.i)&1)!=((mask>>p.j)&1))inc|=1u<<e;
   }
   out.push_back({mask,t,inc});
  }
  return out;
 }
};

int main(int argc,char**argv){
 Machine m;
 if(argc>1&&string(argv[1])=="--sample"){
  vector<State> sample_states{m.initial};
  unordered_map<State,uint32_t,Hash> sample_ids;sample_ids.emplace(m.initial,0);
  for(size_t at=0;at<sample_states.size();at++){
   for(auto const& t:m.transitions(sample_states[at]))
    if(sample_ids.emplace(t.out,uint32_t(sample_states.size())).second)
     sample_states.push_back(t.out);
  }
  assert(sample_states.size()==322564);
  for(int k:{0,1,7,42,100,500,1000,4999,10000,50000,100000,150000,200000,250000,300000,322563}){
   auto s=sample_states[k];auto tr=m.transitions(s);
   cout<<"state ";for(auto x:s.r)cout<<int(x)<<' ';for(auto x:s.order)cout<<int(x)<<' ';cout<<tr.size()<<'\n';
   for(auto const& t:tr){cout<<"edge "<<int(t.mask)<<' '<<t.inc<<' ';
    for(auto x:t.out.r)cout<<int(x)<<' ';for(auto x:t.out.order)cout<<int(x)<<' ';cout<<'\n';
   }
  }
  return 0;
 }
 constexpr size_t STATE_CAP=1000000, EDGE_CAP=30000000;
 constexpr int SECOND_CAP=60;
 auto start=Clock::now();
 vector<State> states{m.initial};
 unordered_map<State,uint32_t,Hash> ids;ids.emplace(m.initial,0);
 vector<uint32_t> parent{0};vector<uint8_t> parentmask{0};vector<uint16_t> depth{0};
 vector<Edge> edges;vector<size_t> offsets{0};
 int terminal_id=m.initial==m.terminal?0:-1;
 bool capped=false;string reason;
 for(size_t at=0;at<states.size();at++){
  if(at%1000==0){
   auto sec=chrono::duration_cast<chrono::seconds>(Clock::now()-start).count();
   if(sec>SECOND_CAP){capped=true;reason="time";break;}
  }
  auto tr=m.transitions(states[at]);
  for(auto const& t:tr){
   auto it=ids.find(t.out);uint32_t dest;
   if(it==ids.end()){
    dest=uint32_t(states.size());ids.emplace(t.out,dest);states.push_back(t.out);
    parent.push_back(uint32_t(at));parentmask.push_back(t.mask);
    depth.push_back(uint16_t(depth[at]+1));
    if(t.out==m.terminal&&terminal_id<0)terminal_id=int(dest);
   }else dest=it->second;
   edges.push_back({dest,t.inc,t.mask});
  }
  offsets.push_back(edges.size());
  if(states.size()>STATE_CAP||edges.size()>EDGE_CAP){capped=true;reason=states.size()>STATE_CAP?"states":"edges";break;}
 }
 cout<<"states "<<states.size()<<" edges "<<edges.size()<<" processed "<<offsets.size()-1
     <<" closed "<<(!capped)<<" reason "<<(capped?reason:"none")<<'\n';
 if(terminal_id>=0){
  vector<int> masks;int x=terminal_id;
  while(x!=0){masks.push_back(parentmask[x]);x=parent[x];}
  reverse(masks.begin(),masks.end());
  cout<<"terminal_depth "<<masks.size()<<" masks";
  for(int mask:masks)cout<<' '<<mask;
  cout<<'\n';
 }else cout<<"terminal_unseen\n";
 if(capped)return 2;
 assert(states.size()==322564&&edges.size()==2903076&&terminal_id>=0&&depth[terminal_id]==11);
 // Iterative Tarjan SCC on the closed graph.
 int n=states.size(),tick=0,components=0,largest=0;
 vector<int> idx(n,-1),low(n),comp(n,-1),node_stack;vector<uint8_t> on(n,0);
 struct Frame{int v;size_t next;};vector<Frame> call;
 for(int root=0;root<n;root++)if(idx[root]<0){
  idx[root]=low[root]=tick++;node_stack.push_back(root);on[root]=1;
  call.push_back({root,offsets[root]});
  while(!call.empty()){
   auto &fr=call.back();int v=fr.v;
   if(fr.next<offsets[v+1]){
    int w=edges[fr.next++].dest;
    if(idx[w]<0){idx[w]=low[w]=tick++;node_stack.push_back(w);on[w]=1;call.push_back({w,offsets[w]});}
    else if(on[w])low[v]=min(low[v],idx[w]);
   }else{
    if(low[v]==idx[v]){
     int size=0,w;
     do{w=node_stack.back();node_stack.pop_back();on[w]=0;comp[w]=components;size++;}while(w!=v);
     largest=max(largest,size);components++;
    }
    call.pop_back();if(!call.empty())low[call.back().v]=min(low[call.back().v],low[v]);
   }
  }
 }
 vector<uint8_t> closed(components,1);vector<int> compsize(components,0);
 for(int v=0;v<n;v++){
  compsize[comp[v]]++;
  for(size_t k=offsets[v];k<offsets[v+1];k++)
   if(comp[v]!=comp[edges[k].dest])closed[comp[v]]=0;
 }
 int closed_count=0,closed_states=0;
 for(int c=0;c<components;c++)if(closed[c]){closed_count++;closed_states+=compsize[c];}
 cout<<"scc_count "<<components<<" largest "<<largest<<" closed_scc "<<closed_count
     <<" closed_states "<<closed_states<<'\n';
 assert(components==5&&largest==322560&&closed_count==1&&closed_states==322560);
 cout<<"initial_scc_size "<<compsize[comp[0]]<<" terminal_scc_size "
     <<(terminal_id>=0?compsize[comp[terminal_id]]:0)<<'\n';
 int giant=-1;for(int c=0;c<components;c++)if(compsize[c]==largest)giant=c;
 assert(giant>=0&&largest==322560&&closed[giant]);
 int binom[9]={1,8,28,56,70,56,28,8,1};
 vector<uint32_t> indegree(size_t(n)*7,0);
 array<int,256> label_count{};
 bool nine_edges=true,one_each=true,permutation=true,uniform_labels=true;
 for(int v=0;v<n;v++){
  if(offsets[v+1]-offsets[v]!=9)nine_edges=false;
  int by_size[9]{};
  for(size_t e=offsets[v];e<offsets[v+1];e++){
   auto x=edges[e];int k=__builtin_popcount(unsigned(x.mask));
   by_size[k]++;
   if(comp[v]!=giant)continue;
   label_count[x.mask]++;
   if(k>=1&&k<=7){
    if(comp[x.dest]!=giant)permutation=false;
    ++indegree[size_t(k-1)*n+x.dest];
   }
  }
  for(int k=1;k<=7;k++)if(by_size[k]!=1)one_each=false;
  if(by_size[0]!=1||by_size[8]!=1)one_each=false;
 }
 for(int v=0;v<n;v++)if(comp[v]==giant)
  for(int k=1;k<=7;k++)if(indegree[size_t(k-1)*n+v]!=1)permutation=false;
 for(int mask=0;mask<256;mask++){
  int k=__builtin_popcount(unsigned(mask));
  if(label_count[mask]!=largest/binom[k])uniform_labels=false;
 }
 cout<<"all_outdegree9 "<<nine_edges<<" all_one_each_size "<<one_each
     <<" giant_each_size_permutation "<<permutation
     <<" giant_uniform_label_counts "<<uniform_labels<<'\n';
 assert(nine_edges&&one_each&&permutation&&uniform_labels);
 // The active-only weights omit both constant columns.  A size-k edge has
 // weight C(8,k)+2 for k=1,7, and C(8,k) otherwise.
 int active_weight[9]={0,10,28,56,70,56,28,10,0};
 int64_t active_length=0;array<int64_t,28> active_distance{};
 for(int v=0;v<n;v++)if(comp[v]==giant)
  for(size_t e=offsets[v];e<offsets[v+1];e++){
   auto x=edges[e];int w=active_weight[__builtin_popcount(unsigned(x.mask))];
   active_length+=w;
   for(int j=0;j<28;j++)if((x.inc>>j)&1)active_distance[j]+=w;
  }
 bool active_isotropic=active_length==int64_t(largest)*258;
 for(auto d:active_distance)active_isotropic&=d==int64_t(largest)*129;
 cout<<"active_only_isotropic "<<active_isotropic<<" length "<<active_length
     <<" pair_distance "<<active_distance[0]<<'\n';
 assert(active_isotropic);
 // For each directed edge u->v in the giant component, join a fixed path
 // root->u and a fixed path v->root.  Their 29-component label vectors span
 // the rational cycle-label space.  Rank 29 modulo a prime certifies full
 // rank over the rationals and avoids floating-point arithmetic.
 using Label=array<int32_t,29>;
 auto add_edge=[](Label &x,Edge const& e){
  x[0]++;for(int j=0;j<28;j++)x[j+1]+=(e.inc>>j)&1;
 };
 vector<Label> from(n),to(n);
 vector<uint8_t> seen_from(n,0),seen_to(n,0);
 vector<uint32_t> from_edge(n,UINT32_MAX),to_edge(n,UINT32_MAX);
 vector<uint32_t> queue;queue.reserve(largest);
 uint32_t root=uint32_t(terminal_id);
 queue.push_back(root);seen_from[root]=1;
 for(size_t at=0;at<queue.size();at++){
  auto u=queue[at];
  for(size_t e=offsets[u];e<offsets[u+1];e++){
   auto x=edges[e];if(comp[x.dest]!=giant||seen_from[x.dest])continue;
   seen_from[x.dest]=1;from[x.dest]=from[u];add_edge(from[x.dest],x);
   from_edge[x.dest]=uint32_t(e);
   queue.push_back(x.dest);
  }
 }
 assert(queue.size()==size_t(largest));
 vector<uint32_t> reverse_offsets(n+1,0);
 for(int u=0;u<n;u++)if(comp[u]==giant)
  for(size_t e=offsets[u];e<offsets[u+1];e++)reverse_offsets[edges[e].dest+1]++;
 for(int u=1;u<=n;u++)reverse_offsets[u]+=reverse_offsets[u-1];
 vector<uint32_t> reverse_edge(reverse_offsets[n]);
 auto cursor=reverse_offsets;
 for(int u=0;u<n;u++)if(comp[u]==giant)
  for(size_t e=offsets[u];e<offsets[u+1];e++)reverse_edge[cursor[edges[e].dest]++]=uint32_t(e);
 vector<uint32_t> edge_source(edges.size());
 for(int u=0;u<n;u++)for(size_t e=offsets[u];e<offsets[u+1];e++)edge_source[e]=u;
 queue.clear();queue.push_back(root);seen_to[root]=1;
 for(size_t at=0;at<queue.size();at++){
  auto u=queue[at];
  for(uint32_t z=reverse_offsets[u];z<reverse_offsets[u+1];z++){
   auto e=reverse_edge[z],v=edge_source[e];
   if(seen_to[v])continue;
   seen_to[v]=1;to[v]=to[u];add_edge(to[v],edges[e]);queue.push_back(v);
   to_edge[v]=e;
  }
 }
 assert(queue.size()==size_t(largest));
 constexpr int64_t prime=1000000007LL;
 array<array<int64_t,29>,29> basis{};array<int,29> has_pivot{};
 int rank=0;
 for(int u=0;u<n&&rank<29;u++)if(comp[u]==giant)
  for(size_t e=offsets[u];e<offsets[u+1]&&rank<29;e++){
   auto x=edges[e];if(x.mask==0||x.mask==255)continue;Label cyc{};
   for(int j=0;j<29;j++)cyc[j]=from[u][j]+to[x.dest][j];
   add_edge(cyc,x);
   array<int64_t,29> row{};
   for(int j=0;j<29;j++)row[j]=(cyc[j]%prime+prime)%prime;
   for(int j=0;j<29;j++)if(row[j]){
    if(has_pivot[j]){
     int64_t mult=row[j];
     for(int k=j;k<29;k++)row[k]=(row[k]-mult*basis[j][k]%prime+prime)%prime;
    }else{
     int64_t power=prime-2,a=row[j],inv=1;
     while(power){if(power&1)inv=inv*a%prime;a=a*a%prime;power>>=1;}
     for(int k=j;k<29;k++)basis[j][k]=row[k]*inv%prime;
     has_pivot[j]=1;rank++;
     vector<uint32_t> cycle_edges;
     for(uint32_t v=u;v!=root;v=edge_source[from_edge[v]]){
      assert(from_edge[v]!=UINT32_MAX);
      cycle_edges.push_back(from_edge[v]);
     }
     reverse(cycle_edges.begin(),cycle_edges.end());
     cycle_edges.push_back(uint32_t(e));
     for(uint32_t v=x.dest;v!=root;v=edges[to_edge[v]].dest){
      assert(to_edge[v]!=UINT32_MAX);
      cycle_edges.push_back(to_edge[v]);
     }
     assert(cycle_edges.size()==size_t(cyc[0]));
     Label replay{};uint32_t at=root;
     for(auto ce:cycle_edges){
      assert(edge_source[ce]==at&&comp[at]==giant);
      assert(edges[ce].mask!=0&&edges[ce].mask!=255);
      add_edge(replay,edges[ce]);at=edges[ce].dest;
     }
     assert(at==root&&replay==cyc);
     cout<<"cycle_witness "<<rank<<" edge "<<e<<" mask "<<int(x.mask)<<" labels";
     for(auto z:cyc)cout<<' '<<z;
     cout<<" path";for(auto ce:cycle_edges)cout<<' '<<int(edges[ce].mask);
     cout<<'\n';break;
    }
   }
  }
 cout<<"active_cycle_label_rank_mod_1000000007 "<<rank<<'\n';
 assert(rank==29);
}
