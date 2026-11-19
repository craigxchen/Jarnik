// read: K k, then k lines of K digits (F_9 as a+3b) for phi and K digits for phi^{-1}; output max class size and #classes
#include <bits/stdc++.h>
using namespace std; int MUL[9][9],ADD[9][9]; int K;
typedef array<uint8_t,32> S;
struct H{ size_t operator()(const S&s)const{ size_t h=0; for(int i=0;i<K;i++) h=h*1000003u+s[i]; return h;} };
S mul(const S&a,const S&b){ S r{}; for(int i=0;i<K;i++) if(a[i]) for(int j=0;i+j<K;j++) r[i+j]=ADD[r[i+j]][MUL[a[i]][b[j]]]; return r; }
int main(){ for(int x=0;x<9;x++)for(int y=0;y<9;y++){int a=x%3,b=x/3,c=y%3,d=y/3;MUL[x][y]=((a*c-b*d)%3+3)%3+3*((a*d+b*c)%3);ADD[x][y]=(a+c)%3+3*((b+d)%3);}
  int k; cin>>K>>k; unordered_map<S,long long,H> m; S one{}; one[0]=1; m[one]=1;
  for(int p=0;p<k;p++){ S f{},g{}; for(int i=0;i<K;i++){int x;cin>>x;f[i]=x;} for(int i=0;i<K;i++){int x;cin>>x;g[i]=x;}
    unordered_map<S,long long,H> n; n.reserve(m.size()*2+16); for(auto&e:m){ n[mul(e.first,f)]+=e.second; n[mul(e.first,g)]+=e.second; } m.swap(n); }
  long long mx=0,c=0; for(auto&e:m){ if(e.second>mx){mx=e.second;c=1;} else if(e.second==mx) c++; } printf("%lld %lld %zu\n",mx,c,m.size()); }
