#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <string>
using namespace std;
bool allowed(int w,int M){
 for(int a=0;a<M;a++)for(int b=a+1;b<M;b++)for(int c=b+1;c<M;c++)for(int d=c+1;d<M;d++)for(int e=d+1;e<M;e++){
   int x=(w>>a)&1,y=(w>>b)&1,z=(w>>c)&1,t=(w>>d)&1,u=(w>>e)&1;
   if(x==u&&y==z&&z==t&&x!=y)return false;
 }return true;
}
bool search(const vector<int>&g,int M,int at,int left,array<int,4>&sel,vector<int>&ans){
 if(left==0){
  int d[200]={};
  for(int q=0;q<4;q++) for(int i=0;i<M;i++) for(int j=0;j<i;j++) if(((sel[q]>>i)&1)!=((sel[q]>>j)&1)) d[i* M+j]++;
  for(int i=0;i<M;i++)for(int j=0;j<i;j++)if(d[i*M+j]<2)return false;
  ans.assign(sel.begin(),sel.end());return true;
 }
 // Keep at=q on recursion: this enumerates multisets (repeated columns)
 // rather than only four-element subsets.
 for(int q=at;q<(int)g.size();q++){
  sel[4-left]=g[q];if(search(g,M,q,left-1,sel,ans))return true;
 }
 return false;
}
int main(){for(int M=5;M<=11;M++){vector<int>g;for(int w=0;w<(1<<M);w++)if(allowed(w,M))g.push_back(w);array<int,4>s;vector<int>a;bool ok=search(g,M,0,4,s,a);cout<<"M="<<M<<" words="<<g.size()<<" k4_example="<<ok; if(ok){cout<<" cols=";for(int w:a)cout<<" "<<string([&](){string x;for(int i=0;i<M;i++)x+=char('0'+((w>>i)&1));return x;}());}cout<<"\n";}}
