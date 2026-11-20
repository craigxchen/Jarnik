// Exhaustive: k pairwise non-parallel vectors from [-R,R]^2 \ {0}; maximise the number of
// zero-sum subsets (empty set included).  r = 2, so E = k/(k-3).
#include <stdio.h>
#include <stdlib.h>
int R,k,np; int px[400],py[400]; int sel[16]; int best=0; int bsel[16];
long long cnt_calls=0;
static int par(int a,int b){ return px[a]*py[b]-py[a]*px[b]==0; }
static void eval(){
  int c=0; int n=1<<k;
  for(int m=0;m<n;m++){ int sx=0,sy=0; for(int j=0;j<k;j++) if(m>>j&1){sx+=px[sel[j]];sy+=py[sel[j]];} if(!sx&&!sy) c++; }
  if(c>best){best=c; for(int j=0;j<k;j++) bsel[j]=sel[j];
    printf("best=%d :",best); for(int j=0;j<k;j++) printf(" (%d,%d)",px[bsel[j]],py[bsel[j]]); printf("\n"); fflush(stdout);}
}
static void rec(int pos,int start){
  if(pos==k){eval();return;}
  for(int a=start;a<np;a++){ int ok=1; for(int j=0;j<pos;j++) if(par(a,sel[j])){ok=0;break;} if(!ok) continue; sel[pos]=a; rec(pos+1,a+1);} }
int main(int argc,char**argv){ R=atoi(argv[1]); k=atoi(argv[2]); np=0;
  // canonical: first vector restricted to x>0 or (x==0,y>0) handled by symmetry z->-z of whole set? not used.
  for(int x=-R;x<=R;x++) for(int y=-R;y<=R;y++) if(x||y){px[np]=x;py[np]=y;np++;}
  rec(0,0); printf("final best=%d\n",best); return 0; }
