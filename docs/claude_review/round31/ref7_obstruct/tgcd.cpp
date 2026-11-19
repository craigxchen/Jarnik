#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
int main(){ F f; f.init(3,{2,0,2,0,1,1,1}); mt19937 rng(5);
 for(int it=0;it<2000;it++){ int da=rng()%6, db=rng()%6, dc=rng()%6; auto rp=[&](int d){P a(d+1); for(auto&x:a) x=(rng()%5==0)?-1:(E)(rng()%(f.Q-1)); a[d]=0; return a;};
  P a=rp(da),b=rp(db),c=rp(dc); P A=pmul(f,a,c),B=pmul(f,b,c); P g=pgcd(f,A,B); P q,r; pdivmod(f,A,g,q,r); if(!r.empty()){printf("fail it %d da %d db %d dc %d degg %d\n",it,da,db,dc,deg(f,g)); return 0;}
  P q2,r2; pdivmod(f,g,c,q2,r2); if(!r2.empty() ) {printf("c does not divide g %d\n",it);return 0;} 
  // check (x-1)^2 derivative
 }
 printf("gcd ok\n"); }
