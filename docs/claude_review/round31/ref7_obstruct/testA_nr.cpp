// VII.2(a): specialisations X=T(t) of the VI.4 marks on V = i F_q[t]_{<k}; is the reduced common norm ever squarefree?
#include "gf.h"
typedef GFlog F; typedef F::E E; typedef Poly<F> P;
F f; int p,r,k; ll q; E I; vector<E> Fq; 
int chi(E y){ if(f.isz(y)) return 0; ll e=(ll)y*((q-1)/2)%(f.Q-1); return e==0?1:-1; }
E conjE(E a){ return a<0?a:(E)(((ll)a*q)%(f.Q-1)); }
P conjP(const P&a){P b(a.size());for(size_t i=0;i<a.size();i++) b[i]=conjE(a[i]);return b;}
P powP(P a,int e){P r={f.one()};while(e--) r=pmul(f,r,a);return r;}
int main(int argc,char**argv){
  // args: p r k seed ntests maxdegT  [modpoly coeffs...]
  p=atoi(argv[1]); r=atoi(argv[2]); k=atoi(argv[3]); int seed=atoi(argv[4]); int ntests=atoi(argv[5]); int maxdT=atoi(argv[6]);
  vector<int> m; for(int a=7;a<argc;a++) m.push_back(atoi(argv[a]));
  f.init(p,m); q=1; for(int a=0;a<r;a++) q*=p; if(f.Q!=q*q){fprintf(stderr,"bad field\n");return 1;}
  I=(E)((f.Q-1)/4); // i^2=-1
  Fq.push_back(f.zero()); for(ll a=0;a<f.Q-1;a+=q+1) Fq.push_back((E)a);
  // forms over F_p up to sign: first nonzero coeff =1
  vector<vector<int>> forms; { vector<int> c(k,0); ll tot=1; for(int a=0;a<k;a++) tot*=p; for(ll idx=1;idx<tot;idx++){ ll t=idx; for(int a=0;a<k;a++){c[a]=t%p;t/=p;} int fn=-1; for(int a=0;a<k;a++) if(c[a]){fn=a;break;} if(c[fn]==1) forms.push_back(c);} }
  int nf=forms.size(); fprintf(stderr,"q=%lld k=%d forms=%d\n",q,k,nf);
  ll n=1; for(int a=0;a<k;a++) n*=q;
  // labels
  vector<vector<E>> X(n, vector<E>(k)); for(ll idx=0;idx<n;idx++){ ll t=idx; for(int a=0;a<k;a++){X[idx][a]=Fq[t%q]; t/=q;} }
  mt19937 rng(seed);
  int nsf=0, ntot=0;
  for(int test=0;test<ntests;test++){
    // choose m forms (random subset size m in [2,nf]) and signs, and T in F_p[t] nonzero random degree <= maxdT
    int mm = 2 + rng()%(nf-1); vector<int> perm(nf); iota(perm.begin(),perm.end(),0); shuffle(perm.begin(),perm.end(),rng); perm.resize(mm);
    vector<int> sg(mm); for(auto&s:sg) s=(rng()&1)?1:-1;
    int dT = rng()%(maxdT+1); P T(dT+1); bool nz=false; for(int a=0;a<=dT;a++){ T[a]=f.fromint(rng()%p); if(!f.isz(T[a])) nz=true;} if(!nz) T[0]=f.one(); trim(f,T);
    if(test==0){ T={f.zero(),f.one()}; } // T=t
    if(getenv("NONREAL")){ int d2=rng()%(maxdT+1); P T2(d2+1); for(int a=0;a<=d2;a++) T2[a]=f.mul(I,f.fromint(rng()%p)); trim(f,T2); if(T2.empty()) T2={I}; T.resize(max(T.size(),T2.size()),f.zero()); for(size_t a=0;a<T2.size();a++) T[a]=f.add(T[a],T2[a]); trim(f,T);}
    // exponents
    vector<vector<int>> e(mm, vector<int>(n)); vector<int> mn(n,2); ll mixed=0;
    for(ll idx=0;idx<n;idx++){ bool pos=false,ng=false; for(int j=0;j<mm;j++){ E val=f.zero(); for(int a=0;a<k;a++) if(forms[perm[j]][a]) val=f.add(val,f.mul(f.fromint(forms[perm[j]][a]),X[idx][a])); int c=chi(val)*sg[j]; e[j][idx]=1+c; mn[idx]=min(mn[idx],1+c); if(c>0)pos=true; if(c<0)ng=true;} if(pos&&ng) mixed++; }
    // G_formal and residuals
    P Gf={f.one()}; vector<P> R(mm,P{f.one()}); P PV={f.one()};
    for(ll idx=0;idx<n;idx++){ P lin=T; lin.resize(max((int)lin.size(),k),f.zero()); for(int a=0;a<k;a++) lin[a]=f.sub(lin[a],f.mul(I,X[idx][a])); trim(f,lin);
      PV=pmul(f,PV,lin); for(int c=0;c<mn[idx];c++) Gf=pmul(f,Gf,lin); for(int j=0;j<mm;j++) for(int c=0;c<e[j][idx]-mn[idx];c++) R[j]=pmul(f,R[j],lin); }
    P g=R[0]; for(int j=1;j<mm;j++) g=pgcd(f,g,R[j]);
    P G=pmul(f,Gf,g);
    // norm check equal for all j: Z_j = Gf*R_j/G ; N = Z conj(Z)
    P qq,rr; P Z0; pdivmod(f,R[0],g,Z0,rr); if(!rr.empty()){fprintf(stderr,"div err\n");return 1;}
    P N0=pmul(f,Z0,conjP(Z0)); 
    bool normok=true; for(int j=1;j<mm && j<3;j++) if(!getenv("NONREAL")){ P Zj; pdivmod(f,R[j],g,Zj,rr); P Nj=pmul(f,Zj,conjP(Zj)); if(psub(f,Nj,N0).size()) normok=false; }
    // full norm = (PV^2/(G conjG)) ; N0 above is the residual part times...; compute full N
    if(getenv("NONREAL")){ P Fa=pmul(f,Gf,R[0]), Fb=pmul(f,Gf,R[1]); P Na=pmul(f,Fa,conjP(Fa)), Nb=pmul(f,Fb,conjP(Fb)); printf("nonreal test %d m=%d: raw norms equal=%d\n",test,mm,(int)psub(f,Na,Nb).empty()); continue; }
    P PV2=pmul(f,PV,PV); P GG=pmul(f,G,conjP(G)); P N; pdivmod(f,PV2,GG,N,rr); if(!rr.empty()){fprintf(stderr,"div2 err\n");return 1;}
    P dN=pderiv(f,N); P h=pgcd(f,N,dN); bool sf = deg(f,h)==0;
    // also is N in F_p[t]?
    bool base=true; for(auto c:N){ if(!f.isz(c)){ ll d=f.todig(c); if(d>=p) base=false; } }
    ntot++; if(sf) nsf++;
    printf("test %d m=%d degT=%d mixed=%lld/%lld degG=%d degN=%d sqfree=%d normsequal=%d NinFp=%d degsqpart=%d\n",test,mm,deg(f,T),mixed,n,deg(f,G),deg(f,N),sf,normok,base,deg(f,h));
  }
  printf("SUMMARY q=%lld k=%d tests=%d squarefree=%d\n",q,k,ntot,nsf);
}
