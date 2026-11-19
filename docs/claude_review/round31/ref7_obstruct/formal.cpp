// Formal-X brute force for VII.2 (d), (e) and the translate/fibre class F1.
// K = F_{q^{2r}}, conj x -> x^{q^r}; V = i*span_Fq(1,theta,..,theta^{k-1}); tau = i*theta^k.
#include "gf.h"
template<class F> struct Ctx {
  F f; typedef typename F::E E; int q,r,k; ll n,h; E I,theta; vector<E> V; vector<vector<int>> X; E tau;
  E conj(E a){ ll e=1; for(int j=0;j<r;j++) e*=q; return f.pw_(a,e); }
  int chi(int y){ y%=q; if(y<0) y+=q; if(!y) return 0; ll e=1; for(int j=0;j<(q-1)/2;j++) e=e*y%q; return e==1?1:-1; }
  void setup(int q_,int r_,int k_, vector<int> m){ q=q_;r=r_;k=k_; f.init(q,m); ll qr=1; for(int j=0;j<r;j++) qr*=q;
    E g=f.gen(); I=f.pw_(g,(qr+1)/2); theta=f.pw_(g,qr+1);
    // checks
    if(conj(I)!=f.neg(I)){fprintf(stderr,"i not anti-invariant\n");exit(1);} 
    n=1; for(int j=0;j<k;j++) n*=q; h=n/q;
    vector<E> basis; E t=I; for(int b=0;b<k;b++){ basis.push_back(t); t=f.mul(t,theta);} tau=t;
    for(ll idx=0;idx<n;idx++){ vector<int> x(k); ll s=idx; E v=f.zero(); for(int b=0;b<k;b++){x[b]=s%q; s/=q; v=f.add(v,f.mul(f.fromint(x[b]),basis[b]));} V.push_back(v); X.push_back(x);} }
  int form(const vector<int>&l, ll idx){ int s=0; for(int b=0;b<k;b++) s+=l[b]*X[idx][b]; return ((s%q)+q)%q; }
  ll negidx(ll idx){ ll r_=0,mul=1; for(int b=0;b<k;b++){ r_+=((q-X[idx][b])%q)*mul; mul*=q;} return r_; }
};
template<class F> Poly<F> prodroots(const F&f, const vector<typename F::E>&rts){ Poly<F> a={f.one()}; for(auto c:rts) mullin(f,a,c); return a; }
template<class F> int excess(const F&f, const Poly<F>&A,const Poly<F>&B, int D, int g){ Poly<F> d=psub(f,A,B); if(d.empty()) return INT_MIN; return 2*deg(f,d)-(D+g); }
template<class F> void run(Ctx<F>&C, int mode, ll maxpairs, int seed){
  typedef typename F::E E; auto&f=C.f; int q=C.q; ll n=C.n,h=C.h; mt19937 rng(seed);
  // forms: all nonproportional pairs (normalised first nonzero coefficient 1) -- take a few
  vector<vector<int>> forms; { ll tot=n; for(ll idx=1;idx<tot;idx++){ vector<int> c=C.X[idx]; int fn=-1; for(int b=0;b<C.k;b++) if(c[b]){fn=b;break;} if(c[fn]==1) forms.push_back(c);} }
  int best=INT_MAX; ll ncrit=0, npairs=0;
  int npf=0;
  for(size_t a=0;a<forms.size();a++) for(size_t b=a+1;b<forms.size();b++){ if(npf++>=6) break; // a few form pairs
   for(int si=-1;si<=1;si+=2) for(int sj=-1;sj<=1;sj+=2){ if(mode==2 && (si!=1||sj!=1)) continue;
    auto build=[&](const vector<int>&l,int eps, ll fill, vector<ll>&rootids)->Poly<F>{
      vector<E> rts; rootids.clear();
      if(mode==0){ // (d): roots V\0, positive char off ker, kernel pairs filled by bits of fill
        int bit=0; vector<char> done(n,0);
        for(ll idx=1;idx<n;idx++){ int c=C.chi(eps*C.form(l,idx)); if(c==1){rts.push_back(C.V[idx]); rootids.push_back(idx);} 
          else if(c==0 && !done[idx]){ ll ni=C.negidx(idx); done[idx]=done[ni]=1; ll pick=((fill>>bit)&1)?idx:ni; bit++; rts.push_back(C.V[pick]); rootids.push_back(pick);} }
      } else if(mode==2){ // F1: whole fibres, arbitrary S subset F_q encoded by bits of fill (eps ignored)
        for(ll idx=0;idx<n;idx++){ int lam=C.form(l,idx); if((fill>>lam)&1){ rts.push_back(f.add(C.V[idx],C.tau)); rootids.push_back(idx);} }
        for(ll idx=0;idx<n;idx++){ int lam=C.form(l,idx); int ml=(q-lam)%q; if(!((fill>>ml)&1)){ rts.push_back(f.sub(C.V[idx],C.tau)); rootids.push_back(n+idx);} }
      } else { // (e): roots V+tau (id idx) and V-tau (id n+idx); positive labels both; kernel pairs {w+tau, -w-tau}
        int bit=0;
        for(ll idx=0;idx<n;idx++){ int c=C.chi(eps*C.form(l,idx)); if(c==1){ rts.push_back(f.add(C.V[idx],C.tau)); rootids.push_back(idx); rts.push_back(f.sub(C.V[idx],C.tau)); rootids.push_back(n+idx);} 
          else if(c==0){ ll ni=C.negidx(idx); // pair (idx,+) with (ni,-)
            if(((fill>>bit)&1)){ rts.push_back(f.add(C.V[idx],C.tau)); rootids.push_back(idx);} else { rts.push_back(f.sub(C.V[ni],C.tau)); rootids.push_back(n+ni);} bit++; } }
      }
      return prodroots(f,rts); };
    ll nfill = (mode==0)? (1LL<<((h-1)/2)) : (mode==2? (1LL<<q) : (1LL<<h));
    // build all fills (or sample)
    ll cap=min(nfill,(ll)2048); vector<ll> fills; if(nfill<=cap){ for(ll x=0;x<nfill;x++) fills.push_back(x);} else { for(ll x=0;x<cap;x++) fills.push_back(((ll)rng()<<20 ^ rng())%nfill);} 
    vector<Poly<F>> Pa,Pb; vector<vector<ll>> Ra,Rb; for(ll x:fills){ vector<ll> ids; Pa.push_back(build(forms[a],si,x,ids)); sort(ids.begin(),ids.end()); Ra.push_back(ids); Pb.push_back(build(forms[b],sj,x,ids)); sort(ids.begin(),ids.end()); Rb.push_back(ids);} 
    int D=deg(f,Pa[0]);
    ll lim=maxpairs; 
    for(size_t x=0;x<fills.size();x++) for(size_t y=0;y<fills.size();y++){ if((ll)(x*fills.size()+y)>=lim) break; vector<ll> inter; set_intersection(Ra[x].begin(),Ra[x].end(),Rb[y].begin(),Rb[y].end(),back_inserter(inter)); int g=inter.size(); int ex=excess(f,Pa[x],Pb[y],D,g); npairs++; if(ex==INT_MIN) continue; if(ex<=0) ncrit++; best=min(best,ex);} 
    printf("forms (%d,%d)x(%d,%d) signs %d %d: D=%d fills=%zu best excess=%d (critical iff <=0)\n",forms[a][0],forms[a][1],forms[b][0],forms[b][1],si,sj,D,fills.size(),best);
   }
  }
  printf("mode %d q=%d k=%d n=%lld: pairs=%lld critical=%lld min excess=%d\n",mode,q,C.k,n,npairs,ncrit,best);
}
int main(int argc,char**argv){ int q=atoi(argv[1]), r=atoi(argv[2]), k=atoi(argv[3]), mode=atoi(argv[4]); ll maxpairs=atoll(argv[5]); bool dir=atoi(argv[6]); vector<int> m; for(int a=7;a<argc;a++) m.push_back(atoi(argv[a]));
  if(dir){ Ctx<GFdir> C; C.setup(q,r,k,m); run(C,mode,maxpairs,1);} else { Ctx<GFlog> C; C.setup(q,r,k,m); run(C,mode,maxpairs,1);} }
