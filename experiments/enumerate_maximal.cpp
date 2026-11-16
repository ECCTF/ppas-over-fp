// Complete (T^2+p)^2 Jacobian census by reduced binary-sextic models.
// Every degree d=5,6 model can be translated to kill x^(d-1), then
// scaled so that the next coefficient is 0, 1, or the least nonsquare.
// Both square classes of the leading coefficient are retained.
#define main legacy_census_main
#include "enumerate_curves.cpp"
#undef main
#include <unordered_set>
int main(int argc,char**argv){
 if(argc!=3)return 2;
 p=stoi(argv[1]);assert(p>=7&&p<100);pd[0]=1;
 for(int i=1;i<8;i++)pd[i]=pd[i-1]*p;
 invs.resize(p);chi.assign(p,-1);chi[0]=0;
 for(int i=1;i<p;i++){chi[i*i%p]=1;for(int j=1;j<p;j++)if(i*j%p==1)invs[i]=j;}
 nu=2;while(chi[nu]!=-1)nu++;
 auto ts=actions();unordered_set<long long> seen;
 string prefix=argv[2];ofstream out(prefix+"_maximal_curves.csv");
 out<<"class_id,c0,c1,c2,c3,c4,c5,c6,p_rank,a_number,aut_order,orbit_size,n1,n2,weil_a1,weil_a2\n";
 int cid=0,sp=0;long long checked=0;
 vector<pair<long long,int>> maximal_orbits;
 for(int d: {5,6})for(int a:{0,1,nu}){
  for(long long ix=0;ix<pd[d-2];ix++){
   Poly f{};f[d]=1;f[d-2]=a;long long v=ix;
   for(int i=0;i<d-2;i++){f[i]=v%p;v/=p;} checked++;
   auto m=cm(f);if(rank2(m)||!smooth(f))continue;
   for(int scalar: {1,nu}){
    Poly g=f;for(int&i:g)i=mod(i*scalar);
    if(seen.count(encode(g)))continue;
    long long minimum=encode(g);int stab=0;
    for(auto&t:ts){long long q=encode(transform_form(t,g));seen.insert(q);minimum=min(minimum,q);if(q==encode(g))stab++;}
    g=decode(minimum);sp++;auto[n1,n2]=pointcounts(g);
    if(n1!=p+1||n2!=p*p+1+4*p)continue;
    for(auto&t:ts)maximal_orbits.emplace_back(encode(transform_form(t,g)),cid);
    out<<cid++;for(int c:g)out<<','<<c;
    out<<",0,2,"<<2*stab<<','<<ts.size()/stab<<','<<n1<<','<<n2<<",0,"<<2*p<<'\n';
   }
  }
  cerr<<"p="<<p<<" degree="<<d<<" normalized_coefficient="<<a<<" checked="<<checked<<" maximal="<<cid<<endl;
 }
 sort(maximal_orbits.begin(),maximal_orbits.end());
 maximal_orbits.erase(unique(maximal_orbits.begin(),maximal_orbits.end()),maximal_orbits.end());
 ofstream orbitout(prefix+"_maximal_orbits.csv");
 for(auto [code,i]:maximal_orbits)orbitout<<code<<','<<i<<'\n';
 cout<<"p="<<p<<" maximal_Jacobians="<<cid<<" superspecial_all_blocks="<<sp<<" normalized_candidates="<<checked<<endl;
}
