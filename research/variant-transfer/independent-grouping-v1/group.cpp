// Independently implemented metadata-only connected components and exposure scan.
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;
struct DNA {array<uint64_t,7> words{};bool operator==(const DNA&)const=default;bool operator<(const DNA&o)const{return words<o.words;}};
struct DH {size_t operator()(const DNA&d)const{uint64_t h=0x9e3779b97f4a7c15ULL;for(auto w:d.words)h^=w+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
int base(char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;default:throw runtime_error("nonliteral DNA");}}
DNA pack(const string&s){if(s.size()!=200)throw runtime_error("DNA length");DNA d;for(int i=0;i<200;++i)d.words[i/32]|=uint64_t(base(s[i]))<<(2*(i%32));return d;}
string rev(const string&s){string r;for(auto i=s.rbegin();i!=s.rend();++i)r+=string("TGCA")[base(*i)];return r;}
int dist(const DNA&a,const DNA&b,int bound){int n=0;for(int i=0;i<7;++i){uint64_t x=a.words[i]^b.words[i];n+=popcount((x|(x>>1))&0x5555555555555555ULL);if(n>bound)return n;}return n;}
vector<string> split(const string&s,char sep='\t'){vector<string> r;stringstream ss(s);string t;while(getline(ss,t,sep))r.push_back(t);return r;}
struct UF {vector<uint32_t> p,n;UF(size_t z):p(z),n(z,1){iota(p.begin(),p.end(),0);}uint32_t find(uint32_t x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}bool join(uint32_t a,uint32_t b){a=find(a);b=find(b);if(a==b)return false;if(n[a]<n[b])swap(a,b);p[b]=a;n[a]+=n[b];return true;}};
struct New {string id,chr,parent,ref,alt;uint32_t start,end;array<DNA,4> dna;};
struct Old {string id,s;DNA dna;vector<pair<int,string>> mutations;};
struct Guard {string collection,chr,id;uint32_t start,end;};
struct Entry {uint32_t key,group,target;bool operator<(const Entry&o)const{return tie(key,group,target)<tie(o.key,o.group,o.target);}};
uint32_t seed(const string&s,int j,int k){uint32_t w=0;for(int i=j*200/k;i<(j+1)*200/k;++i)w=4*w+base(s[i]);return (j<<20)|w;}
vector<Entry> index_entries(const vector<pair<string,uint32_t>>&seqs,const vector<uint32_t>&groups,int k){vector<Entry>x;x.reserve(seqs.size()*k);for(uint32_t i=0;i<seqs.size();++i)for(int j=0;j<k;++j)x.push_back({seed(seqs[i].first,j,k),groups.empty()?0:groups[seqs[i].second],i});sort(x.begin(),x.end());return x;}
pair<size_t,size_t> bounds(const vector<Entry>&x,uint32_t key){auto lo=lower_bound(x.begin(),x.end(),key,[](const Entry&a,uint32_t k){return a.key<k;});auto hi=upper_bound(lo,x.end(),key,[](uint32_t k,const Entry&a){return k<a.key;});return {lo-x.begin(),hi-x.begin()};}
struct Budget {uint64_t candidate=0,limit=1000000000ULL;void add(uint64_t n){candidate+=n;if(candidate>limit)throw runtime_error("independent candidate budget exceeded without dropping candidates");}};
int main(int argc,char**argv){try{
 if(argc!=7)throw runtime_error("usage: binary new.tsv old.tsv guards.tsv controls.tsv out.tsv events.tsv");
 vector<New> rows;string line;ifstream f(argv[1]);while(getline(f,line)){auto a=split(line);if(a.size()!=7)throw runtime_error("new field count");New r;r.id=a[0];r.chr=a[1];r.start=stoul(a[2]);r.end=stoul(a[3]);r.parent=a[4];r.ref=a[5];r.alt=a[6];if(r.end-r.start!=200)throw runtime_error("interval length");r.dna={pack(r.ref),pack(r.alt),pack(rev(r.ref)),pack(rev(r.alt))};rows.push_back(move(r));}if(rows.empty())throw runtime_error("empty new population");
 UF uf(rows.size());map<pair<string,uint32_t>,uint32_t> bins;map<string,uint32_t> parents;unordered_map<DNA,uint32_t,DH> exact;uint64_t physical=0,exactjoins=0;
 for(uint32_t i=0;i<rows.size();++i){auto key=make_pair(rows[i].chr,rows[i].start/1000000U);auto [p,added]=bins.emplace(key,i);if(!added)physical+=uf.join(i,p->second);auto [q,newparent]=parents.emplace(rows[i].parent,i);if(!newparent)physical+=uf.join(i,q->second);for(int a=0;a<2;++a){DNA canonical=min(rows[i].dna[a],rows[i].dna[a+2]);auto [t,newdna]=exact.emplace(canonical,i);if(!newdna)exactjoins+=uf.join(i,t->second);}}
 vector<uint32_t> coordinate_order(rows.size());iota(coordinate_order.begin(),coordinate_order.end(),0);sort(coordinate_order.begin(),coordinate_order.end(),[&](auto a,auto b){return tie(rows[a].chr,rows[a].start,rows[a].end)<tie(rows[b].chr,rows[b].start,rows[b].end);});string chr;uint32_t end=0,owner=0;for(auto i:coordinate_order){if(rows[i].chr!=chr){chr=rows[i].chr;end=rows[i].end;owner=i;}else{if(rows[i].start<end)physical+=uf.join(i,owner);if(rows[i].end>end){end=rows[i].end;owner=i;}}}
 vector<uint32_t> initial(rows.size());for(uint32_t i=0;i<rows.size();++i)initial[i]=uf.find(i);exact.clear();bins.clear();parents.clear();vector<pair<string,uint32_t>> strings;strings.reserve(rows.size()*2);for(uint32_t i=0;i<rows.size();++i){strings.emplace_back(rows[i].ref,i);strings.emplace_back(rows[i].alt,i);}auto idx=index_entries(strings,initial,21);vector<uint32_t> seen(rows.size(),UINT32_MAX);Budget budget;uint64_t checks=0,joins=0;ofstream events(argv[6]);events<<"stage\tnew_id\told_or_new_id\tdistance\n";
 for(uint32_t i=0;i<rows.size();++i){vector<uint32_t> peers;array<string,4> queries={rows[i].ref,rows[i].alt,rev(rows[i].ref),rev(rows[i].alt)};for(auto&s:queries)for(int j=0;j<21;++j){auto [lo,hi]=bounds(idx,seed(s,j,21));for(size_t p=lo;p<hi;){auto group=idx[p].group;auto stop=upper_bound(idx.begin()+p,idx.begin()+hi,group,[](uint32_t g,const Entry&e){return g<e.group;})-idx.begin();if(uf.find(group)!=uf.find(i))for(size_t z=p;z<size_t(stop);++z){uint32_t peer=strings[idx[z].target].second;if(peer<i&&seen[peer]!=i&&uf.find(peer)!=uf.find(i)){seen[peer]=i;peers.push_back(peer);}}p=stop;}}budget.add(peers.size());sort(peers.begin(),peers.end());for(auto p:peers){if(uf.find(i)==uf.find(p))continue;++checks;int best=201;for(int a=0;a<4;++a)for(int b=0;b<2;++b)best=min(best,dist(rows[i].dna[a],rows[p].dna[b],20));if(best<=20){uf.join(i,p);++joins;events<<"new_new\t"<<rows[i].id<<'\t'<<rows[p].id<<'\t'<<best<<'\n';}}if(i%10000==0)cerr<<"new_new "<<i<<" candidates "<<budget.candidate<<" joins "<<joins<<'\n';}
 vector<Entry>().swap(idx);vector<pair<string,uint32_t>>().swap(strings);vector<uint32_t>().swap(initial);
 vector<set<string>> reasons(rows.size());vector<Guard> guards;ifstream gf(argv[3]);while(getline(gf,line)){auto a=split(line);if(a.size()!=5)throw runtime_error("guard field count");guards.push_back({a[0],a[1],a[4],uint32_t(stoul(a[2])),uint32_t(stoul(a[3]))});}map<pair<string,string>,vector<Guard>> grouped;for(auto&g:guards)grouped[{g.collection,g.chr}].push_back(g);for(auto&[key,v]:grouped)sort(v.begin(),v.end(),[](auto&a,auto&b){return a.start<b.start;});
 map<pair<string,string>,vector<pair<int64_t,int64_t>>> expanded;for(auto&[key,v]:grouped){auto&spans=expanded[key];for(auto&g:v){int64_t left=int64_t(g.start)-10000,right=int64_t(g.end)+10000;if(!spans.empty()&&left<=spans.back().second)spans.back().second=max(spans.back().second,right);else spans.emplace_back(left,right);}}for(uint32_t i=0;i<rows.size();++i)for(auto&[key,v]:expanded){if(key.second!=rows[i].chr)continue;auto hi=upper_bound(v.begin(),v.end(),int64_t(rows[i].end),[](int64_t n,const pair<int64_t,int64_t>&g){return n<g.first;});if(hi!=v.begin()&&prev(hi)->second>=rows[i].start)reasons[uf.find(i)].insert("within_inclusive_10000bp_of_"+key.first);}
 set<string> controls;ifstream cf(argv[4]);while(getline(cf,line))if(!line.empty())controls.insert(line);vector<bool> control_component(rows.size(),false);for(uint32_t i=0;i<rows.size();++i)if(controls.count(rows[i].id))control_component[uf.find(i)]=true;
 vector<Old> explicitold,factor;ifstream of(argv[2]);while(getline(of,line)){auto a=split(line);if(a.size()<3||a.size()>4)throw runtime_error("old field count");Old o;o.id=a[1];o.s=a[2];o.dna=pack(o.s);if(a.size()==4&&!a[3].empty())for(auto&m:split(a[3],',')){auto t=split(m,':');o.mutations.emplace_back(stoi(t[0]),t[1]);}if(a[0]=="factor")factor.push_back(move(o));else if(a[0]=="explicit")explicitold.push_back(move(o));else throw runtime_error("old stage");}
 for(int stage=0;stage<2;++stage){auto&old=stage?factor:explicitold;int k=stage?22:21;vector<pair<string,uint32_t>> oldstrings;for(uint32_t p=0;p<old.size();++p)oldstrings.emplace_back(old[p].s,p);auto oi=index_entries(oldstrings,{},k);vector<uint32_t> marked(old.size(),UINT32_MAX);uint32_t epoch=0;uint64_t cleared=0;for(uint32_t i=0;i<rows.size();++i){auto root=uf.find(i);if(!reasons[root].empty())continue;bool found=false;array<string,4> qs={rows[i].ref,rows[i].alt,rev(rows[i].ref),rev(rows[i].alt)};for(auto&q:qs){++epoch;vector<uint32_t> peers;for(int j=0;j<k;++j){auto [lo,hi]=bounds(oi,seed(q,j,k));for(size_t p=lo;p<hi;++p){auto target=oi[p].target;if(marked[target]!=epoch){marked[target]=epoch;peers.push_back(target);}}}budget.add(peers.size());sort(peers.begin(),peers.end());auto packedq=pack(q);for(auto p:peers){int d=dist(packedq,old[p].dna,stage?21:20);bool match=d<=20;if(stage&&d==21)for(auto&[pos,alts]:old[p].mutations)if(q[pos]!=old[p].s[pos]&&alts.find(q[pos])!=string::npos){match=true;break;}if(match){reasons[root].insert("full200nt_ungapped_ge90pct_prior_DNA_match");events<<(stage?"old_factorized":"old_explicit")<<'\t'<<rows[i].id<<'\t'<<old[p].id<<'\t'<<(d==21?20:d)<<'\n';found=true;break;}}if(found)break;}if(!found)++cleared;}cerr<<"oldstage "<<stage<<" individually_cleared "<<cleared<<" cumulative_candidates "<<budget.candidate<<'\n';}
 ofstream out(argv[5]);out<<"id\troot_index\tretained_prior\treasons\ttechnical_control_component\n";for(uint32_t i=0;i<rows.size();++i){auto r=uf.find(i);out<<rows[i].id<<'\t'<<r<<'\t'<<int(reasons[r].empty())<<'\t';bool first=true;for(auto&s:reasons[r]){if(!first)out<<',';out<<s;first=false;}out<<'\t'<<int(control_component[r])<<'\n';}
 cerr<<"COMPLETE variants "<<rows.size()<<" candidate_total "<<budget.candidate<<" new_checks "<<checks<<" homology_joins "<<joins<<" physical_joins "<<physical<<" exact_joins "<<exactjoins<<'\n';
 }catch(const exception&e){cerr<<"FAILED "<<e.what()<<'\n';return 2;}}
