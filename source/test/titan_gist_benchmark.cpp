// Comparable workloads: exact neighbor search vs exhaustive pairs, exact SO3
// KD search vs exhaustive quaternion dot products. No old/new physics mixing.
#include "titan_gist.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
using Clock=std::chrono::steady_clock;
template<class F> double median_seconds(F f){std::vector<double> samples;for(int i=0;i<3;++i){auto start=Clock::now();f();samples.push_back(std::chrono::duration<double>(Clock::now()-start).count());}std::sort(samples.begin(),samples.end());return samples[1];}
int main(){try{
 std::cout<<std::setprecision(12)<<"kind,n,fast_median_s,exhaustive_median_s,speedup,max_absolute_error\n";
 for(int n:{1024,4096,16384}){std::mt19937_64 rng(981);std::normal_distribution<double> normal;std::vector<titan_gist::Quaternion> q(n);for(auto& p:q){double norm=0;for(auto& a:p){a=normal(rng);norm+=a*a;}for(auto& a:p)a/=std::sqrt(norm);}std::vector<double> fast(n),slow(n);
  auto kd=[&]{titan_gist::QuaternionTree tree(q);for(int i=0;i<n;++i)fast[i]=tree.nearest_angle(i);};auto brute=[&]{for(int i=0;i<n;++i){double best=0;for(int j=0;j<n;++j)if(i!=j){double dot=0;for(int k=0;k<4;++k)dot+=q[i][k]*q[j][k];best=std::max(best,std::abs(dot));}slow[i]=2*std::acos(std::min(1.,best));}};
  kd();brute();double error=0;for(int i=0;i<n;++i)error=std::max(error,std::abs(fast[i]-slow[i]));if(error>1e-9)throw std::runtime_error("Quaternion benchmark correctness failed");double tf=median_seconds(kd),ts=median_seconds(brute);std::cout<<"orientation,"<<n<<','<<tf<<','<<ts<<','<<ts/tf<<','<<error<<std::endl;
 }
 std::vector<mol_info> ri,li;Multi_MOL_parameters rp(0),lp(1);MOL lig;lig.pos.resize(1);
 for(int side:{6,10,14}){int n=side*side*side;TitanGIST g;g.init(lig,1,.5,300,9,1,"TIP3P",rp,ri,lp,li);g.periodic_box_.emplace(box(side*3.2,side*3.2,side*3.2,90.,90.,90.));g.wat_buf_.resize(3*n);std::vector<int> voxels(n,0);
  for(int i=0;i<n;++i){double x=(i%side)*3.2,y=((i/side)%side)*3.2,z=(i/(side*side))*3.2;g.wat_buf_[3*i]={x,y,z};g.wat_buf_[3*i+1]={x+.7,y+.6,z};g.wat_buf_[3*i+2]={x-.7,y+.6,z};}
  auto fast=g.water_water_energies(voxels),slow=g.water_water_energies(voxels,false);double error=0;for(int i=0;i<n;++i)error=std::max(error,std::abs(fast[i]-slow[i]));if(error>1e-8)throw std::runtime_error("Eww benchmark correctness failed");double tf=median_seconds([&]{fast=g.water_water_energies(voxels);}),ts=median_seconds([&]{slow=g.water_water_energies(voxels,false);});std::cout<<"Eww,"<<n<<','<<tf<<','<<ts<<','<<ts/tf<<','<<error<<std::endl;
 }
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
