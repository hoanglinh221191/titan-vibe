#include "titan_gb_model.h"
#include <array>
#include <iostream>
#ifdef TITAN_TEST_ENGINE
#include "gbsa_intrin.h"
#endif
int main(){try {
 std::array<double,2> x{0,100},y{},z{},r{2,0};
 auto size=[&]{return titan_gb::electrostatic_size_adet_nm(x,y,z,r);};
 if(!std::isfinite(size()) || std::abs(size()-.2)>1e-12)throw std::runtime_error("zero-radius spectator changed geometric size");
 x[1]=-1e6;if(!std::isfinite(size()) || std::abs(size()-.2)>1e-12)throw std::runtime_error("zero-volume site position affected size");
 r={0,0};if(std::isfinite(size()))throw std::runtime_error("zero total volume accepted");
 r={2,-1};if(std::isfinite(size()))throw std::runtime_error("negative radius accepted");
#ifdef TITAN_TEST_ENGINE
 MOL rec,lig;rec.pos.resize(1);lig.pos.resize(1);Multi_MOL_parameters rp(1),lp(1);std::vector<mol_info> ri{mol_info(1,"H","H","A","NUM")},li=ri;
 GBSASasaCache cache;cache.sasa_rec_free={0};cache.sasa_lig_free={0};cache.sasa_complex={0,0};arguments args;args.flag_alpb=true;args.int_diel=1;args.ext_diel=78.5;args.igb=5;
 cache_sigma_epsilon paircache;std::vector<double> unused;AlignedDoubleVector born(1,.1);GBSA_energy_snap energy;bool rejected=false;
 try{cal_energy_gbsa_avx2(rec,lig,ri,li,paircache,0,1,0,1,energy,0,args,rp,lp,unused,unused,unused,unused,born,born,born,born,nullptr,&cache,false);}catch(const std::invalid_argument&){rejected=true;}
 if(!rejected)throw std::runtime_error("invalid ALPB cache did not raise explicit error");
#endif
 std::cout<<"PASS ALPB zero-radius sites, zero-volume/negative-size rejection";
#ifdef TITAN_TEST_ENGINE
 std::cout<<", production invalid-cache error";
#endif
 std::cout<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
