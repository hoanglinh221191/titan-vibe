#include "titan_topology_energy.h"
#include <array>
#include <iomanip>
#include <iostream>
int main(int argc,char **argv) {
  try {
    if(argc<2 || argc>3)return 2;
    cache_sigma_epsilon cache;cache.nonbonded_topology=titan_nb::read_topology(argv[1]);
    auto normal_top=std::make_shared<titan_nb::Topology>(*cache.nonbonded_topology);
    normal_top->listed_pairs.clear();
    cache_sigma_epsilon normal_cache;normal_cache.nonbonded_topology=normal_top;
    arguments args;args.force_f="AMBER99SB";args.epsilon_in=2;
    std::array<std::array<double,3>,5> xyz{{{1,1,1},{1.15,1,1},{1.30,1,1},{1.45,1,1},{1.1,1.15,1.02}}};
    if(argc==3)for(auto &point:xyz)for(double &x:point)x=static_cast<float>(x);
    std::array<std::array<double,3>,5> forces{};
    std::vector<mol_info> info;
    for(int i=0;i<5;++i){info.emplace_back(1,"C","C"+std::to_string(i+1),"A","T");info.back().set_original_atom_id(i);}
    titan_nb::validate_selection(*cache.nonbonded_topology,info);
    double elec=0,vdw=0,normal_elec=0,normal_vdw=0;
    for(int i=0;i<5;++i)for(int j=i+1;j<5;++j) {
      double d2=0;for(int k=0;k<3;++k)d2+=std::pow(xyz[i][k]-xyz[j][k],2);
      auto pair=titan_nb::evaluate_pair(cache,args,info[i],info[j],std::sqrt(d2),4);
      auto normal=titan_nb::evaluate_pair(normal_cache,args,info[i],info[j],std::sqrt(d2),4);
      elec+=pair.elec;vdw+=pair.vdw;
      normal_elec+=normal.elec;normal_vdw+=normal.vdw;
      for(int k=0;k<3;++k){double f=(pair.force_elec_prefactor+pair.force_vdw_prefactor)*(xyz[i][k]-xyz[j][k]);forces[i][k]+=f;forces[j][k]-=f;}
    }
    std::cout<<std::setprecision(17)<<normal_elec<<' '<<normal_vdw<<' '<<elec-normal_elec<<' '<<vdw-normal_vdw;
    for(const auto &f:forces)for(double v:f)std::cout<<' '<<v;
    std::cout<<'\n';return 0;
  }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
