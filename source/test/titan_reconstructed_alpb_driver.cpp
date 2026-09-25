#include "gbsa_intrin.h"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){try {
 if(argc!=2)return 2;std::ifstream input(argv[1]);int nr,nl,igb;double declared_salt;input>>nr>>nl>>igb>>declared_salt;
 if(!input||nr<=0||nl<=0)return 2;
 MOL rec,lig;rec.pos.resize(nr);lig.pos.resize(nl);Multi_MOL_parameters rp(nr),lp(nl);std::vector<mol_info> ri,li;
 auto load=[&](MOL& mol,Multi_MOL_parameters& params,std::vector<mol_info>& info){for(std::size_t i=0;i<mol.size();++i){std::string name,type;double x,y,z,q,s,e,r,gb,screen;if(!(input>>name>>type>>x>>y>>z>>q>>s>>e>>r>>gb>>screen))throw std::runtime_error("Invalid reconstructed fixture");mol.pos.x[i]=x;mol.pos.y[i]=y;mol.pos.z[i]=z;params[i].set_data(q,r,s,e,12.);info.emplace_back(1,type,name,"A","TST");}};
 load(rec,rp,ri);load(lig,lp,li);AlignedDoubleVector rf,lf,rc,lc;
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rf,igb,0,nr);cal_effective_born_radii_avx2(lig,lp.radius,li,lf,igb,0,nl);
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rc,igb,0,nr,&lig,&lp.radius,&li);cal_effective_born_radii_avx2(lig,lp.radius,li,lc,igb,0,nl,&rec,&rp.radius,&ri);
 GBSASasaCache cache_sasa;
 // Exercise production A_det cache construction. The low-resolution SASA
 // values are not reported or used for the ALPB polar comparison.
 build_gbsa_sasa_cache(rec,lig,rp,lp,cache_sasa,32,true,nullptr);
 arguments args;args.igb=igb;args.int_diel=1;args.ext_diel=78.5;args.force_f="AMBER99SB";args.flag_alpb=true;args.flag_use_common_topo=true;
 cache_sigma_epsilon cache;std::vector<double> unused;std::cout<<std::setprecision(17);
 for(double salt:{0.,.15}){args.saltcon=salt;GBSA_energy_snap e;if(cal_energy_gbsa_avx2(rec,lig,ri,li,cache,0,nr,0,nl,e,0,args,rp,lp,unused,unused,unused,unused,rf,lf,rc,lc,nullptr,&cache_sasa,false)!=0)throw std::runtime_error("ALPB kernel rejected fixture");
  std::cout<<salt<<' '<<e.E_elec<<' '<<e.E_vdw<<' '<<e.G_polar<<' '<<10*cache_sasa.alpb_rec_electrostatic_size_nm<<' '<<10*cache_sasa.alpb_lig_electrostatic_size_nm<<' '<<10*cache_sasa.alpb_complex_electrostatic_size_nm<<'\n';}
 return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
