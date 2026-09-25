#include "gbsa_intrin.h"
#include "titan_sasa.h"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){try{
 if(argc!=3)return 2;std::ifstream input(argv[1]);int nr,nl,igb;double declared_salt;input>>nr>>nl>>igb>>declared_salt;
 if(!input||nr<=0||nl<=0)return 2;
 MOL rec,lig;rec.pos.resize(nr);lig.pos.resize(nl);Multi_MOL_parameters rp(nr),lp(nl);std::vector<mol_info> ri,li;
 auto load=[&](MOL& mol,Multi_MOL_parameters& params,std::vector<mol_info>& info){for(std::size_t i=0;i<mol.size();++i){std::string name,type;double x,y,z,q,s,e,r,gb,screen;if(!(input>>name>>type>>x>>y>>z>>q>>s>>e>>r>>gb>>screen))throw std::runtime_error("Invalid reconstructed fixture");mol.pos.x[i]=x;mol.pos.y[i]=y;mol.pos.z[i]=z;params[i].set_data(q,r,s,e,12.);info.emplace_back(1,type,name,"A","TST");}};
 load(rec,rp,ri);load(lig,lp,li);AlignedDoubleVector rf,lf,rc,lc;
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rf,igb,0,nr);
 cal_effective_born_radii_avx2(lig,lp.radius,li,lf,igb,0,nl);
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rc,igb,0,nr,&lig,&lp.radius,&li);
 cal_effective_born_radii_avx2(lig,lp.radius,li,lc,igb,0,nl,&rec,&rp.radius,&ri);
 GBSASasaCache sasa;build_gbsa_sasa_cache(rec,lig,rp,lp,sasa,std::stoi(argv[2]),false,nullptr);
 double ar=TitanSASA::sum_sasa(sasa.sasa_rec_free),al=TitanSASA::sum_sasa(sasa.sasa_lig_free),ac=TitanSASA::sum_sasa(sasa.sasa_complex);
 arguments args;args.igb=igb;args.int_diel=1;args.ext_diel=78.5;args.force_f="AMBER99SB";args.flag_alpb=false;args.flag_use_common_topo=true;
 cache_sigma_epsilon cache;std::vector<double> unused;
 std::cout<<std::setprecision(17);
 for(double salt:{0.,.15}) {args.saltcon=salt;GBSA_energy_snap energy;
  cal_energy_gbsa_avx2(rec,lig,ri,li,cache,0,nr,0,nl,energy,0,args,rp,lp,unused,unused,unused,unused,rf,lf,rc,lc,nullptr,&sasa,false);
  std::cout<<salt<<' '<<energy.E_elec<<' '<<energy.E_vdw<<' '<<energy.G_polar<<' '<<energy.G_nonpolar<<' '<<energy.E_elec+energy.E_vdw+energy.G_polar+energy.G_nonpolar<<' '<<ar<<' '<<al<<' '<<ac<<'\n';
 }
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
