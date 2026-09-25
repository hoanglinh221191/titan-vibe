#include "gbsa_intrin.h"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv){try{
 // argv[2] (optional) is rgbmax in Angstrom; omitted means the uncut sum, which
 // corresponds to running sander with a cutoff larger than the system.
 if(argc<2||argc>3)return 2;
 double rgbmax=std::numeric_limits<double>::infinity();
 if(argc==3){rgbmax=std::stod(argv[2]);if(!(rgbmax>0.0))return 2;}
 std::ifstream in(argv[1]);int nr,nl,igb;if(!(in>>nr>>nl>>igb)||nr<1||nl<1) return 2;
 MOL rec,lig;rec.pos.resize(nr);lig.pos.resize(nl);
 Multi_MOL_parameters rp(nr),lp(nl);std::vector<mol_info> ri,li;
 auto read=[&](MOL& mol,Multi_MOL_parameters& params,std::vector<mol_info>& info){
  for(std::size_t i=0;i<mol.size();++i){std::string name,element;double x,y,z,q,sigma,epsilon,radius,screen;
   if(!(in>>name>>element>>x>>y>>z>>q>>sigma>>epsilon>>radius>>screen))throw std::runtime_error("bad fixture");
   mol.pos.x[i]=x;mol.pos.y[i]=y;mol.pos.z[i]=z;params[i].set_data(q,radius,sigma,epsilon,12.0);
   info.emplace_back(1,"TST",name,"A",element);
  }
 };
 read(rec,rp,ri);read(lig,lp,li);
 AlignedDoubleVector rf,lf,rc,lc;
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rf,igb,0,nr,nullptr,nullptr,nullptr,rgbmax);
 cal_effective_born_radii_avx2(lig,lp.radius,li,lf,igb,0,nl,nullptr,nullptr,nullptr,rgbmax);
 cal_effective_born_radii_avx2(rec,rp.radius,ri,rc,igb,0,nr,&lig,&lp.radius,&li,rgbmax);
 cal_effective_born_radii_avx2(lig,lp.radius,li,lc,igb,0,nl,&rec,&rp.radius,&ri,rgbmax);
 arguments args;args.igb=igb;args.int_diel=1;args.ext_diel=78.5;args.saltcon=0;args.force_f="AMBER99SB";args.flag_alpb=false;
 cache_sigma_epsilon cache;GBSA_energy_snap energy;GBSASasaCache sasa;
 sasa.sasa_rec_free.assign(nr,0);sasa.sasa_lig_free.assign(nl,0);sasa.sasa_complex.assign(nr+nl,0);
 std::vector<double> unused;
 cal_energy_gbsa_avx2(rec,lig,ri,li,cache,0,nr,0,nl,energy,0,args,rp,lp,unused,unused,unused,unused,rf,lf,rc,lc,nullptr,&sasa,false);
 std::cout<<std::setprecision(17)<<energy.E_elec<<' '<<energy.E_vdw<<' '<<energy.G_polar<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
