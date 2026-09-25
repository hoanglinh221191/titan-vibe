#include "gbsa_intrin.h"
#include "titan_gbnsr6.h"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) { try {
  if (argc != 2) return 2;
  std::ifstream in(argv[1]);
  size_t nr,nl; int igb; double salt,cut,spacing;
  if (!(in>>nr>>nl>>igb>>salt>>cut>>spacing) || !nr || !nl) return 2;
  MOL rec,lig; rec.pos.resize(nr); lig.pos.resize(nl);
  Multi_MOL_parameters rp(nr),lp(nl); std::vector<mol_info> ri,li;
  const auto read=[&](MOL &m,Multi_MOL_parameters &p,std::vector<mol_info> &info) {
    for(size_t i=0;i<m.size();++i) {
      std::string name,type,res; size_t resnr; double x,y,z,q,radius;
      if(!(in>>name>>type>>res>>resnr>>x>>y>>z>>q>>radius)) throw std::runtime_error("bad fixture");
      m.pos.x[i]=x; m.pos.y[i]=y; m.pos.z[i]=z;
      p[i].set_data(q,radius,.3,.1,12.);
      info.emplace_back(resnr,type,name,"A",res);
    }
  };
  read(rec,rp,ri); read(lig,lp,li);
  AlignedDoubleVector rf,lf,rc,lc;
  const double cutoff=cut>0 ? cut : std::numeric_limits<double>::infinity();
  if (igb==66) {
    calculate_gbnsr6_binding_radii(rec,lig,rp.radius,lp.radius,ri,li,rf,lf,rc,lc,spacing);
  } else {
  cal_effective_born_radii_avx2(rec,rp.radius,ri,rf,igb,0,nr,nullptr,nullptr,nullptr,cutoff,spacing);
  cal_effective_born_radii_avx2(lig,lp.radius,li,lf,igb,0,nl,nullptr,nullptr,nullptr,cutoff,spacing);
  cal_effective_born_radii_avx2(rec,rp.radius,ri,rc,igb,0,nr,&lig,&lp.radius,&li,cutoff,spacing);
  cal_effective_born_radii_avx2(lig,lp.radius,li,lc,igb,0,nl,&rec,&rp.radius,&ri,cutoff,spacing);
  }
  arguments args; args.igb=igb; args.int_diel=1.; args.ext_diel=78.5;
  args.saltcon=salt; args.force_f="AMBER99SB"; args.flag_alpb=igb==66;
  cache_sigma_epsilon cache; GBSA_energy_snap energy; GBSASasaCache sasa;
  sasa.sasa_rec_free.assign(nr,0.); sasa.sasa_lig_free.assign(nl,0.);
  sasa.sasa_complex.assign(nr+nl,0.); std::vector<double> unused;
  cal_energy_gbsa_avx2(rec,lig,ri,li,cache,0,nr,0,nl,energy,0,args,rp,lp,
    unused,unused,unused,unused,rf,lf,rc,lc,nullptr,&sasa,false);
  std::cout<<std::setprecision(17)<<"ENERGY "<<energy.E_elec<<' '<<energy.G_polar<<'\n';
  for (const auto *r : {&rf,&lf,&rc,&lc}) {
    std::cout<<"RADII"; for(double v:*r) std::cout<<' '<<v*10.; std::cout<<'\n';
  }
  return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<'\n'; return 1; } }
