#include "titan_topology_energy.h"
#include "titan_energy_decomposition.h"
#include "functions.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

int force_energy_add(const cache_sigma_epsilon&,size_t,double&,double&,force&,force&,
                     const Multi_MOL_parameters&,const Multi_MOL_parameters&,
                     const Multi_MOL&,const Multi_MOL&,size_t,const arguments&,size_t,size_t);
int force_energy_add_intrin(const cache_sigma_epsilon&,size_t*,double&,double&,force&,force&,
                            const Multi_MOL_parameters&,const Multi_MOL_parameters&,
                            const Multi_MOL&,const Multi_MOL&,size_t,const arguments&,size_t,size_t);

void close(double a,double b,const char *message){if(!std::isfinite(a)||std::abs(a-b)>1e-9*std::max(1.0,std::abs(b)))throw std::runtime_error(message);}
int main(int argc,char **argv){try{
  if(argc!=2)return 2;
  Multi_MOL rec,lig;rec.m_mol.resize(1);lig.m_mol.resize(1);
  auto &r=rec.m_mol[0];auto &l=lig.m_mol[0];r.pos.resize(5);
  const double xyz[5][3]={{10,10,10},{11.5,10,10},{13,10,10},{14.5,10,10},{11,11.5,10.2}};
  for(int i=0;i<5;++i){r.pos.x[i]=xyz[i][0];r.pos.y[i]=xyz[i][1];r.pos.z[i]=xyz[i][2];
    rec.info.emplace_back(1,"C","C"+std::to_string(i+1),"A","T");rec.info.back().set_original_atom_id(i);}
  l=r;lig.info=rec.info;rec.parameters.resize(5);lig.parameters.resize(5);
  cache_sigma_epsilon cache;cache.nonbonded_topology=titan_nb::read_topology(argv[1]);
  arguments args;args.self_int=true;args.epsilon_in=2;args.distance_energy_cutoff=4;
  int calls=0;
  for(const std::string ff:{"AMBER99SB","CHARMM36","OPLSAA","GROMOS54A7","MARTINI"}) {
    args.force_f=ff;cache.force_field_name=ff;
    double ee=0,ev=0;double fx[5]{},fy[5]{},fz[5]{};
    for(int i=0;i<5;++i)for(int j=i+1;j<5;++j) {
      const double dx=.1*(xyz[i][0]-xyz[j][0]),dy=.1*(xyz[i][1]-xyz[j][1]),dz=.1*(xyz[i][2]-xyz[j][2]);
      auto v=titan_nb::evaluate_pair(cache,args,rec.info[i],rec.info[j],std::sqrt(dx*dx+dy*dy+dz*dz),4);
      ee+=v.elec;ev+=v.vdw;const double f=v.force_elec_prefactor+v.force_vdw_prefactor;
      fx[i]+=f*dx;fy[i]+=f*dy;fz[i]+=f*dz;fx[j]-=f*dx;fy[j]-=f*dy;fz[j]-=f*dz;
    }
    for(bool cached:{false,true}) {
      cache.elec_ener.clear();cache.vdw_ener.clear();cache.opt_forces.reset();
      if(cached){cache.elec_ener.assign(25,0);cache.vdw_ener.assign(25,0);cache.opt_forces.emplace();cache.opt_forces->assign(25);}
      cache_energy energy;cache_force force;
      auto check_energy=[&](){close(energy.elec_ener,2*ee,"dispatch electro energy");close(energy.vdw_ener,2*ev,"dispatch LJ energy");++calls;};
#define CHECK_ENERGY(function) function(r,l,rec.info,lig.info,cache,0,5,energy,0,args,rec.parameters,lig.parameters);check_energy();
      CHECK_ENERGY(cal_energy_amber)
      CHECK_ENERGY(cal_energy_amber_screening)
      CHECK_ENERGY(cal_energy_opls)
      CHECK_ENERGY(cal_energy_gromos)
      CHECK_ENERGY(cal_energy_martini)
      CHECK_ENERGY(cal_energy_amber_intrin)
      CHECK_ENERGY(cal_energy_amber_screening_intrin)
      CHECK_ENERGY(cal_energy_opls_intrin)
      CHECK_ENERGY(cal_energy_gromos_intrin)
      CHECK_ENERGY(cal_energy_martini_intrin)
      CHECK_ENERGY(cal_energy_amber_intrin_nofilter)
      CHECK_ENERGY(cal_energy_amber_screening_intrin_nofilter)
      CHECK_ENERGY(cal_energy_opls_intrin_nofilter)
      CHECK_ENERGY(cal_energy_gromos_intrin_nofilter)
#undef CHECK_ENERGY
      auto check_forces=[&](){check_energy();if(cached)for(int i=0;i<5;++i){double x=0,y=0,z=0;for(int j=0;j<5;++j){auto f=cache.get_force(5*j+i);x+=f->forces[0].Fx+f->forces[1].Fx;y+=f->forces[0].Fy+f->forces[1].Fy;z+=f->forces[0].Fz+f->forces[1].Fz;}close(x,fx[i],"cached atom force x");close(y,fy[i],"cached atom force y");close(z,fz[i],"cached atom force z");}};
#define CHECK_FORCE(function) function(r,l,rec.info,lig.info,cache,0,5,energy,force,0,args,rec.parameters,lig.parameters);check_forces();
      CHECK_FORCE(cal_energy_force_amber)
      CHECK_FORCE(cal_energy_force_opls)
      CHECK_FORCE(cal_energy_force_gromos)
      CHECK_FORCE(cal_energy_force_martini)
#undef CHECK_FORCE
      EnergyDecompositionProvider provider;provider.configure(rec,lig,cache,args,4,true);
      double de=0,dv=0;
      for(int i=0;i<5;++i)for(int j=0;j<5;++j){auto v=provider.pair_energy(r,l,i,j);de+=v.elec;dv+=v.vdw;}
      close(de,ee,"fallback decomposition self double-count correction");close(dv,ev,"fallback decomposition LJ closure");
      double intra=0;for(int i=0;i<5;++i)for(int j=i+1;j<5;++j){auto v=provider.rec_pair_energy(r,i,j);intra+=v.elec+v.vdw;}
      close(intra,ee+ev,"intra-group unique pair energy");
      if(cached){for(auto &v:cache.elec_ener)v*=0.5;for(auto &v:cache.vdw_ener)v*=0.5;}
      double scalar_e=0,scalar_v=0,simd_e=0,simd_v=0;
      ::force scalar_fe(0.0,0.0,0.0),scalar_fv(0.0,0.0,0.0),simd_fe(0.0,0.0,0.0),simd_fv(0.0,0.0,0.0);
      for(size_t j=0;j<4;++j)force_energy_add(cache,5*j,scalar_e,scalar_v,scalar_fe,scalar_fv,rec.parameters,lig.parameters,rec,lig,0,args,0,j);
      size_t indices[4]={0,5,10,15};
      force_energy_add_intrin(cache,indices,simd_e,simd_v,simd_fe,simd_fv,rec.parameters,lig.parameters,rec,lig,0,args,0,0);
      close(scalar_e,simd_e,"legacy scalar/SIMD subgroup electro");close(scalar_v,simd_v,"legacy scalar/SIMD subgroup LJ");
      close(scalar_fe.Fx,simd_fe.Fx,"legacy subgroup electro force");close(scalar_fv.Fx,simd_fv.Fx,"legacy subgroup LJ force");
      auto expected=titan_nb::evaluate_pair(cache,args,rec.info[0],rec.info[3],0.45,4);
      close(scalar_e,0.5*expected.elec,"legacy subgroup 1-4 energy factor");
      close(scalar_v,0.5*expected.vdw,"legacy subgroup 1-4 LJ factor");
      calls+=2;
    }
  }
  cache.nonbonded_topology.reset();bool rejected=false;
  try{cache_energy e;cal_energy_amber(r,l,rec.info,lig.info,cache,0,5,e,0,args,rec.parameters,lig.parameters);}catch(const std::exception&){rejected=true;}
  if(!rejected)throw std::runtime_error("unbound self kernel must reject");
  std::cout<<"PASS "<<calls<<" topology dispatch routes, per-atom forces, cache/fallback decomposition and unbound rejection\n";return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
