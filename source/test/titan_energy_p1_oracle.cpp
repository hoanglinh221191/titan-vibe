#include "titan_gist.h"
#include "titan_llps.h"
#include "functions.h"
#include "gromos_nonbonded.h"
#include <iostream>
#include <stdexcept>

void check(double got,double expected,const char *label,double tol=1e-10) {
  if(!std::isfinite(got) || std::abs(got-expected)>tol*std::max(1.0,std::abs(expected))) {
    std::cerr<<label<<" got="<<got<<" expected="<<expected<<'\n';
    throw std::runtime_error(label);
  }
}

void gist_oracle() {
  MOL ligand;ligand.pos.resize(1);ligand.pos.x[0]=100;
  Multi_MOL_parameters rec(2),lig(1);
  rec.charge[0]=1;rec.sigma[0]=0.3;rec.epsilon[0]=0.25;
  rec.charge[1]=-0.2;
  std::vector<mol_info> ri{mol_info(1,"C","CA","A","ALA"),mol_info(2,"H","H","A","ALA")};
  std::vector<mol_info> li{mol_info(1,"C","C","B","LIG")};
  TitanGIST gist;
  gist.init(ligand,1,0.5,300,9,1,"TIP3P",rec,ri,lig,li);
  if(gist.n_solute_real_!=3)throw std::runtime_error("GIST must include hydrogen charges");
  gist.solute_x_[1]=1;gist.solute_x_[2]=100;
  const double sigma=0.5*(3.1507+3.0),eps=std::sqrt(0.636*0.25);
  check(gist.solute_sigma_half_[0],sigma,"GIST nm-to-Angstrom mixing");
  check(gist.solute_eps_mixed_[0],eps,"GIST epsilon mixing");
  const double electro=1389.3546*(-0.834/4.0+0.834/std::sqrt(17.0)
      -0.2*(-0.834/3.0+0.834/std::sqrt(10.0)));
  const double sr6=std::pow(sigma/4.0,6);
  const double lj=4*eps*(sr6*sr6-sr6);
  check(gist.compute_esw(4,0,0,4,1,0,4,-1,0),electro+lj,"GIST water/solute analytic pair");
  std::cout<<"PASS GIST Coulomb units, hydrogen charges, sigma conversion, epsilon and LJ energy\n";
}

void gromos_oracle(bool sparse) {
  MOL rec,lig;rec.pos.resize(5);lig.pos.resize(1);
  Multi_MOL_parameters rp(5),lp(1);
  for(int i=0;i<5;++i) {
    rec.pos.x[i]=3+i;rp.set_type_id(i,0);rp.charge[i]=0.4;
    rp.sigma[i]=0.02;rp.epsilon[i]=0.03;
  }
  lp.set_type_id(0,1);lp.charge[0]=-0.5;lp.sigma[0]=0.01;lp.epsilon[0]=0.02;
  std::vector<mol_info> ri(5),li(1);
  arguments args;args.force_f="GROMOS54A7";args.self_int=false;
  args.distance_energy_cutoff=0;args.epsilon_in=2;
  cache_sigma_epsilon cache;cache.force_field_name=args.force_f;
  const double c6=0.004,c12=0.000012;
  if(sparse)cache.set_charmm_pair_sparse({charmm::make_pair_override_key(0,1)},{c6},{c12});
  else cache.set_charmm_pair_table({-1,c6,c6,-1},{-1,c12,c12,-1},2);
  double ev=0,ee=0,fx=0;
  for(int i=0;i<5;++i) {
    const double r=0.1*rec.pos.x[i];
    ev+=c12/std::pow(r,12)-c6/std::pow(r,6);
    ee+=138.935457839*0.4*(-0.5)/(2*r);
    fx+=12*c12/std::pow(r,13)-6*c6/std::pow(r,7);
  }
  for(bool cached:{false,true}) {
    if(cached){cache.elec_ener.assign(5,0);cache.vdw_ener.assign(5,0);}
    cache_energy e;
    cal_energy_gromos(rec,lig,ri,li,cache,0,1,e,0,args,rp,lp);
    check(e.vdw_ener,ev,"GROMOS scalar override");check(e.elec_ener,ee,"GROMOS scalar electro");
    cal_energy_gromos_intrin(rec,lig,ri,li,cache,0,1,e,0,args,rp,lp);
    check(e.vdw_ener,ev,"GROMOS AVX + tail override");
    cal_energy_gromos_intrin_nofilter(rec,lig,ri,li,cache,0,1,e,0,args,rp,lp);
    check(e.vdw_ener,ev,"GROMOS nofilter AVX + tail override");
  }
  cache.elec_ener.clear();cache.vdw_ener.clear();
  cache_energy e;cache_force f;
  cal_energy_force_gromos(rec,lig,ri,li,cache,0,1,e,f,0,args,rp,lp);
  check(e.vdw_ener,ev,"GROMOS force-path energy");check(f.forces[1].Fx,fx,"GROMOS force derivative");
  std::cout<<"PASS GROMOS explicit C6/C12 scalar/AVX/tail/nofilter/cache/force\n";
}

void llps_oracle() {
  check(TitanLLPS::pbc_dist2(20.2f,-30.3f,40.4f,10,10,10),0.29,
        "LLPS translation-invariant distance",1e-5);
  TitanLLPS llps(1,{0},{1});
  MOL frame;frame.pos.resize(1);
  frame.Box=box(10.0,10.0,10.0,90.0,90.0,180.0);
  std::vector<mol_info> info(1);
  bool sync=false,async=false;
  try{llps.process_frame(frame,info,1,1,1,0,0);}catch(const std::invalid_argument&){sync=true;}
  llps.start_async();
  try{llps.process_frame_async(frame,info,1,1,1,0,0);}catch(const std::invalid_argument&){async=true;}
  if(!sync||!async)throw std::runtime_error("LLPS invalid-cell rejection must cover sync and async");
  TitanLLPS spherical(2,{0,1},{1,1});
  std::vector<mol_info> beads{mol_info(1,"C","BB","A","ALA"),mol_info(2,"C","BB","B","ALA")};
  spherical.set_sticker_beads(beads,"BB");
  spherical.configure_phase({"droplet",2,.5,1,2});
  MOL pair;pair.pos.resize(2);pair.pos.x[0]=-1;pair.pos.x[1]=1;
  pair.Box=box(100.0,100.0,100.0,90.0,90.0,90.0);
  spherical.process_frame(pair,beads,10,10,10,0,0);
  const auto &r=spherical.result().frames.back();
  const double radius=.5;
  check(r.gcc_rg_nm,0.1,"LLPS two-point Rg",1e-6);
  check(r.density_dense,2/(4.0/3.0*std::acos(-1.0)*std::pow(radius,3)),"LLPS declared-window density",1e-5);
  spherical.configure_phase({});
  pair.pos.x[1]=30;
  spherical.process_frame(pair,beads,10,10,10,1,1);
  if(!std::isnan(spherical.result().frames.back().density_dense))
    throw std::runtime_error("LLPS unconfigured phase must have no density estimate");
  std::cout<<"PASS LLPS multiple images, invalid-cell rejection and declared phase windows\n";
}

void cli_oracle() {
  arguments args;
  if(args.ml_predict_experimental)throw std::runtime_error("RF prediction must require opt-in");
  std::vector<std::string> tokens{"titan","-qha-temp","315","-ml-predict-experimental","true","-tica-backend","exact_streaming"};
  std::vector<char*> argv;for(auto &token:tokens)argv.push_back(token.data());
  if(args.parse_arguments(static_cast<int>(argv.size()),argv.data())!=0 ||
     args.qha_temperature!=315 || !args.ml_predict_experimental || !args.flag_ml ||
     args.msm.tica_backend!="exact_streaming")
    throw std::runtime_error("P1 options must reach effective arguments");
  std::cout<<"PASS effective CLI QHA temperature, exact streaming and RF opt-in\n";
}

int main() {
  std::cout<<std::unitbuf;std::cerr<<std::unitbuf;
  try{gist_oracle();gromos_oracle(false);gromos_oracle(true);llps_oracle();cli_oracle();return 0;}
  catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
