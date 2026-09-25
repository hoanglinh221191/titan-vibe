#include "titan_topology_energy.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

void require(bool good,const char *what){if(!good)throw std::runtime_error(what);}
void close(double a,double b,const char *what,double tol=1e-10){require(std::isfinite(a)&&std::abs(a-b)<=tol*std::max(1.0,std::abs(b)),what);}
std::shared_ptr<const titan_nb::Topology> load(const std::string &name,const std::string &text){std::ofstream(name)<<text;return titan_nb::read_topology(name);}
std::string make(int rule,std::string definitions={},std::string local={},bool gen=true,int nrexcl=3) {
  return "[ defaults ]\n1 "+std::to_string(rule)+(gen?" yes":" no")+" 0.5 0.8333333333333333\n[ atomtypes ]\n"+
    (rule==1?"X 6 12 0 A 0.004 0.000012\nY 6 12 0 A 0.009 0.000027\n":"X 6 12 0 A 0.3 0.2\nY 6 12 0 A 0.4 0.8\n")+
    definitions+"[ moleculetype ]\nChain "+std::to_string(nrexcl)+"\n[ atoms ]\n1 X 1 T C1 1 0.2 12\n2 X 1 T C2 1 0.1 12\n3 Y 1 T C3 1 -0.1 12\n4 Y 1 T C4 1 -0.4 12\n5 X 1 T C5 1 0.5 12\n"+local+"[ system ]\nSynthetic\n[ molecules ]\nChain 2\n";
}
void rejects(const std::string &text,const char *why){bool failed=false;try{load("invalid.top",text);}catch(const std::exception&){failed=true;}require(failed,why);}
int main(){try{
  const std::string bonds="[ bonds ]\n1 2 1\n2 3 1\n3 4 1\n";
  for(int rule:{1,2,3}){
    auto top=load("chain.top",make(rule,{},bonds+"[ pairs ]\n1 4 1\n"));
    require(top->atoms.size()==10,"molecule copies");
    require(top->exclusions.contains(titan_nb::pair_key(0,1))&&top->exclusions.contains(titan_nb::pair_key(0,2)),"1-2/1-3 exclusion");
    require(top->pair(0,3).normal.c6==0,"1-4 excluded normal pair");
    require(!top->exclusions.contains(titan_nb::pair_key(0,5)),"copies must not share exclusions");
    require(!top->exclusions.contains(titan_nb::pair_key(0,4)),"near nonbonded atom must remain");
    const double sigma=rule==2?0.35:std::sqrt(0.3*0.4);
    const double c6=rule==1?0.006:4*0.4*std::pow(sigma,6);
    const double c12=rule==1?0.000018:4*0.4*std::pow(sigma,12);
    close(top->pair(0,3).listed.c6,0.5*c6,"generated 1-4 C6");
    close(top->pair(0,3).listed.c12,0.5*c12,"generated 1-4 C12");
    close(top->pair(0,3).listed_charge,-0.08/1.2,"fudgeQQ");
    cache_sigma_epsilon cache;cache.nonbonded_topology=top;
    arguments args;args.force_f="AMBER99SB";args.self_int=true;args.epsilon_in=2;
    std::vector<mol_info> info;
    for(int i=0;i<5;++i){info.emplace_back(1,"C","C"+std::to_string(i+1),"A","T");info.back().set_original_atom_id(i);}
    titan_nb::validate_selection(*top,info);
    const double r=0.45;
    auto value=titan_nb::evaluate_pair(cache,args,info[0],info[3],r,0.2);
    close(value.elec,138.935457839*(-0.08/1.2)/(2*r),"listed pair ignores normal cutoff");
    close(value.vdw,0.5*c12/std::pow(r,12)-0.5*c6/std::pow(r,6),"1-4 LJ energy");
    const double h=1e-6;
    auto left=titan_nb::evaluate_pair(cache,args,info[0],info[3],r-h,0.2);
    auto right=titan_nb::evaluate_pair(cache,args,info[0],info[3],r+h,0.2);
    close(value.force_vdw_prefactor*r,-(right.vdw-left.vdw)/(2*h),"LJ force derivative",1e-7);
    close(value.force_elec_prefactor*r,-(right.elec-left.elec)/(2*h),"Coulomb force derivative",1e-7);
    auto near=titan_nb::evaluate_pair(cache,args,info[0],info[4],0.15,-1);
    require(near.vdw>0 && near.elec>0,"distance threshold must not remove nonbonded short contact");
    args.flag_elec_screening=true;args.salt_concentration=0.15;
    close(titan_nb::evaluate_pair(cache,args,info[0],info[3],r,-1).elec,value.elec,"1-4 always plain Coulomb");
  }
  auto typed=load("typed.top",make(1,"[ pairtypes ]\nX Y 1 0.02 0.00007\n",bonds+"[ pairs ]\n1 4 1\n",false));
  close(typed->pair(0,3).listed.c6,0.02,"explicit pairtypes must not be scaled twice");
  auto nbfix=load("nbfix.top",make(1,"[ nonbond_params ]\nX Y 1 0.02 0.00007\n",bonds+"[ pairs ]\n1 4 1\n"));
  close(nbfix->pair(0,3).listed.c6,0.01,"generated 1-4 must use resolved nonbonded override");
  auto explicit_pair=load("explicit.top",make(1,{},bonds+"[ pairs ]\n1 4 1 0.04 0.00009\n",false));
  close(explicit_pair->pair(0,3).listed.c12,0.00009,"explicit pair parameters take precedence");
  auto listed=load("listed.top",make(1,{},"[ pairs ]\n1 4 2 0.3 0.6 -0.7 0.04 0.00009\n",false,0));
  close(listed->pair(0,3).listed_charge,-0.126,"function 2 explicit charges");
  require(listed->pair(0,3).normal.c6>0,"listed pairs are additive when normal pair is not excluded");
  auto nb=load("pairs_nb.top",make(1,{},"[ pairs_nb ]\n1 4 1 0.6 -0.7 0.04 0.00009\n[ exclusions ]\n1 4 5\n",false,0));
  close(nb->pair(0,3).listed_charge,-0.42,"pairs_nb unscaled charges");
  require(nb->pair(0,3).normal.c6==0 && nb->pair(0,4).normal.c6==0,"manual exclusions");
  auto connections=load("connect.top",make(1,{},"[ bonds ]\n1 2 6\n2 3 5\n[ constraints ]\n3 4 2\n4 5 1\n",true,1));
  require(!connections->exclusions.contains(titan_nb::pair_key(0,1)),"bond type 6 is not chemical");
  require(connections->exclusions.contains(titan_nb::pair_key(1,2)),"bond type 5 generates exclusion");
  require(!connections->exclusions.contains(titan_nb::pair_key(2,3)),"constraint type 2 is not chemical");
  require(connections->exclusions.contains(titan_nb::pair_key(3,4)),"constraint type 1 generates exclusion");
  rejects(make(1,{},bonds+"[ pairs ]\n1 4 1\n",false),"missing pairtype must reject");
  rejects(make(1,{},"[ exclusions ]\n1 80\n"),"out-of-range exclusions must reject");
  rejects(make(1,{},"[ virtual_sites2 ]\n5 1 2 1 0.5\n"),"unsupported virtual-site exclusions must reject");
  rejects("#if 1\n"+make(1)+"#endif\n","unknown conditional expression must reject");
  std::ofstream("included.itp")<<make(1,{},"[ exclusions ]\n1 5\n");
  auto included=load("include.top","#define USE_A\n#ifdef USE_A\n#include \"included.itp\"\n#else\n#error wrong_branch\n#endif\n");
  require(included->exclusions.contains(titan_nb::pair_key(0,4)),"include/conditional source ordering");
  std::cout<<"PASS topology exclusions, 3 combination rules, generated/explicit/pairtypes 1-4, listed pair additivity, function 2/pairs_nb, force derivatives, near nonbonded contacts and reject guards\n";
  return 0;
}catch(const std::exception &error){std::cerr<<"FAIL "<<error.what()<<'\n';return 1;}}
