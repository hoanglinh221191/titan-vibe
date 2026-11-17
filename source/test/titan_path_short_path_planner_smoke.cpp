#include "titan_path/short_path_planner.h"
#include <iostream>
#include <random>
#include <stdexcept>
using namespace titan_path;
SelectionCandidate candidate(double energy,double progress,double nav=0){
  SelectionCandidate c;c.valid=true;c.energy_key=energy;c.delta_progress=progress;
  c.pull_distance=.25;c.minimum_clearance=2.;c.angle_deg=10.;c.path_distance=1.;
  c.navigation_enabled=true;c.navigation_cost=nav;return c;
}
int main(){
  int checks=0;auto require=[&](bool ok){++checks;if(!ok)throw std::runtime_error("check "+std::to_string(checks));};
  try{
    std::vector<SelectionCandidate> flat{candidate(1,0)};
    require(!select_epsilon_progress_candidate(flat,RouteMode::Advance,.05,false,false,.02).found);
    require(select_epsilon_nonnegative_prefix_candidate(flat,.05,.02).found);
    require(rank_short_path_frontier(flat,8,false).indices.size()==1);
    require(rank_short_path_frontier(flat,8,true).indices.empty());
    flat[0].delta_progress=-.01;require(rank_short_path_frontier(flat,8,false).indices.empty());
    std::vector<SelectionCandidate> band{candidate(1,.1),candidate(1.04,.2),candidate(1.08,10)};
    auto b=rank_short_path_frontier(band,8,true);
    require(b.indices==std::vector<std::size_t>({1,0}));
    band={candidate(1,.1,0),candidate(1,.2,.019),candidate(1,10,.03)};
    require(rank_short_path_frontier(band,8,true).indices==std::vector<std::size_t>({1,0}));
    band={candidate(1,.1),candidate(1,.2)};band[1].minimum_clearance=1.49;
    require(rank_short_path_frontier(band,8,true).indices==std::vector<std::size_t>({0}));
    band[1].minimum_clearance=2.;band[1].angle_deg=21;
    require(rank_short_path_frontier(band,8,true).indices==std::vector<std::size_t>({0}));
    require(!rank_short_path_frontier(band,1,true).md_authorized);
    std::mt19937 gen(210921);std::uniform_real_distribution<double> random(0,1);
    for(int trial=0;trial<100;++trial){
      std::vector<SelectionCandidate> c;
      for(int i=0;i<30;++i)c.push_back(candidate(random(gen),.01+random(gen),random(gen)));
      auto legacy=select_epsilon_progress_candidate(c,RouteMode::Advance,.05,false,false,.02);
      auto planned=rank_short_path_frontier(c,8,true);
      require(legacy.found&&!planned.indices.empty()&&planned.indices[0]==legacy.index);
    }
    std::cout<<"SHORT_PATH_NATIVE_RANK_PASS checks="<<checks<<"\n";
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
