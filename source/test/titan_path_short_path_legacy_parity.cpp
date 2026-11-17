// Compile once with frozen R3 control.cpp, once with the new shared backend.
#include "titan_path/control.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
int main(){
  using namespace titan_path;std::mt19937 rng(210921);std::uniform_real_distribution<double> r(0,1);
  std::cout<<std::setprecision(17);
  for(int trial=0;trial<1000;++trial){
    std::vector<SelectionCandidate> c(30);
    for(auto &x:c){x.valid=r(rng)>.1;x.energy_key=r(rng)*10;x.delta_progress=r(rng)*2-.5;
      x.pull_distance=.25;x.minimum_clearance=1.5+r(rng);x.angle_deg=r(rng)*20;
      x.probe_length=.25;x.fixed_probe=true;x.path_distance=r(rng)*4;
      x.reference_path_distance=r(rng)*4;x.navigation_enabled=true;x.navigation_cost=r(rng);}
    auto mode=trial%2?RouteMode::Advance:RouteMode::Recovery;
    auto v=select_epsilon_progress_candidate(c,mode,.05,trial%3==0,trial%5==0,
        mode==RouteMode::Advance&&trial%7==0?std::optional<double>(.02):std::nullopt);
    std::cout<<trial<<' '<<v.found<<' '<<v.index<<' '<<v.minimum_energy_key<<' '
             <<v.valid_count<<' '<<v.mode_eligible_count<<' '<<v.selectable_count<<' '
             <<v.eligible_count<<' '<<v.navigation_eligible_count<<' '<<v.minimum_navigation_cost<<'\n';
  }
}
