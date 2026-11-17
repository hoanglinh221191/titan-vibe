#include "titan_path/control.h"
#include <array>
#include <cassert>
#include <iostream>
int main(){
  std::array<titan_path::SelectionCandidate,3> c;
  for(auto &x:c){x.valid=true;x.energy_key=1.;x.delta_progress=-.1;x.path_distance=2.;x.pull_distance=.25;x.minimum_clearance=2.;x.angle_deg=20.;x.navigation_enabled=true;x.navigation_cost=0.;}
  assert(!titan_path::select_epsilon_progress_candidate(c,titan_path::RouteMode::Advance,.05,false).found);
  auto r=titan_path::select_epsilon_reroute_candidate(c,.05,.02);assert(r.found&&r.index==0);
  c[0].delta_progress=-.6;r=titan_path::select_epsilon_reroute_candidate(c,.05,.02);assert(r.found&&r.index==1);
  c[2].energy_key=.9;r=titan_path::select_epsilon_reroute_candidate(c,.05,.02);assert(r.index==2);
  c[2].energy_key=1.;c[2].navigation_cost=.1;r=titan_path::select_epsilon_reroute_candidate(c,.05,.02);assert(r.index==1);
  c[2].navigation_cost=0.;c[2].delta_progress=.01;r=titan_path::select_epsilon_reroute_candidate(c,.05,.02);assert(r.index==2);
  for(auto &x:c){x.delta_progress=-1.;x.navigation_cost=0.;}
  assert(!titan_path::select_epsilon_reroute_candidate(c,.05,.02).found);
  assert(titan_path::select_epsilon_exploratory_candidate(c,.05,.02).found);
  for(auto &x:c)x.valid=false;
  assert(!titan_path::select_epsilon_exploratory_candidate(c,.05,.02).found);
  std::cout<<"REROUTE_RANK_PASS checks=9\n";
}
