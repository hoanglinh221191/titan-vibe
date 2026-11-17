// DHAM input mapping: TitanDHAM::compute_state_energies maps frames to
// (state, replica) means and rejects bad indices and lengths. Since KIN-13
// the -dham estimator is TRAM; its checks are in test/titan_tram_oracle.cpp.
#include "titan_dham.h"
#include <iostream>
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
 std::unordered_map<int,std::vector<int>> traj{{0,{0,1,0}},{1,{1,0,1}}};
 std::unordered_map<int,std::vector<double>> values{{0,{2,4,6}},{1,{10,20,30}}};
 auto means=TitanDHAM::compute_state_energies(traj,values,2,2,2);
 require(means(0,0)==4&&means(1,0)==4&&means(0,1)==20&&means(1,1)==20,"state/replica mapping");
 int failures=0;
 try{auto bad=values;bad[0].pop_back();TitanDHAM::compute_state_energies(traj,bad,2,2,2);}catch(const std::invalid_argument&){++failures;}
 try{auto bad=traj;bad[0][0]=2;TitanDHAM::compute_state_energies(bad,values,2,2,2);}catch(const std::out_of_range&){++failures;}
 require(failures==2,"invalid mapping accepted");
 std::cout<<"PASS DHAM input mapping: threaded state/replica means, invalid indices and lengths rejected\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
