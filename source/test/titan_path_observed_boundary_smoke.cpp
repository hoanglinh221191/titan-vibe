#include "titan_path/observed_boundary.h"
#include <iostream>
#include <sstream>

int main() {
  int cases=0;
  const std::string good="SCHEMA=TITAN_OBSERVED_BOUNDARY_V1\nSTATE=OBSERVATION_VALID\nROUTE_ID=test\nROUTE_FILE_SHA256="+std::string(64,'a')+
      "\nREFERENCE_SHA256="+std::string(64,'b')+"\nOBSERVATION_SHA256="+std::string(64,'c')+
      "\nLAST_APPLIED_MODE=Recovery\nOBSERVED_PROGRESS_A=3.5\nTIME_PS=5\nPATH_TOTAL_A=10\n";
  auto parse=[&](std::string s) {std::istringstream in(s);return titan_path::read_observed_boundary(in,"test",std::string(64,'a'),std::string(64,'b'),10,5);};
  auto check=[&](bool ok){++cases;if(!ok)throw std::runtime_error("case "+std::to_string(cases));};
  try {
    check(parse(good).progress_angstrom==3.5);
    check(parse(good).mode==titan_path::RouteMode::Recovery);
    for (auto change : {std::pair{"OBSERVATION_VALID","HOLD"},std::pair{"TIME_PS=5","TIME_PS=4"},
                        std::pair{"OBSERVED_PROGRESS_A=3.5","OBSERVED_PROGRESS_A=11"},
                        std::pair{"OBSERVED_PROGRESS_A=3.5","OBSERVED_PROGRESS_A=nan"},
                        std::pair{"ROUTE_ID=test","ROUTE_ID=other"},
                        std::pair{"LAST_APPLIED_MODE=Recovery","LAST_APPLIED_MODE=Transit"},
                        std::pair{"TITAN_OBSERVED_BOUNDARY_V1","TITAN_PATH_STATUS_V1"}}) {
      auto s=good; s.replace(s.find(change.first),std::string(change.first).size(),change.second);
      bool rejected=false;try{(void)parse(s);}catch(const std::exception&){rejected=true;}check(rejected);
    }
    bool rejected=false;try{(void)parse(good+"TIME_PS=5\n");}catch(const std::exception&){rejected=true;}check(rejected);
    std::cout<<"OBSERVED_BOUNDARY_SMOKE_PASS cases="<<cases<<" no_md=true\n";
  }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
