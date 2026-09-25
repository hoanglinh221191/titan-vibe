#include "titan_exact_periodic_box.h"
#include <iostream>
#include <random>
int main(){try{
 std::mt19937 rng(33);std::uniform_real_distribution<double> u(-3,3);int checked=0;
 for(double gamma:{25.,65.,90.,120.}) {
  box b(10.,12.,14.,90.,90.,double(gamma));TitanExactPeriodicBox cell(b);
  for(int n=0;n<200;++n){auto p=cell.metric.scaled_to_cartesian(u(rng),u(rng),u(rng));double best=1e100;
   for(int i=-5;i<=5;++i)for(int j=-5;j<=5;++j)for(int k=-5;k<=5;++k){auto q=cell.metric.scaled_to_cartesian(i,j,k);double d2=0;for(int v=0;v<3;++v)d2+=(p[v]-q[v])*(p[v]-q[v]);best=std::min(best,d2);}
   cell.apply_mic(p[0],p[1],p[2]);if(std::abs(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]-best)>1e-8)throw std::runtime_error("exact MIC mismatch");++checked;
  }
 }
 TitanExactPeriodicBox cell(box(20.,20.,20.,90.,90.,90.));std::vector<std::array<double,3>> lifted;
 if(titan_lift_radius_graph(cell,{{0,0,0},{5,0,0},{10,0,0},{15,0,0}},5.1,lifted))throw std::runtime_error("winding not rejected");
 if(!titan_lift_radius_graph(cell,{{0,0,0},{25,0,0},{50,0,0}},5.1,lifted)||std::abs(lifted[2][0]-10)>1e-9)throw std::runtime_error("finite graph lift");
 bool rejected=false;try{cell.apply_mic(lifted[0][0]=std::numeric_limits<double>::quiet_NaN(),lifted[0][1],lifted[0][2]);}catch(const std::runtime_error&){rejected=true;}
 if(!rejected)throw std::runtime_error("invalid coordinate accepted");
 std::cout<<"PASS exact periodic adapter "<<checked<<" exhaustive MIC checks, finite/winding graph and nonfinite rejection\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
