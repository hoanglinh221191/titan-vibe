#include "titan_llps_profile_geometry.h"
#include <iostream>
#ifdef TITAN_TEST_ENGINE
#include "titan_llps_cluster.h"
#include "arguments.h"
#include "info.h"
#endif
void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}
int main(){try {
 int cases=0;
 for(double gamma:{25.,65.,90.,120.}) {
  box Box(100.,120.,140.,90.,90.,double(gamma));TitanExactPeriodicBox cell(Box);
  for(int axis=0;axis<3;++axis) {
   std::vector<std::array<double,3>> points;std::vector<bool> dense{true,true,true,true,false};
   for(int i=0;i<5;++i){std::array<double,3> f{.15+i*.13,.25+i*.07,.12+i*.09};f[axis]=.98+(i==4?.3:(i%2?.03:-.03));points.push_back(cell.metric.scaled_to_cartesian(f[0],f[1],f[2]));}
   auto slab=TitanLLPSProfileGeometry::build(Box,points,dense,"slab",axis,10.,.02);
   require(slab.valid,"slab rejected");
   require(std::abs(slab.distance_nm[4]-.3*cell.heights[axis]/10.)<1e-9,"slab normal is not reciprocal face distance");
   require(std::abs(slab.dense_envelope_nm3-.06*Box.volume()/1000.)<1e-8,"slab envelope volume");
   double volume=0;for(double r=0;r<slab.limit_nm;r+=.17)volume+=slab.shell_volume(r,r+.17);
   require(std::abs(volume-Box.volume()/1000.)<1e-8,"partial slab shell conservation");
   for(std::size_t i=0;i<points.size();++i){auto shift=cell.metric.scaled_to_cartesian(3*(i+1),-2.*i,4.*i);for(int d=0;d<3;++d)points[i][d]+=shift[d];}
   auto moved=TitanLLPSProfileGeometry::build(Box,points,dense,"slab",axis,10.,.02);
   require(moved.valid,"image shifted slab rejected");for(std::size_t i=0;i<points.size();++i)require(std::abs(moved.distance_nm[i]-slab.distance_nm[i])<1e-8,"slab not invariant to atom images");++cases;
  }
  std::vector<std::array<double,3>> points{{0,0,0},{2,0,0},{4,0,0},{6,0,0},{3,10,0}};std::vector<bool> dense{true,true,true,true,false};
  auto drop=TitanLLPSProfileGeometry::build(Box,points,dense,"droplet",2,2.1,.02);require(drop.valid,"finite droplet rejected");
  require(std::abs(drop.radius_nm-.3)<1e-10&&std::abs(drop.distance_nm[4]-1.)<1e-10,"finite droplet center/shape");
  double volume=0;for(double r=0;r<drop.limit_nm;r+=.17)volume+=drop.shell_volume(r,r+.17);
  require(std::abs(volume-4*std::numbers::pi*std::pow(drop.limit_nm,3)/3)<1e-8,"truncated radial sphere volume");
  require(drop.shell_volume(drop.limit_nm,drop.limit_nm+1)==0,"shell beyond safe domain");++cases;
#ifdef TITAN_TEST_ENGINE
  LLPSPaperAccumulator accumulator;LLPSPaperConfig config;config.enabled=true;config.geometry="droplet";config.profile_bin_nm=.17;config.connectivity_cutoff_angstrom=2.1;accumulator.configure(config);
  std::vector<LLPSNode> nodes;
  for(std::size_t i=0;i<points.size();++i){LLPSNode n;n.id=i;n.type=1;n.x=points[i][0];n.y=points[i][1];n.z=points[i][2];nodes.push_back(n);}
  accumulator.process_frame(0,0,nodes,dense,1,Box);
  double count=0,actual_volume=0;
  for(const auto& b:accumulator.profile_bins_){count+=b.protein_count;actual_volume+=b.volume_nm3_sum;}
  require(count==5&&std::abs(actual_volume-volume)<1e-8,"production profile mass/volume");
  require(accumulator.frames_.size()==1&&accumulator.geometry_failures_.empty(),"production frame accounting");
#endif
 }
 box b(20.,20.,20.,90.,90.,90.);std::vector<std::array<double,3>> winding{{0,0,0},{5,0,0},{10,0,0},{15,0,0}};std::vector<bool> dense(4,true);
 require(!TitanLLPSProfileGeometry::build(b,winding,dense,"droplet",2,5.1,.02).valid,"winding droplet accepted");
 box tilted(20.,20.,20.,90.,90.,65.);require(TitanLLPSProfileGeometry::build(tilted,{{0,0,0}},{true},"auto",2,5.,.02).reason=="triclinic_requires_explicit_geometry","implicit tilted geometry accepted");
 require(!TitanLLPSProfileGeometry::build(b,{{0,0,0},{0,0,10}},{true,true},"slab",2,5.,.02).valid,"ambiguous slab center accepted");
 std::cout<<"PASS paper geometry "<<cases<<" lattice cases, image invariance, shell conservation, finite/winding/ambiguous domains";
#ifdef TITAN_TEST_ENGINE
 std::cout<<", production mass/volume";
#endif
 std::cout<<'\n';return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
