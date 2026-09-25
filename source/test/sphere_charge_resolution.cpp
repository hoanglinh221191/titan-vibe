#include "tabi/tabi_charge_resolution.h"
#include "titan_executor.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace TitanPBE::Tabi;
namespace {
double edge(const Panel&p){
  auto d=[](double x,double y,double z){return std::sqrt(x*x+y*y+z*z);};
  return std::max({d(p.v0x-p.v1x,p.v0y-p.v1y,p.v0z-p.v1z),d(p.v0x-p.v2x,p.v0y-p.v2y,p.v0z-p.v2z),d(p.v1x-p.v2x,p.v1y-p.v2y,p.v1z-p.v2z)});
}
void print(const TitanPBE::ChargeResolutionDiagnostic &r,std::size_t n){
  std::cout<<std::setprecision(17)<<"{\"status\":\""<<r.status<<"\",\"min_d_over_h\":"<<r.min_ratio
    <<",\"distance\":"<<r.distance_a<<",\"edge\":"<<r.local_edge_a
    <<",\"min_d_over_atom_radius\":"<<r.min_distance_over_atom_radius
    <<",\"inside_atom_count\":"<<r.surface_inside_atom_count<<",\"charged_atoms\":"<<r.charged_atoms
    <<",\"panels\":"<<n<<",\"patch_evaluations\":"<<r.patch_evaluations<<",\"ms\":"<<r.milliseconds<<"}\n";
}
}
int main(int argc,char**argv){try{
  TitanExecutor executor(4);Options o;o.surface_mode="vdw";o.surface_only=true;o.surface_point_count=96;
  if(argc>1){
    std::ifstream input(argv[1]);std::size_t n{};if(!(input>>n))return 2;
    std::vector<Atom>a(n);for(auto &x:a)if(!(input>>x.x>>x.y>>x.z>>x.radius>>x.q))return 2;
    if(argc>2)o.surface_point_count=std::stoi(argv[2]);
    auto p=build_surface(a,o,nullptr,&executor);print(diagnose_charge_resolution(a,p,o),p.size());return 0;
  }
  int checks=0;auto require=[&](bool pass,const char*message){++checks;if(!pass)throw std::runtime_error(message);};
  Atom sphere{0,0,0,0,2};Panel octant;octant.has_triangle_geometry=true;
  octant.v0x=2;octant.v1y=2;octant.v2z=2;
  for(int reverse=0;reverse<2;++reverse){
    require(std::abs(charge_to_spherical_panel_distance({0,0,0,1,0},sphere,octant)-2)<1e-13,"curved center distance is radius, not flat-plane distance");
    require(std::abs(charge_to_spherical_panel_distance({1,1,1,1,0},sphere,octant)-(2-std::sqrt(3.)))<1e-13,"interior radial projection");
    require(std::abs(charge_to_spherical_panel_distance({-1,1,1,1,0},sphere,octant)-std::sqrt(1+2*std::pow(std::sqrt(2.)-1,2)))<1e-13,"great-circle edge minimum");
    require(std::abs(charge_to_spherical_panel_distance({3,-1,-1,1,0},sphere,octant)-std::sqrt(3.))<1e-13,"vertex minimum");
    std::swap(octant.v1y,octant.v2y);std::swap(octant.v1z,octant.v2z);
  }
  auto measure=[&](const std::vector<Atom>&a){auto p=build_surface(a,o,nullptr,&executor);return diagnose_charge_resolution(a,p,o);};
  auto centered=measure({{0,0,0,1,1.7}});
  require(centered.min_ratio>3 && centered.min_ratio<3.5 && !centered.warning,"centered default separation");
  require(std::abs(centered.min_distance_over_atom_radius-1)<1e-12 && centered.surface_inside_atom_count==0,"centered clearance");
  auto scaled=measure({{0,0,0,-2,3.4}});
  require(std::abs(centered.min_ratio-scaled.min_ratio)<1e-12,"uniform scaling invariance");
  auto shifted=measure({{100,-50,25,1,1.7}});
  require(std::abs(centered.min_ratio-shifted.min_ratio)<1e-10,"translation invariance");
  std::vector<Atom> offset={{0,0,0,0,1.7},{1.36,0,0,1,.017}};
  auto coarse=measure(offset);require(coarse.warning && coarse.min_ratio<1,"near-boundary coarse warning");
  o.surface_point_count=4096;auto fine=measure(offset);
  require(!fine.warning && fine.min_ratio>2.4,"fine near-boundary separation");
  o.surface_point_count=96;
  std::vector<Atom> overlap={{0,0,0,1,1.7},{2.5,0,0,-1,1.7}};
  auto overlap_metric=measure(overlap);
  require(overlap_metric.min_distance_over_atom_radius>=1-1e-12 && overlap_metric.surface_inside_atom_count==0,"clipped patches respect ideal-union clearance");
  overlap.insert(overlap.end(),{{.7,.4,.3,.3,0},{3.1,-.5,.2,-.7,0},{-1.7,0,0,.2,0},{10,2,3,1,0}});
  auto panels=build_surface(overlap,o,nullptr,&executor);auto accelerated=diagnose_charge_resolution(overlap,panels,o);
  double exhaustive=INFINITY,own=INFINITY;
  for(const auto &a:overlap){
    double dmin=INFINITY;
    for(const auto &p:panels){double d=charge_to_spherical_panel_distance(a,overlap[p.atom_index],p);exhaustive=std::min(exhaustive,d/edge(p));dmin=std::min(dmin,d);}
    if(a.radius>0)own=std::min(own,dmin/a.radius);
  }
  require(std::abs(accelerated.min_ratio-exhaustive)<1e-12,"BVH matches exhaustive normalized distance");
  require(std::abs(accelerated.min_distance_over_atom_radius-own)<1e-12,"BVH matches exhaustive clearance");
  for(std::size_t i=0;i<overlap.size();++i){
    auto one=overlap;for(auto &a:one)a.q=0;one[i].q=1;
    auto indexed=diagnose_charge_resolution(one,panels,o);
    double reference=INFINITY;
    for(const auto &p:panels)reference=std::min(reference,charge_to_spherical_panel_distance(one[i],one[p.atom_index],p)/edge(p));
    require(std::abs(indexed.min_ratio-reference)<1e-12,"per-charge BVH agrees with exhaustive search");
  }
  auto neutral=measure({{0,0,0,0,1.7}});require(neutral.status=="no_charges" && !neutral.warning,"neutral mesh not falsely certified");
  auto invalid=measure({{0,0,0,std::numeric_limits<double>::quiet_NaN(),1.7}});require(invalid.status=="invalid_charge","invalid charge status");
  Options unsupported=o;unsupported.surface_mode="sas";
  require(diagnose_charge_resolution(overlap,panels,unsupported).status=="unsupported_surface","unsupported surface status");
  std::vector<Atom> many;for(int i=0;i<8;++i)for(int j=0;j<8;++j)many.push_back({10.*i,10.*j,0,1,1.7});
  o.surface_point_count=80;auto mesh=build_surface(many,o,nullptr,&executor);auto metric=diagnose_charge_resolution(many,mesh,o);
  require(metric.patch_evaluations<many.size()*mesh.size()/10,"spatial pruning avoids all-pairs scan");
  require(!metric.warning && metric.surface_inside_atom_count==0,"separated centered spheres");
  std::cout<<"PASS "<<checks<<" charge-resolution checks\n";print(centered,320);print(coarse,320);print(fine,5120);print(overlap_metric,panels.size());print(metric,mesh.size());
  return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
