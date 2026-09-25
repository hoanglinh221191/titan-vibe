#include "titan_llps_geometry.h"
#include "titan_neighbor_grid.h"
#include <iostream>
#include <random>
#include <set>
using namespace titan_llps_geometry;
void require(bool yes,const char* msg) {if(!yes)throw std::runtime_error(msg);}
void close(double a,double b,double tol=1e-9) {require(std::isfinite(a)&&std::abs(a-b)<tol*std::max(1.0,std::abs(b)),"numeric oracle mismatch");}
int main(){try {
 std::mt19937 rng(88);std::uniform_real_distribution<double> u(-2,2);
 int checks=0;
 for(const auto b:{std::array<double,9>{8,0,0,0,9,0,0,0,10}, {8,5,2,0,7,3,0,0,9}, {3,2.8,1,0,1.1,.5,0,0,4}}) {
  Cell cell(b);std::vector<Point> points;
  for(int n=0;n<48;++n)points.push_back(cell.cartesian({u(rng),u(rng),u(rng)}));
  auto brute=[&](Point d){double best=1e100;for(int x=-8;x<=8;++x)for(int y=-8;y<=8;++y)for(int z=-8;z<=8;++z){auto q=cell.cartesian({double(x),double(y),double(z)});double v=0;for(int k=0;k<3;++k)v+=std::pow(d[k]-q[k],2);best=std::min(best,v);}return best;};
  std::set<std::pair<int,int>> expected,actual;
  constexpr double cutoff=.9;
  for(int i=0;i<48;++i)for(int j=i+1;j<48;++j){Point d{};for(int k=0;k<3;++k)d[k]=points[i][k]-points[j][k];auto q=cell.delta(d);double d2=q[0]*q[0]+q[1]*q[1]+q[2]*q[2];close(d2,brute(d));++checks;if(d2<cutoff*cutoff)expected.insert({i,j});}
  titan_neighbor::HalfShellGrid grid;titan_neighbor::OrthorhombicGridSpec spec;
  spec.periodic=true;spec.cutoff=cutoff;spec.len_x=cell.heights[0];spec.len_y=cell.heights[1];spec.len_z=cell.heights[2];
  grid.for_each_pair(points.size(),spec,[&](std::size_t i,double& x,double& y,double& z){x=points[i][0];y=points[i][1];z=points[i][2];cell.grid_coordinate(x,y,z);},[](std::size_t){return true;},[&](std::size_t i,std::size_t j){Point d{};for(int k=0;k<3;++k)d[k]=points[i][k]-points[j][k];auto q=cell.delta(d);if(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]<cutoff*cutoff)require(actual.insert(std::minmax(int(i),int(j))).second,"duplicate contact");});
  require(actual==expected,"triclinic grid missed contacts");
 }
 Cell cell({20,5,2,0,20,3,0,0,20});
 PhaseConfig drop{"droplet",2,1,2,3};
 std::vector<Point> all{{.2,0,0},{-.2,0,0},{2.5,0,0},{7,0,0}},gcc{all[0],all[1]};
 auto a=measure(cell,drop,all,gcc,Point{0,0,0});require(a.status=="ok"&&a.dense_count==2&&a.dilute_count==1,"droplet counts");
 close(a.dense,2/(4*std::acos(-1.)/3));close(a.dilute,1/(4*std::acos(-1.)/3*19));
 auto translation=cell.cartesian({3,-2,1});for(auto& p:all)for(int k=0;k<3;++k)p[k]+=translation[k];
 auto shifted=measure(cell,drop,all,gcc,Point{0,0,0});close(a.dense,shifted.dense);close(a.dilute,shifted.dilute);
 for(int axis=0;axis<3;++axis) {
  PhaseConfig slab{"slab",axis,1,2,3};all.clear();gcc.clear();
  for(double d:{-.2,.2,2.5,7.}){Point f{.31,.4,.6};f[axis]=d/cell.heights[axis];all.push_back(cell.cartesian(f));}
  gcc={all[0],all[1]};auto r=measure(cell,slab,all,gcc,{});
  require(r.status=="ok"&&r.dense_count==2&&r.dilute_count==1,"slab counts");
  close(r.dense,2/(cell.volume/cell.heights[axis]*2));close(r.dilute,1/(cell.volume/cell.heights[axis]*2));
  for(int i=0;i<4;++i) {auto d=cell.cartesian({double(i-2),double(2*i),double(-i)});for(int k=0;k<3;++k)all[i][k]+=d[k];}
  gcc={all[0],all[1]};auto t=measure(cell,slab,all,gcc,{});close(r.dense,t.dense);close(r.dilute,t.dilute);
 }
 drop.dilute_outer=11;require(measure(cell,drop,all,gcc,Point{0,0,0}).status=="window_exceeds_cell","invalid window accepted");
 bool rejected=false;try{PhaseConfig{"slab",2,2,1,3}.validate();}catch(...){rejected=true;}require(rejected,"overlapping windows accepted");
 rejected=false;try{Cell::from_box(box{},-1,2,3);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"negative fallback cell accepted");
 require(measure(cell,PhaseConfig{"droplet",2,1,2,3},all,gcc,{}).status=="gcc_embedding_ambiguous","ambiguous droplet accepted");
 std::cout<<"PASS LLPS geometry: "<<checks<<" exhaustive lattice distance comparisons; contact completeness/multiplicity; droplet/slab a,b,c counts and volumes; independent image invariance; invalid/ambiguous domains\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
