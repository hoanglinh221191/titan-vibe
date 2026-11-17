#include "titan_llps_geometry.h"
#include "titan_exact_periodic_box.h"
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
 // STR-01 (1): Cell::from_box must reproduce the box lattice. Rows of the
 // cell matrix are the lattice vectors (the PbctopoLatticeMetric and
 // BoxMatrix convention); from_box used to pass the box vectors as columns,
 // which is the transposed lattice for a triclinic box. Every box vector
 // must be a period, nearest-image distances must equal those of the exact
 // box metric (TitanExactPeriodicBox, A), heights and volume as well.
 int box_checks=0;
 for(const auto angles:{std::array<double,6>{150,60,60,70,60,75},{70,70,70,60,60,90},{70,70,70,70.5288,109.4712,70.5288},{80,95,60,90,90,90}}) {
  box b(angles[0],angles[1],angles[2],angles[3],angles[4],angles[5]);
  const auto from=Cell::from_box(b,0,0,0);TitanExactPeriodicBox exact(b);
  double v[3][3];b.get_v1(v[0]);b.get_v2(v[1]);b.get_v3(v[2]);
  for(int k=0;k<3;++k){auto q=from.delta({.1*v[k][0],.1*v[k][1],.1*v[k][2]});require(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]<1e-20,"box vector is not a period of Cell::from_box");close(from.heights[k],.1*exact.heights[k],1e-12);}
  close(from.volume,1e-3*b.volume(),1e-12);
  for(int n=0;n<200;++n){
   const Point d{u(rng)*angles[0]*.1,u(rng)*angles[1]*.1,u(rng)*angles[2]*.1};
   auto q=from.delta(d);double x=10*d[0],y=10*d[1],z=10*d[2];exact.apply_mic(x,y,z);
   close(q[0]*q[0]+q[1]*q[1]+q[2]*q[2],.01*(x*x+y*y+z*z),1e-10);++box_checks;
  }
 }
 // STR-01 (2): slab normal and lateral (in-plane) squared displacement of
 // the -llps-dynamics lateral MSD. Axis = lattice vector with the largest
 // plane spacing (near ties: c, then b); normal = unit normal of the planes
 // of the two other lattice vectors. Cells (rows a, b, c): orthorhombic
 // along z and x (the old Cartesian rule, to rounding), triclinic slab along
 // a (normal not x: the old rule fails), triclinic slab along c (GROMACS
 // lower-triangular box: normal exactly z), rhombic dodecahedron (heights
 // 0.8165 d along a and b, 0.7071 d along c: b by the tie rule).
 int lateral_checks=0;
 const double r2=std::sqrt(2.)/2;
 const std::array<std::array<double,9>,5> lateral_cells{{{8,0,0,0,9,0,0,0,30},{30,0,0,0,9,0,0,0,8},{30,0,0,5,9,0,2,3,8},{8,0,0,5,7,0,2,3,30},{10,0,0,0,10,0,5,5,10*r2}}};
 const int expected_axis[5]={2,0,0,2,1};
 for(int c=0;c<5;++c) {
  Cell lattice(lateral_cells[c]);const auto sn=slab_normal(lattice);
  require(sn.axis==expected_axis[c],"slab axis = largest cell height");
  for(int k=0;k<3;++k)if(k!=sn.axis)require(lattice.heights[k]<=lattice.heights[sn.axis]*(1+1e-12),"slab axis not the largest height");
  close(sn.normal[0]*sn.normal[0]+sn.normal[1]*sn.normal[1]+sn.normal[2]*sn.normal[2],1.,1e-14);
  Point in_plane{};
  for(int k=0;k<3;++k) if(k!=sn.axis) {
   const Point v{lateral_cells[c][3*k],lateral_cells[c][3*k+1],lateral_cells[c][3*k+2]};
   close(1+v[0]*sn.normal[0]+v[1]*sn.normal[1]+v[2]*sn.normal[2],1.,1e-13); // normal is perpendicular to the plane vectors
   for(int j=0;j<3;++j)in_plane[j]+=(k==0?.3:k==1?-.7:.45)*v[j];
  }
  close(lateral_sq(in_plane,sn.normal),in_plane[0]*in_plane[0]+in_plane[1]*in_plane[1]+in_plane[2]*in_plane[2],1e-12);
  close(1+lateral_sq(sn.normal,sn.normal),1.,1e-14);
  for(int n=0;n<64;++n) {
   const Point d{u(rng),u(rng),u(rng)};const double along=d[0]*sn.normal[0]+d[1]*sn.normal[1]+d[2]*sn.normal[2];
   close(lateral_sq(d,sn.normal),d[0]*d[0]+d[1]*d[1]+d[2]*d[2]-along*along,1e-12);
   const double cartesian=sn.axis==0?d[1]*d[1]+d[2]*d[2]:sn.axis==1?d[0]*d[0]+d[2]*d[2]:d[0]*d[0]+d[1]*d[1];
   if(c<2||c==3)close(lateral_sq(d,sn.normal),cartesian,1e-15);
   ++lateral_checks;
  }
  if(c==2){const Point d{.3,.4,.5};require(std::abs(lateral_sq(d,sn.normal)-(.4*.4+.5*.5))>1e-3,"triclinic slab: the old Cartesian projection must differ");}
 }
 std::cout<<"PASS LLPS geometry: "<<checks<<" exhaustive lattice distance comparisons; contact completeness/multiplicity; droplet/slab a,b,c counts and volumes; independent image invariance; invalid/ambiguous domains; "<<box_checks<<" from_box/exact-box image distances (triclinic, dodecahedron, octahedron); "<<lateral_checks<<" lateral-MSD projections (orthorhombic, triclinic slab along a and c, dodecahedron)\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
