#ifndef TITAN_GIST_IMPLEMENTATION
#define TITAN_GIST_IMPLEMENTATION "titan_gist.cpp"
#endif
#include TITAN_GIST_IMPLEMENTATION
#include <Eigen/Geometry>
#include <filesystem>
#include <iostream>
void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void close(double a,double b,const char* m){require(std::isfinite(a)&&std::abs(a-b)<1e-7*std::max(1.,std::abs(b)),m);}
MOL transform(const MOL& src,const Eigen::Matrix3d& R,const Eigen::Vector3d& shift) {
 MOL out=src;
 for(size_t i=0;i<out.pos.size();++i){Eigen::Vector3d p=R*Eigen::Vector3d(src.pos.x[i],src.pos.y[i],src.pos.z[i])+shift;out.pos.x[i]=p[0];out.pos.y[i]=p[1];out.pos.z[i]=p[2];}
 return out;
}
double energy_sum(const TitanGIST& g){double e=0;for(const auto& v:g.grid_)e+=v.energy_sw;return e;}
int main(){try {
 std::vector<mol_info> ri{mol_info(1,"C","CA","A","ALA"),mol_info(2,"C","CA","A","ALA"),mol_info(3,"C","CA","A","ALA"),mol_info(4,"C","X","A","ALA")};
 std::vector<mol_info> li{mol_info(1,"C","C","B","LIG")};Multi_MOL_parameters rp(4),lp(1);rp.charge[3]=1;rp.sigma[3]=.3;rp.epsilon[3]=.25;
 for(double gamma:{25.,65.,90.,120.}) {
  box Box(100.,120.,140.,90.,90.,double(gamma));TitanExactPeriodicBox cell(Box);
  MOL pro,lig,wat;pro.pos.resize(4);lig.pos.resize(1);wat.pos.resize(108);pro.Box=Box;lig.Box=Box;wat.Box=Box;
  pro.pos.x[0]=-1;pro.pos.y[0]=-1;pro.pos.x[1]=1;pro.pos.y[1]=-1;pro.pos.y[2]=1;
  // Distinct waters away from voxel faces: exact-face occupancy is not stable
  // under floating-point rigid transforms and is tested separately by bounds.
  for(int i=0;i<36;++i){int w=i*3;const double x=-4.09+2.8*(i%4),y=-2.67+2.8*((i/4)%3),z=-2.63+2.8*(i/12);wat.pos.x[w]=x;wat.pos.y[w]=y;wat.pos.z[w]=z;wat.pos.x[w+1]=x+.7;wat.pos.y[w+1]=y+.6;wat.pos.z[w+1]=z;wat.pos.x[w+2]=x-.7;wat.pos.y[w+2]=y+.6;wat.pos.z[w+2]=z;}
  TitanGIST base,moved;base.init(lig,6,.5,300,9,1,"TIP3P",rp,ri,lp,li);moved.init(lig,6,.5,300,9,4,"TIP3P",rp,ri,lp,li);
  auto frame=[&](TitanGIST& g,const MOL& p,const MOL& l,const MOL& w,bool first){g.process_frame(p,l,w,ri,Eigen::MatrixXd{},Eigen::Matrix3d::Identity(),first,false);};
  // Initialization must work at the first selected frame, even if its source
  // frame number is nonzero and is_frame_zero is false.
  frame(base,pro,lig,wat,false);frame(moved,pro,lig,wat,false);frame(base,pro,lig,wat,false);
  auto R=Eigen::AngleAxisd(.63,Eigen::Vector3d::UnitZ()).toRotationMatrix();Eigen::Vector3d shift(6,2,3);
  auto p2=transform(pro,R,shift),l2=transform(lig,R,shift),w2=transform(wat,R,shift);
  for(size_t i=0;i<w2.pos.size();++i){auto image=cell.metric.scaled_to_cartesian(double(int(i%5)-2),double(int(i%7)-3),double(int(i%3)-1));w2.pos.x[i]+=image[0];w2.pos.y[i]+=image[1];w2.pos.z[i]+=image[2];}
  frame(moved,p2,l2,w2,true); // a later file's local frame 0 must not reset reference
  require(base.processed_frames_==2&&moved.processed_frames_==2,"processed frame tracking");
  for(size_t i=0;i<base.grid_.size();++i){require(base.grid_[i].count==moved.grid_[i].count,"rigid image/rotation occupancy invariance");close(base.grid_[i].energy_sw,moved.grid_[i].energy_sw,"rigid image/rotation energy invariance");close(base.grid_[i].dipole_y,moved.grid_[i].dipole_y,"whole-water dipole invariance");}
  require(base.voxel_idx(base.origin_[0]-1e-6,base.origin_[1]+.1,base.origin_[2]+.1)==-1,"lower-grid outside coordinate counted");
  require(base.voxel_idx(std::numeric_limits<double>::infinity(),0,0)==-1,"nonfinite voxel accepted");
  // A non-alignment atom moves; compare the new energy with a direct pair
  // formula instead of a second copy of the production evaluator.
  TitanGIST deform;deform.init(lig,6,.5,300,9,1,"TIP3P",rp,ri,lp,li);frame(deform,pro,lig,wat,true);double before=energy_sum(deform);pro.pos.x[3]=1;frame(deform,pro,lig,wat,false);
  auto dist=[](double x,double y,double z){return std::sqrt((x-1)*(x-1)+y*y+z*z);};
  double expected=0;
  for(int w=0;w<108;w+=3){double o=dist(wat.pos.x[w],wat.pos.y[w],wat.pos.z[w]),h1=dist(wat.pos.x[w+1],wat.pos.y[w+1],wat.pos.z[w+1]),h2=dist(wat.pos.x[w+2],wat.pos.y[w+2],wat.pos.z[w+2]);double sr6=std::pow(.5*(3.1507+3.)/o,6);expected+=1389.3546*(-.834/o+.417/h1+.417/h2)+4*std::sqrt(.636*.25)*(sr6*sr6-sr6);}
  close(energy_sum(deform)-before,expected,"current receptor energy (not frozen frame zero)");
  if(gamma==90.) {
   moved.export_dx_maps("gist_checked",999,1);
   require(std::filesystem::exists("gist_checked_E_sw.dx")&&!std::filesystem::exists("gist_checked_dG.dx")&&!std::filesystem::exists("gist_checked_TdS_orient.dx"),"unsupported thermodynamic maps emitted");
   std::ifstream f("gist_checked_observables.json");std::string text((std::istreambuf_iterator<char>(f)),{});
   require(text.find("\"processed_frames\": 2")!=std::string::npos&&text.find("\"free_energy_available\": false")!=std::string::npos,"observable status mismatch");
   int rejected=0;
   try{TitanGIST g;g.init(lig,1,.5,300,9,1,"unknown",rp,ri,lp,li);}catch(const std::invalid_argument&){++rejected;}
   try{TitanGIST g;g.init(lig,1e200,.5,300,9,1,"TIP3P",rp,ri,lp,li);}catch(const std::invalid_argument&){++rejected;}
   auto bad=wat;bad.pos.resize(4);try{frame(moved,pro,lig,bad,false);}catch(const std::invalid_argument&){++rejected;}
   bad=wat;bad.pos.x[1]=bad.pos.x[0]+4;try{frame(moved,pro,lig,bad,false);}catch(const std::invalid_argument&){++rejected;}
   auto small=pro;small.Box=box(20.,20.,20.,90.,90.,90.);try{frame(moved,small,lig,wat,false);}catch(const std::invalid_argument&){++rejected;}
   require(rejected==5,"invalid water/model/domain inputs accepted");
  }
 }
 std::cout<<"PASS GIST 4 lattices, image/rotation/worker invariance, current receptor analytic energy, voxel bounds, unavailable outputs and input rejections\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
