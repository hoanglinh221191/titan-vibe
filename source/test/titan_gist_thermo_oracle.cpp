#include "titan_gist_entropy.h"
#ifdef TITAN_TEST_ENGINE
#include "titan_gist.h"
#endif
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
using Q=titan_gist::Quaternion;
void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
std::vector<Q> random_q(size_t n){std::mt19937_64 rng(12349);std::normal_distribution<double> normal;std::vector<Q> q(n);for(auto& p:q){double norm=0;for(auto& a:p){a=normal(rng);norm+=a*a;}for(auto& a:p)a/=std::sqrt(norm);}return q;}
double brute_angle(const std::vector<Q>& q,size_t i){double best=0;for(size_t j=0;j<q.size();++j)if(j!=i){double dot=0;for(int k=0;k<4;++k)dot+=q[i][k]*q[j][k];best=std::max(best,std::abs(dot));}return 2*std::acos(std::min(1.,best));}
int main(int argc,char**argv){try{
 const double pi=std::acos(-1.);require(std::abs(titan_gist::ball_probability(pi)-1)<1e-14,"Haar ball pi normalization");
 require(std::abs(titan_gist::ball_probability(1e-6)/(1e-18/(6*pi))-1)<1e-12,"small angle Haar volume");
 auto q=random_q(4096);titan_gist::QuaternionTree tree(q);for(size_t i=0;i<100;++i)require(std::abs(tree.nearest_angle(i)-brute_angle(q,i))<1e-10,"antipodal KD versus dot-product brute force");
 auto entropy=titan_gist::orientational_entropy(q);require(std::isfinite(entropy.s_over_R)&&std::abs(entropy.s_over_R)<.08,"uniform SO3 entropy");
 auto flipped=q;for(size_t i=0;i<flipped.size();i+=2)for(auto& a:flipped[i])a=-a;require(std::abs(titan_gist::orientational_entropy(flipped).s_over_R-entropy.s_over_R)<1e-12,"quaternion sign invariance");
 require(!std::isfinite(titan_gist::orientational_entropy({q[0]}).s_over_R),"single orientation NaN");auto duplicate=titan_gist::orientational_entropy({q[0],q[0]});require(duplicate.zero_neighbors==2&&!std::isfinite(duplicate.s_over_R),"duplicate orientation NaN");
 Q rotation{.5,.5,.5,.5};Eigen::Quaterniond r(rotation[0],rotation[1],rotation[2],rotation[3]);auto rotated=q;for(auto& p:rotated){auto v=r*Eigen::Quaterniond(p[0],p[1],p[2],p[3]);p={v.w(),v.x(),v.y(),v.z()};}require(std::abs(titan_gist::orientational_entropy(rotated).s_over_R-entropy.s_over_R)<1e-10,"global rotation entropy invariance");
 std::cout<<std::setprecision(15)<<"PASS SO3 analytic Haar, exact KD, sign/rotation, duplicate gates; uniform S/R="<<entropy.s_over_R<<'\n';
#ifdef TITAN_TEST_ENGINE
 std::vector<mol_info> ri{mol_info(1,"C","CA","A","ALA"),mol_info(2,"C","CA","A","ALA"),mol_info(3,"C","CA","A","ALA")},li{mol_info(1,"C","C","B","LIG")};Multi_MOL_parameters rp(3),lp(1);MOL pro,lig;pro.pos.resize(3);lig.pos.resize(1);pro.pos.x[0]=-1;pro.pos.y[0]=-1;pro.pos.x[1]=1;pro.pos.y[1]=-1;pro.pos.y[2]=1;
 for(double gamma:{35.,70.,90.,120.}){
  TitanGIST gist;gist.init(lig,6,.5,300,9,1,"TIP3P",rp,ri,lp,li);gist.periodic_box_.emplace(box(60.,65.,70.,90.,90.,double(gamma)));gist.wat_buf_.resize(64*3);std::vector<int> voxels(64,0);for(int i=0;i<64;++i){double x=(i%4)*3.2,y=((i/4)%4)*3.2,z=(i/16)*3.2;gist.wat_buf_[3*i]={x,y,z};gist.wat_buf_[3*i+1]={x+.7,y+.6,z};gist.wat_buf_[3*i+2]={x-.7,y+.6,z};if(i%3==0)voxels[i]=-1;}
  auto fast=gist.water_water_energies(voxels),slow=gist.water_water_energies(voxels,false);for(size_t i=0;i<fast.size();++i)require(std::isfinite(fast[i])&&std::abs(fast[i]-slow[i])<1e-9,"periodic water pair coverage/half-allocation");
  for(size_t i=0;i<gist.wat_buf_.size();++i){auto d=gist.periodic_box_->metric.scaled_to_cartesian(double(int(i%3)-1),double(int(i%5)-2),0);for(int k=0;k<3;++k)gist.wat_buf_[i][k]+=d[k];}auto shifted=gist.water_water_energies(voxels);for(size_t i=0;i<fast.size();++i)require(std::abs(fast[i]-shifted[i])<1e-8,"water pair arbitrary atom-image invariance");
 }
 // Two-water analytic case; outside-grid endpoint retains its contribution
 // to the on-grid water while receiving no on-grid allocation of its own.
 TitanGIST pair;pair.init(lig,6,.5,300,9,1,"TIP3P",rp,ri,lp,li);pair.wat_buf_={{{0,0,0}},{{.7,.6,0}},{{-.7,.6,0}},{{3,0,0}},{{3.7,.6,0}},{{2.3,.6,0}}};double e=4*.636*(std::pow(3.1507/3,12)-std::pow(3.1507/3,6));double charge[3]={-.834,.417,.417};for(int i=0;i<3;++i)for(int j=0;j<3;++j){double d2=0;for(int k=0;k<3;++k)d2+=std::pow(pair.wat_buf_[i][k]-pair.wat_buf_[j+3][k],2);e+=1389.3546*charge[i]*charge[j]/std::sqrt(d2);}auto v=pair.water_water_energies({0,-1});require(std::abs(v[0]-.5*e)<1e-10&&v[1]==0,"off-grid partner analytic Eww");
 if(argc==3){std::ifstream input(argv[1]);int nf,nw;input>>nf>>nw;require(input&&nf>1&&nw>0,"reference fixture header");MOL wat;wat.pos.resize(nw*3);box box50(50.,50.,50.,90.,90.,90.);pro.Box=box50;lig.Box=box50;wat.Box=box50;TitanGIST run;run.init(lig,6,.5,300,9,4,"TIP3P",rp,ri,lp,li);
  for(int f=0;f<nf;++f){for(size_t i=0;i<wat.pos.size();++i)input>>wat.pos.x[i]>>wat.pos.y[i]>>wat.pos.z[i];require(bool(input),"reference fixture coordinates");run.process_frame(pro,lig,wat,ri,Eigen::MatrixXd{},Eigen::Matrix3d::Identity(),f==0,false);}
  run.export_dx_maps(argv[2],nf,2);std::ofstream out(std::string(argv[2])+"_numerical.csv");out<<std::setprecision(17)<<"voxel,x,y,z,population,Eww_kj,minus_TSorient_kj\n";
  for(int i=0;i<run.n_voxels();++i)if(run.grid_[i].count){int ix=i/(run.ny_*run.nz_),iy=(i/run.nz_)%run.ny_,iz=i%run.nz_;auto entropy=titan_gist::orientational_entropy(run.orientations_.at(i));out<<i<<','<<run.origin_[0]+(ix+.5)*run.voxel_size_<<','<<run.origin_[1]+(iy+.5)*run.voxel_size_<<','<<run.origin_[2]+(iz+.5)*run.voxel_size_<<','<<run.grid_[i].count<<','<<run.grid_[i].energy_ww/run.grid_[i].count<<','<<-.008314462*300*entropy.s_over_R<<'\n';}
 }
 std::cout<<"PASS Eww analytic pair, off-grid partners, triclinic coverage/image gates\n";
#endif
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
