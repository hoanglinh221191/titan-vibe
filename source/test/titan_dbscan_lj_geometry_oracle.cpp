#include "titan_exact_periodic_box.h"
#include "titan_llps_cluster.h"
#include "titan_lj.h"
#include "titan_executor.h"
#include <iostream>
#include <random>
#include <set>
#include <numeric>
void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}
void close(double a,double b,double tol=1e-7){require(std::isfinite(a)&&std::abs(a-b)<=tol*std::max(1.,std::abs(b)),"numeric mismatch");}
int main(){try{
 std::mt19937 rng(892);std::uniform_real_distribution<double> uniform(-1,2);
 int cases=0,pairs=0;
 for(auto angles:{std::array<double,3>{90,90,90},{90,90,120},{75,85,55},{90,90,25}}) {
  box b(20.,24.,28.,double(angles[0]),double(angles[1]),double(angles[2]));TitanExactPeriodicBox cell(b);
  const int n=256;const double eps=3.7;MOL frame;frame.pos.resize(n);frame.Box=b;
  std::vector<LLPSNode> nodes;
  for(int i=0;i<n;++i){auto p=cell.metric.scaled_to_cartesian(uniform(rng),uniform(rng),uniform(rng));frame.pos.x[i]=p[0];frame.pos.y[i]=p[1];frame.pos.z[i]=p[2];nodes.push_back({i+1,1,i,p[0],p[1],p[2],0,false});}
  std::vector<std::vector<int>> near(n);std::set<std::pair<int,int>> expected;
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j) {
   double best=1e100;double dx=frame.pos.x[i]-frame.pos.x[j],dy=frame.pos.y[i]-frame.pos.y[j],dz=frame.pos.z[i]-frame.pos.z[j];
   for(int x=-6;x<=6;++x)for(int y=-6;y<=6;++y)for(int z=-6;z<=6;++z){double qx=dx-cell.m[0][0]*x-cell.m[1][0]*y-cell.m[2][0]*z,qy=dy-cell.m[1][1]*y-cell.m[2][1]*z,qz=dz-cell.m[2][2]*z;best=std::min(best,qx*qx+qy*qy+qz*qz);}
   cell.apply_mic(dx,dy,dz);close(dx*dx+dy*dy+dz*dz,best);++pairs;
   if(best<=eps*eps){expected.insert({i,j});near[i].push_back(j);near[j].push_back(i);}
  }
  std::vector<int> root(n);std::iota(root.begin(),root.end(),0);
  auto find=[&](int i){while(root[i]!=i)i=root[i];return i;};
  for(auto [i,j]:expected)if(near[i].size()+1>=4&&near[j].size()+1>=4){int x=find(i),y=find(j);root[std::max(x,y)]=std::min(x,y);}
  std::set<int> roots;for(int i=0;i<n;++i)if(near[i].size()+1>=4)roots.insert(find(i));
  std::vector<int> labels(n,-1);for(int i=0;i<n;++i){int r=n;if(near[i].size()+1>=4)r=find(i);else for(int j:near[i])if(near[j].size()+1>=4)r=std::min(r,find(j));if(r<n)labels[i]=1+std::distance(roots.begin(),roots.find(r));}
  for(int threads:{1,4}) {
   TitanExecutor executor(threads);TitanLLPSCluster db;db.eps_=eps;db.eps_sq_=eps*eps;db.min_pts_=4;auto actual=nodes;db.run_dbscan(executor,actual,b);
   for(int i=0;i<n;++i){require(db.neighbor_count_[i]==int(near[i].size()),"DBSCAN neighbor multiplicity");require(actual[i].is_core==(near[i].size()+1>=4),"DBSCAN core classification");require(actual[i].cluster_id==labels[i],"DBSCAN cluster/noise oracle");}
   ++cases;
  }
  TitanLJAnalyzer lj;lj.cutoff_a_=eps;lj.prepare_atom_buffers(n);lj.reset_union_find(n);lj.build_contact_graph(frame,n,true);
  std::set<std::pair<int,int>> actual;for(auto [i,j]:lj.edges_)actual.insert(std::minmax(i,j));require(actual==expected&&actual.size()==lj.edges_.size(),"LJ edge coverage/multiplicity");
  for(int i=0;i<n;++i)require(lj.atom_degree_[i]==int(near[i].size()),"LJ degree");
 }
 box b(100.,100.,100.,80.,90.,65.);TitanExactPeriodicBox cell(b);std::vector<std::array<double,3>> points;
 for(int i=0;i<8;++i){auto shift=cell.metric.scaled_to_cartesian(i-4,2-i,i%2);points.push_back({i*4.+shift[0],shift[1],shift[2]});}
 std::vector<std::array<double,3>> lifted;require(titan_lift_radius_graph(cell,points,4.5,lifted),"finite long chain unwrap");
 for(int i=1;i<8;++i)close(lifted[i][0]-lifted[0][0],4.*i);
 box ringbox(20.,20.,20.,90.,90.,90.);TitanExactPeriodicBox ring(ringbox);
 require(!titan_lift_radius_graph(ring,{{0,0,0},{5,0,0},{10,0,0},{15,0,0}},5.1,lifted),"winding cluster must abstain");
 Multi_MOL system;system.m_mol.resize(1);auto& frame=system.m_mol[0];frame.pos.resize(640);frame.Box=box(60.,60.,60.,75.,80.,65.);TitanExactPeriodicBox large(frame.Box);
 for(int i=0;i<640;++i){auto shift=large.metric.scaled_to_cartesian(i%3-1,i%5-2,i%7-3);frame.pos.x[i]=i%8+shift[0];frame.pos.y[i]=(i/8)%10+shift[1];frame.pos.z[i]=i/80+shift[2];system.info.emplace_back(i+1,"C","CA","A","ALA");system.info.back().set_chainnum(i);}
 double previous_q6=0;
 for(int threads:{1,4}) {
  TitanExecutor executor(threads);TitanLJAnalyzer lj;lj.atom_count_=640;lj.initialized_=true;lj.prepare_atom_buffers(640);LjArgs args;args.flag_lj=true;args.lj_cluster_cutoff=1.1;args.lj_rho_bin=1;args.lj_rho_max=10;lj.configure(args);lj.process_frame(system,0,executor);
  require(lj.frames_.size()==1&&lj.frames_[0].largest_cluster_size==640,"LJ full pipeline cluster");close(lj.frames_[0].largest_cluster_rg_a,std::sqrt(225./12));require(lj.frames_[0].mean_dense_q6>0,"LJ local order not calculated");
  if(threads==1)previous_q6=lj.frames_[0].mean_dense_q6;else close(lj.frames_[0].mean_dense_q6,previous_q6);
  TitanLLPSCluster db;db.init(1.1,3,"","");MOL empty;std::vector<mol_info> noinfo;db.process_frame(executor,frame,system.info,empty,noinfo,0,0);
  require(db.droplet_shape_rows_.size()==1,"DBSCAN shape row missing");close(db.droplet_shape_rows_[0].rg_angstrom,std::sqrt(225./12));
 }
 std::vector<double> first_weights;
 for(int threads:{1,4}) {
  TitanExecutor executor(threads);TitanLLPSCluster db;db.init(1.,4,"","");std::vector<LLPSNode> nodes;
  for(int i=0;i<384;++i){auto shift=large.metric.scaled_to_cartesian(i%3-1,i%5-2,i%7-3);nodes.push_back({i+1,1,i,shift[0]+(i<192?0:8)+.01*(i%8),shift[1]+.01*(i%9),shift[2]+.01*(i%7),0,false});}
  db.run_dbscan(executor,nodes,frame.Box);std::vector<bool> dense;std::vector<double> weights,occupancy;db.fuzzy_refine(executor,nodes,frame.Box,1,dense,weights,occupancy);
  require(weights.front()>.99&&weights.back()<.01,"fuzzy weights not calculated");
  if(threads==1)first_weights=weights;else for(std::size_t i=0;i<weights.size();++i)close(weights[i],first_weights[i]);
 }
 std::cout<<"PASS DBSCAN/LJ geometry "<<pairs<<" exhaustive lattice pair distances; "<<cases<<" DBSCAN serial/parallel label oracles; LJ complete unique edges/degrees; finite and winding embeddings\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
