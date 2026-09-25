// Streaming CA-only descriptors; no trajectory fit, kinetic estimate or MD.
#include "xdrfile_xtc.h"
#include "pbctopo/lattice_math.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <algorithm>
using Vec=std::array<double,3>;
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("trajectory selection output required");
 std::ifstream selection(argv[2]);std::array<std::vector<int>,2> ca;int chain,index;
 while(selection>>chain>>index){if(chain<0||chain>1||index<0)throw std::runtime_error("Bad CA index");ca[chain].push_back(index);}
 if(ca[0].size()!=42||ca[1].size()!=42)throw std::runtime_error("Expected two 42-residue chains");
 int n=0;std::vector<char> path(argv[1],argv[1]+std::char_traits<char>::length(argv[1])+1);
 if(read_xtc_natoms(path.data(),&n)!=exdrOK||n!=1254)throw std::runtime_error("Expected protein projection of 1254 atoms");
 XDRFILE* f=xdrfile_open(argv[1],"r");if(!f)throw std::runtime_error("Cannot open XTC");
 std::vector<std::array<float,3>> xyz(n);matrix cell;int step;float t,prec;long long frame=0;
 std::ofstream out(argv[3]);if(!out)throw std::runtime_error("Cannot create output");out<<std::setprecision(12)<<"frame,step,time_ps,valid,com_distance_nm,rg_min_nm,rg_max_nm,end_min_nm,end_max_nm,max_ca_step_nm\n";
 for(;;){int status=read_xtc(f,n,&step,&t,cell,reinterpret_cast<rvec*>(xyz.data()),&prec);if(status==exdrENDOFFILE)break;if(status!=exdrOK)throw std::runtime_error("XTC decoding failed");
  // LatticeMetric stores basis vectors in rows, matching the XTC matrix.
  std::array<double,9> mat;for(int i=0;i<3;++i)for(int j=0;j<3;++j)mat[3*i+j]=cell[i][j];
  titan_pbctopo::PbctopoLatticeMetric lattice{std::span<const double,9>(mat)};
  bool valid=lattice.valid();double maxstep=0;std::array<Vec,2> com{};std::array<double,2> rg{},ends{};
  for(int c=0;c<2&&valid;++c){std::vector<Vec> lift(42);for(int k=0;k<3;++k)lift[0][k]=xyz[ca[c][0]][k];
   for(int a=1;a<42;++a){auto &p=xyz[ca[c][a]],&q=xyz[ca[c][a-1]];auto d=lattice.nearest_image_cartesian(p[0]-q[0],p[1]-q[1],p[2]-q[2]);
    double distance=std::sqrt(d.d2);maxstep=std::max(maxstep,distance);if(!d.valid||d.ambiguous||d.node_budget_exhausted||distance<.2||distance>.5){valid=false;break;}
    lift[a]={lift[a-1][0]+d.x,lift[a-1][1]+d.y,lift[a-1][2]+d.z};}
   if(!valid)break;for(auto p:lift)for(int k=0;k<3;++k)com[c][k]+=p[k]/42.;
   for(auto p:lift)for(int k=0;k<3;++k)rg[c]+=(p[k]-com[c][k])*(p[k]-com[c][k])/42.;rg[c]=std::sqrt(rg[c]);
   for(int k=0;k<3;++k)ends[c]+=(lift.back()[k]-lift.front()[k])*(lift.back()[k]-lift.front()[k]);ends[c]=std::sqrt(ends[c]);
  }
  double distance=0;if(valid){auto d=lattice.nearest_image_cartesian(com[1][0]-com[0][0],com[1][1]-com[0][1],com[1][2]-com[0][2]);valid=d.valid&&!d.ambiguous&&!d.node_budget_exhausted;distance=std::sqrt(d.d2);}
  out<<frame++<<','<<step<<','<<t<<','<<valid;
  if(valid)out<<','<<distance<<','<<std::min(rg[0],rg[1])<<','<<std::max(rg[0],rg[1])<<','<<std::min(ends[0],ends[1])<<','<<std::max(ends[0],ends[1])<<','<<maxstep;
  else out<<",nan,nan,nan,nan,nan,"<<maxstep;out<<'\n';
 }
 xdrfile_close(f);if(frame!=21240)throw std::runtime_error("Expected 21240 audited frames");out.close();if(!out)throw std::runtime_error("Output failed");std::cout<<"PASS decoded "<<frame<<" frames\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
