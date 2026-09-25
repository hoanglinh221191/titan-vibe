#include "xdrfile_xtc.h"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
struct Frame {std::vector<std::array<float,3>> xyz;int step;float time,precision;matrix cell;};
Frame read(const std::string& path,long long offset){int n=0;std::vector<char> file(path.begin(),path.end());file.push_back(0);if(read_xtc_natoms(file.data(),&n)!=exdrOK||n<1)throw std::runtime_error("Cannot read atom count");XDRFILE* input=xdrfile_open(path.c_str(),"r");if(!input)throw std::runtime_error("Cannot open XTC");Frame f;f.xyz.resize(n);if(xdr_seek(input,offset,SEEK_SET)!=exdrOK){xdrfile_close(input);throw std::runtime_error("Cannot seek XTC");}int status=read_xtc(input,n,&f.step,&f.time,f.cell,reinterpret_cast<rvec*>(f.xyz.data()),&f.precision);xdrfile_close(input);if(status!=exdrOK)throw std::runtime_error("Cannot decode XTC");return f;}
int main(int argc,char**argv){try{if(argc!=2)return 2;std::ifstream in(argv[1]);std::string full,protein;long long a,b;int cases=0;double maximum=0;
 while(in>>std::quoted(full)>>a>>std::quoted(protein)>>b){auto x=read(full,a),y=read(protein,b);if(y.xyz.size()!=1254||x.xyz.size()!=29037||x.step!=y.step||x.time!=y.time)throw std::runtime_error("Projection metadata mismatch");for(std::size_t i=0;i<y.xyz.size();++i)for(int k=0;k<3;++k)maximum=std::max(maximum,std::abs(double(x.xyz[i][k])-y.xyz[i][k]));++cases;}
 std::cout<<std::setprecision(12)<<"projection_cases="<<cases<<" max_coordinate_difference_nm="<<maximum<<" tolerance_nm=0.002\n";
 if(!cases||maximum>.002)return 1;std::cout<<"PASS decoded first-1254-atom projection identity\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
