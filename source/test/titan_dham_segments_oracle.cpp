#include "titan_dham_segments.h"
#include <iostream>
using namespace titan_dham_segments;
void require(bool b,const char* message){if(!b)throw std::runtime_error(message);}
void save(const std::string& data){std::ofstream f("fixture.txt");f<<data;}
int main(int argc,char** argv){try {
 if(argc==2) {
  auto m=Metadata::read(argv[1],48);
  const std::size_t expected[]={0,262,3,0};int lags[]={1,2,4,10};
  for(int v=0;v<4;++v) {
   std::size_t total=0;
   for(int k=0;k<48;++k) {
    std::vector<int> states(m.replicas[k].frames.size(),0);std::vector<double> times;
    for(const auto& f:m.replicas[k].frames)times.push_back(f.time);
    auto a=audit_counts(m,k,states,times,lags[v],1);total+=a.decisions[Decision::eligible];
   }
   require(total==expected[v],"real-data gate differs from independent Python eligibility");
   std::cout<<"PASS real metadata lag_ps="<<lags[v]*50<<" eligible="<<total<<'\n';
  }
  return 0;
 }
 const std::string fixture=
 "TITAN_DHAM_SEGMENTS_V1\nreplica 0 300 1\nsegment 0 1 0 300\nsegment 0 2 0 100\n"
 "frame 0 0 0 1 0 1 2 exact\nframe 0 1 50 1 50 NA 2 absent\nframe 0 2 100 1 100 2 2 exact\n"
 "frame 0 3 150 1 150 NA 2 absent\nframe 0 4 200 1 200 3 2 exact\nframe 0 5 250 1 250 NA 2 absent\n"
 "frame 0 6 300 1 300 4 2 exact\nframe 0 7 350 2 0 5 2 exact\nframe 0 8 400 2 50 NA 2 absent\nframe 0 9 450 2 100 6 2 exact\n"
 "exchange 0 1 150\nexchange 0 1 175\n";
 save(fixture);auto m=Metadata::read("fixture.txt",1);std::vector<int> states{0,1,1,0,0,1,1,0,1,1};std::vector<double> times{0,50,100,150,200,250,300,350,400,450};
 auto a=audit_counts(m,0,states,times,2,2);
 require(a.decisions[Decision::eligible]==3&&a.decisions[Decision::exchange]==3&&a.decisions[Decision::segment_boundary]==2,"segment/exchange transition oracle");
 require(a.counts[std::pair{0,1}]==3,"state counts after masks");
 auto b=audit_counts(m,0,states,times,1,2);require(b.decisions[Decision::eligible]==0&&b.decisions[Decision::missing_energy]==6,"missing energy never substituted");
 require(m.transition(0,m.match(0,0),m.match(0,100),50)==Decision::time_gap,"physical lag mismatch");
 require(m.transition(0,m.match(0,100),m.match(0,150),50)==Decision::exchange,"exchange exactly at end");
 require(m.transition(0,m.match(0,150),m.match(0,200),50)==Decision::exchange,"exchange exactly at start");
 int rejected=0;
 try{m.match(0,40);}catch(const std::invalid_argument&){++rejected;}
 try{m.match(0,std::numeric_limits<double>::quiet_NaN());}catch(const std::invalid_argument&){++rejected;}
 for(const auto& extra:{"exchange 0 1 150\n","exchange 0 1 301\n","frame 0 10 500 1 200 7 2 exact\n","frame 0 10 500 2 150 7 2 exact\n","bogus 0\n"}) {
  save(fixture+extra);try{Metadata::read("fixture.txt",1);}catch(const std::invalid_argument&){++rejected;}
 }
 require(rejected==7,"invalid provenance accepted");
 std::cout<<"PASS segment/exchange masks, exact energy holes, closed boundaries, lag mismatch, 7 malformed/missing input rejections\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
