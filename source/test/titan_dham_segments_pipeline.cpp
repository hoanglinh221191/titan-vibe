#include "postproc/engine_postproc.h"
#include "titan_msm.h"
#include "titan_executor.h"
#include "titan_validation_core.h"
#include <filesystem>
#include <iostream>
void require(bool b,const char* m){if(!b)throw std::runtime_error(m);}
int main(){try {
 std::ofstream("temps.txt")<<"300\n";
 std::ofstream f("segments.txt");f<<"TITAN_DHAM_SEGMENTS_V1\nreplica 0 300 1\nsegment 0 1 0 300\nsegment 0 2 0 100\n";
 const int states[]={0,1,1,0,0,1,1,0,1,1};
 for(int i=0;i<10;++i)f<<"frame 0 "<<i<<' '<<i*50<<' '<<(i<7?1:2)<<' '<<(i<7?i:i-7)*50<<' '<<((i<7?i:i-7)%2?"NA":"1")<<" 2 "<<((i<7?i:i-7)%2?"absent":"exact")<<'\n';
 f<<"exchange 0 1 150\nexchange 0 1 175\n";f.close();
 auto msm=std::make_unique<TitanMSM>(2,1,.05);TitanExecutor executor(1);
 for(int i=0;i<10;++i){Eigen::VectorXd x(1);x[0]=states[i]*10.;msm->add_projected_frame(x,0,i,2.);msm->record_frame_time(i*50.,0);}
 msm->train_kmeans(executor,50,42);msm->build_markov_model(2);
 arguments args;args.dham.flag_dham=true;args.dham.segment_map_file="segments.txt";args.dham.temps_file="temps.txt";args.dham.ensemble="npt";
 info output(args);std::fstream log("pipeline.log",std::ios::out);std::vector<double> weights{1};Eigen::MatrixXd tpm(1,1);Eigen::VectorXd pi(1);bool converged=true;TitanValidation::DhamRebuildContext context;context.available=true;
 run_dham_postproc(args,msm,executor,output,log,weights,&tpm,&pi,&converged,&context);log.flush();
 // KIN-13: segment mode feeds TRAM: samples are the frames with an exact
 // energy (6), transitions the eligible windows (3). These few frames do
 // not determine a model, so TRAM reports non-convergence and no weights,
 // and no bootstrap context is offered in segment mode.
 require(weights.empty()&&tpm.size()==0&&pi.size()==0&&!converged&&!context.available,"non-converged TRAM leaked outputs");
 {std::ifstream summary("dham_tram_summary.dat");std::string s((std::istreambuf_iterator<char>(summary)),{});
  require(s.find("input_mode segment_map")!=std::string::npos&&s.find("\nsamples 6\n")!=std::string::npos&&s.find("\ntransitions 3\n")!=std::string::npos,"TRAM did not receive the exact-energy samples and eligible windows");}
 std::ifstream decisions("dham_segment_decisions.csv");std::string text((std::istreambuf_iterator<char>(decisions)),{});decisions.close();
 require(text.find("0,eligible,3")!=std::string::npos&&text.find("0,exchange,3")!=std::string::npos&&text.find("0,segment_boundary,2")!=std::string::npos,"pipeline did not apply metadata masks");
 require(!std::filesystem::exists("dham_free_energy_profile.dat"),"non-converged TRAM wrote a FES");
 // An incompatible trajectory timestamp must reject the metadata and emit
 // no new count artifact, even though a nearest match would be within 10 ps.
 std::filesystem::rename("dham_segment_decisions.csv","accepted_decisions.csv");std::filesystem::rename("dham_segment_counts.csv","accepted_counts.csv");
 auto wrong=std::make_unique<TitanMSM>(2,1,.05);
 for(int i=0;i<10;++i){Eigen::VectorXd x(1);x[0]=states[i]*10.;wrong->add_projected_frame(x,0,i,2.);wrong->record_frame_time(i*50.+10.,0);}
 wrong->train_kmeans(executor,50,42);wrong->build_markov_model(2);
 run_dham_postproc(args,wrong,executor,output,log,weights,&tpm,&pi,&converged,&context);log.flush();
 require(!std::filesystem::exists("dham_segment_counts.csv")&&!converged&&!context.available,"bad timestamp accepted by pipeline");
 std::cout<<"PASS production DHAM segment pipeline: TRAM samples and eligible windows, no outputs without convergence, bad timestamp rejected\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
