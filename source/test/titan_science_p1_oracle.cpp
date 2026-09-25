#include "box_type.h"
#include "titan_neighbor_grid.h"
#include "titan_qha_model.h"
#include "titan_kinetics.h"
#include "titan_tica_backend_policy.h"
#include "titan_markov_estimator.h"
#include "titan_dham.h"

#include <array>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

void tica_oracle() {
  TitanKinetics alternating(1,1,1,0.001,TicaBackend::ExactDense,1);
  for(double value:{1.0,-1.0,1.0,-1.0,1.0,-1.0}) {
    Eigen::VectorXd x(1);x[0]=value;alternating.accumulate_frame(x);
  }
  alternating.solve_tica();
  require(alternating.is_solved() &&
          std::abs(alternating.eigenvalues()[0]+1.0/(1+1e-6))<1e-10,
          "P0 TICA lag buffer must retain the previous frame before overwrite");
  const std::vector<std::vector<double>> segments{{0,2,1,4,3,-1}, {10,-3,2}, {1e6}};
  for (int lag : {1,2}) {
    double sum=0,squares=0,product=0,pairs=0;
    for (const auto &segment:segments) for(std::size_t i=lag;i<segment.size();++i) {
      double a=segment[i-lag],b=segment[i];
      sum+=a+b;squares+=a*a+b*b;product+=a*b;++pairs;
    }
    const double mean=sum/(2*pairs);
    const double c0=squares/(2*pairs)-mean*mean;
    const double ct=product/pairs-mean*mean;
    const double expected=ct/(c0*(1+1e-6));
    for (auto backend:{TicaBackend::ExactDense,TicaBackend::ExactDual,TicaBackend::Randomized,TicaBackend::ExactStreaming}) {
      TitanKinetics tica(1,lag,1,0.001,backend,1);
      for(const auto &segment:segments) {
        tica.reset_lag_buffer();
        for(double value:segment) {Eigen::VectorXd x(1);x[0]=value;tica.accumulate_frame(x);}
      }
      tica.solve_tica();
      require(tica.is_solved(),"TICA did not solve");
      require(std::abs(tica.mean()[0]-mean)<1e-12,"TICA must center paired endpoints only");
      require(std::abs(tica.eigenvalues()[0]-expected)<1e-10,"TICA lag covariance differs from explicit pairs");
      if(backend==TicaBackend::ExactStreaming) {
        int replayed=0;
        tica.project_all_frames_to_sink([&](int64_t,const Eigen::Ref<const Eigen::VectorXd>&x){require(x.allFinite(),"streaming projection finite");++replayed;});
        require(replayed==10 && tica.raw_feature_cache_rows()==0,"exact streaming disk replay without raw RAM cache");
      }
    }
  }
  TitanTicaBackendPolicyConfig cfg;cfg.available_ram_bytes=24.0L*1024*1024*1024;
  cfg.estimated_frames=300000;cfg.feature_dim=1500;cfg.need_projection_cache=true;
  cfg.dense_memory_estimate=6.0L*1500*1500*8;cfg.dual_memory_estimate=3.0L*300000*1500*8;
  require(select_tica_backend_for_policy(cfg).selected_backend==TicaBackend::ExactStreaming,
          "large auto job must not silently truncate PCA");
  cfg.requested_backend=TicaBackend::ExactDual;
  bool rejected=false;try{select_tica_backend_for_policy(cfg);}catch(const std::exception&){rejected=true;}
  require(rejected,"explicit exact backend must not silently switch to approximate");
  // Full matrix oracle: use explicit endpoint rows and solve the independent
  // generalized eigensystem. This detects transpose/aliasing errors hidden in 1D.
  Eigen::MatrixXd samples(8,2);samples<<0,1, 2,-1, 1,4, -3,2, 4,1, 2,5, -1,-2, 3,4;
  Eigen::MatrixXd left=samples.topRows(7),right=samples.bottomRows(7);
  Eigen::RowVectorXd mu=(left.colwise().sum()+right.colwise().sum())/14;
  left.rowwise()-=mu;right.rowwise()-=mu;
  Eigen::MatrixXd c0=(left.transpose()*left+right.transpose()*right)/14;
  Eigen::MatrixXd ct=(left.transpose()*right+right.transpose()*left)/14;
  c0.diagonal().array()+=1e-6*c0.trace()/2;
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> oracle(ct,c0);
  for(auto backend:{TicaBackend::ExactDense,TicaBackend::ExactDual,TicaBackend::ExactStreaming}) {
    TitanKinetics model(2,1,2,0.001,backend,2);
    for(int row=0;row<samples.rows();++row)model.accumulate_frame(samples.row(row).transpose());
    model.solve_tica();
    require((model.eigenvalues()-oracle.eigenvalues().reverse()).norm()<1e-9,
            "multivariate TICA generalized eigenvalue oracle");
  }
  std::cout<<"PASS TICA circular-buffer lag, reversible paired estimator, all 4 backends, exact memory policy\n";
}

void markov_oracle() {
  Eigen::MatrixXd c(2,2);c<<9,1,9,1;
  Eigen::MatrixXd t;Eigen::VectorXd p;std::string error;
  require(titan_markov::reversible_mle(c,t,p,error),"reversible MLE convergence");
  require(std::abs(t(0,1)-0.1)<1e-9 && std::abs(t(1,0)-0.9)<1e-9 &&
          std::abs(p[0]-0.9)<1e-9,"2-state reversible MLE must equal unconstrained closed form");
  c<<1-1e-8,1e-8,4e-8,1-4e-8;
  require(titan_markov::stationary(c,p,error) && std::abs(p[0]-0.8)<1e-9,
          "metastable stationary distribution without power truncation");
  c=Eigen::MatrixXd::Identity(2,2);
  require(!titan_markov::reversible_mle(c,t,p,error),"reject disconnected MSM counts");
  c.resize(3,3);c<<80,12,2, 4,40,16, 5,3,30;
  require(titan_markov::reversible_mle(c,t,p,error),"3-state reversible MLE convergence");
  Eigen::MatrixXd flux=p.asDiagonal()*t;
  require((flux-flux.transpose()).norm()<1e-12,"MLE detailed balance");
  auto likelihood=[&](Eigen::MatrixXd x){for(int i=0;i<3;++i)x.row(i)/=x.row(i).sum();return (c.array()*x.array().log()).sum();};
  const double optimum=likelihood(flux);
  for(int i=0;i<3;++i)for(int j=i;j<3;++j)for(double step:{-1e-5,1e-5}) {
    Eigen::MatrixXd trial=flux;trial(i,j)+=step;if(i!=j)trial(j,i)+=step;
    require(likelihood(trial)<=optimum+1e-10,"MLE likelihood must not increase under symmetric flux perturbations");
  }
  std::cout<<"PASS reversible MLE, metastable stationary solve, disconnected rejection\n";
}

void dham_oracle() {
  auto build=[](double temp,double offset) {
    TitanDHAM dham(2,1,{temp},300);
    const int counts[2][2]={{90,10},{10,40}};
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int n=0;n<counts[i][j];++n)dham.accumulate_counts(0,i,j);
    dham.set_state_energy(0,0,offset);dham.set_state_energy(0,1,5+offset);
    require(dham.solve(),"DHAM solve");
    return dham.get_unbiased_transition_matrix();
  };
  auto t=build(300,0);
  require(std::abs(t(0,1)-0.1)<1e-12 && std::abs(t(1,0)-0.2)<1e-12,"DHAM target-T identity");
  auto biased=build(400,0),shifted=build(400,1000);
  require((biased-shifted).cwiseAbs().maxCoeff()<1e-12,"DHAM energy offset invariance");
  const double delta=(1/(DHAM_kB_kJ*400)-1/(DHAM_kB_kJ*300))*5;
  const double a=0.1*std::exp(0.5*delta),b=0.2*std::exp(-0.5*delta);
  require(std::abs(biased(0,1)-a/(0.9+a))<1e-12 &&
          std::abs(biased(1,0)-b/(0.8+b))<1e-12,"DHAM pair bias denominator");
  std::cout<<"PASS DHAM unbiased identity, explicit pair-bias formula, offset invariance\n";
}

void geometry_oracle() {
  for(const auto angles:{std::array<double,3>{90,90,120},{90,105,90},{75,80,65},{90+1e-12,90,90}}) {
    box b(10.0,12.0,14.0,double(angles[0]),double(angles[1]),double(angles[2]));
    require(b.get_box_type()!="unknown","valid lattice classified unknown");
    double a[3],c[3],d[3];b.get_v1(a);b.get_v2(c);b.get_v3(d);
    auto dot=[](const double*x,const double*y){return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];};
    require(std::abs(dot(a,a)-100)<1e-10 && std::abs(dot(c,c)-144)<1e-10 &&
            std::abs(dot(d,d)-196)<1e-10,"lattice vector lengths");
    double gram=1;
    std::array<double,3> cosines{};
    for(int i=0;i<3;++i){cosines[i]=std::cos(angles[i]*std::acos(-1.0)/180);gram-=cosines[i]*cosines[i];}
    gram+=2*cosines[0]*cosines[1]*cosines[2];
    require(std::abs(b.volume()-10*12*14*std::sqrt(gram))<1e-8,"cell volume determinant");
  }
  std::mt19937 rng(19);std::uniform_real_distribution<double> uniform(-4,4);
  for(double lx:{1.5,2.5,3.5}) for(double ly:{1.5,2.5,3.5}) for(double lz:{1.5,2.5,3.5}) {
    std::vector<std::array<double,3>> points(50);
    for(auto &p:points)p={uniform(rng)*lx,uniform(rng)*ly,uniform(rng)*lz};
    titan_neighbor::OrthorhombicGridSpec spec;spec.cutoff=1;spec.len_x=lx;spec.len_y=ly;spec.len_z=lz;
    auto near=[&](std::size_t i,std::size_t j){double sum=0;for(int k=0;k<3;++k){double l=std::array<double,3>{lx,ly,lz}[k];double d=points[i][k]-points[j][k];d-=std::round(d/l)*l;sum+=d*d;}return sum<=1;};
    std::map<std::pair<std::size_t,std::size_t>,int> actual,expected;
    for(std::size_t i=0;i<points.size();++i)for(std::size_t j=i+1;j<points.size();++j)if(near(i,j))expected[{i,j}]=1;
    titan_neighbor::HalfShellGrid grid;
    grid.for_each_pair(points.size(),spec,[&](std::size_t i,double&x,double&y,double&z){x=points[i][0];y=points[i][1];z=points[i][2];},
       [](std::size_t){return true;},[&](std::size_t i,std::size_t j){require(i!=j,"neighbor grid emitted self pair");if(near(i,j))++actual[std::minmax(i,j)];});
    require(actual==expected,"periodic short-cell pair multiplicity/coverage");
  }
  std::cout<<"PASS monoclinic/hexagonal cell vectors/volume, 27 short-cell neighbor oracles\n";
}

void qha_oracle() {
  constexpr double hbar=1.054571817e-34,kb=1.380649e-23,amu_a2=1.66053906660e-47;
  for(double temp:{200.0,300.0,500.0}) {
    // Choose a mode with hbar*omega/(k_B*T)=1. Sum its occupation entropy
    // directly from the normalized geometric quantum-state probabilities.
    const double lambda=hbar*hbar/(kb*temp*amu_a2);
    double entropy=0;
    for(int n=0;n<100;++n){double p=(1-std::exp(-1.0))*std::exp(-double(n));entropy-=p*std::log(p);}
    const double reference=0.00831446261815324*temp*entropy;
    require(std::abs(titan_qha_ts_kj_mol({lambda},temp)-reference)<1e-12,"QHA oscillator entropy/units");
  }
  require(std::isfinite(titan_qha_ts_kj_mol({1e-8,1e8},300)),"QHA stable limits");
  std::cout<<"PASS QHA quantum occupation entropy and extreme-mode limits\n";
}

int main(int argc,char**argv) {
  try {
    const std::string which=argc>1?argv[1]:"all";
    if(which=="all"||which=="tica")tica_oracle();
    if(which=="all"||which=="geometry")geometry_oracle();
    if(which=="all"||which=="qha")qha_oracle();
    if(which=="all"||which=="markov")markov_oracle();
    if(which=="all"||which=="dham")dham_oracle();
    return 0;
  }catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
