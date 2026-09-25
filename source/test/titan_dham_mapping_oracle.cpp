#include "titan_dham.h"
#include <iostream>
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
 const int counts[3][2][2]={{{90,10},{8,42}},{{60,20},{15,55}},{{40,30},{20,60}}};
 const double temps[3]={300,310,320},energy[3][2]={{-100,-99},{-101,-99.5},{-103,-100}};
 auto solve=[&](double offset){TitanDHAM d(2,3,{300,310,320},300);for(int k=0;k<3;++k)for(int i=0;i<2;++i){d.set_state_energy(k,i,energy[k][i]+offset);for(int j=0;j<2;++j)for(int n=0;n<counts[k][i][j];++n)d.accumulate_counts(k,i,j);}require(d.solve(),"multi-replica solve");return d.get_unbiased_transition_matrix();};
 const auto actual=solve(0),shifted=solve(10000);
 Eigen::MatrixXd expected(2,2);
 for(int i=0;i<2;++i){for(int j=0;j<2;++j){double numerator=0,denominator=0;for(int k=0;k<3;++k){numerator+=counts[k][i][j];double beta=1/(DHAM_kB_kJ*temps[k])-1/(DHAM_kB_kJ*300);denominator+=(counts[k][i][0]+counts[k][i][1])*std::exp(-.5*beta*(energy[k][j]-energy[k][i]));}expected(i,j)=numerator/denominator;}expected.row(i)/=expected.row(i).sum();}
 require((actual-expected).norm()<1e-12 && (actual-shifted).norm()<1e-12,"independent multi-replica denominator/offset oracle");
 std::unordered_map<int,std::vector<int>> traj{{0,{0,1,0}},{1,{1,0,1}}};
 std::unordered_map<int,std::vector<double>> values{{0,{2,4,6}},{1,{10,20,30}}};
 auto means=TitanDHAM::compute_state_energies(traj,values,2,2,2);
 require(means(0,0)==4&&means(1,0)==4&&means(0,1)==20&&means(1,1)==20,"state/replica mapping");
 int failures=0;
 try{auto bad=values;bad[0].pop_back();TitanDHAM::compute_state_energies(traj,bad,2,2,2);}catch(const std::invalid_argument&){++failures;}
 try{auto bad=traj;bad[0][0]=2;TitanDHAM::compute_state_energies(bad,values,2,2,2);}catch(const std::out_of_range&){++failures;}
 TitanDHAM bad(2,1,{400},300);
 try{bad.accumulate_counts(1,0,1);}catch(const std::out_of_range&){++failures;}
 try{bad.set_state_energy(0,0,std::numeric_limits<double>::infinity());}catch(const std::invalid_argument&){++failures;}
 for(int i=0;i<2;++i)for(int j=0;j<2;++j)bad.accumulate_counts(0,i,j);
 require(!bad.solve(),"missing energy accepted");require(failures==4,"invalid mapping accepted");
 std::cout<<"PASS DHAM three replicas, independent bias denominator, offset invariance, threaded state mapping, invalid indices/lengths/nonfinite/missing energy rejection\n";return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
