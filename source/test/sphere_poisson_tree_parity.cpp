#include "tabi/tabi_types.h"
#include "titan_executor.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace TitanPBE::Tabi;
int main() {
  try {
    TitanExecutor executor(4);Options o;
    o.salt_m=0;o.eps_out=78.5;o.surface_point_count=1280;
    o.sphere_mesh_mode="icosphere";o.near_integration="auto";
    o.residual_tol=1e-10;o.max_iter=300;o.preconditioner="diagonal";
    o.sum_mode="direct";
    std::vector<Atom> atoms={{0,0,0,0,1.7},{1.36,0,0,1,.017}};
    auto direct=solve_linear_pb(atoms,o,executor);
    if(!direct.converged)throw std::runtime_error("direct failed");
    o.formulation="jcp2013-treecode";o.sum_mode="treecode";
    o.tree_builder="topdown-aabb";o.tree_leaf_size=16;o.tree_leaf_size_explicit=true;
    o.tree_opening_theta_explicit=true;
    std::cout<<"builder,theta,far_interactions,energy,relative_to_direct,near_policy\n"<<std::setprecision(17);
    for(const char *builder:{"grouped-morton","topdown-aabb"}) {
      o.tree_builder=builder;
      for(double theta:{.05,.25}) {
        o.tree_opening_theta=theta;
        auto r=solve_linear_pb(atoms,o,executor);
        double error=std::abs(r.energy_kcal_mol/direct.energy_kcal_mol-1);
        std::cout<<builder<<','<<theta<<','<<r.far_interactions<<','<<r.energy_kcal_mol<<','<<error<<','<<r.near_kernel_policy<<'\n';
        if(!r.converged || error>1e-4 || r.near_duffy_pairs==0 ||
           r.near_kernel_policy!="sphere_poisson_curved")
          throw std::runtime_error("tree parity failed");
        if(std::string(builder)=="topdown-aabb" && theta==.25 && r.far_interactions==0)
          throw std::runtime_error("far path was not exercised");
      }
    }
    return 0;
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
