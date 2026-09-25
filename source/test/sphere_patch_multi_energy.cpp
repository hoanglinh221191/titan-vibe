#include "tabi/tabi_types.h"
#include "titan_executor.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace TitanPBE::Tabi;
int main(int argc, char **argv) {
  try {
    std::cout << std::unitbuf;
    if (argc < 4 || argc > 11)
      return 2;
    std::ifstream in(argv[1]);
    std::size_t n{};
    if (!(in >> n) || !n)
      return 2;
    std::vector<Atom> atoms(n);
    for (auto &a : atoms)
      if (!(in >> a.x >> a.y >> a.z >> a.radius >> a.q))
        return 2;
    Options o;
    o.discretization = argc > 10 ? argv[10] : "centroid";
    o.surface_mode = "vdw";
    o.sphere_mesh_mode = argc > 9 ? argv[9] : "icosphere";
    o.surface_point_count = std::stoi(argv[2]);
    o.highres_surface = true;
    o.eps_in = 1;
    o.eps_out = 78.5;
    o.salt_m = argc > 6 ? std::stod(argv[6]) : 0;
    o.near_integration = "auto";
    o.residual_tol = argc > 4 ? std::stod(argv[4]) : 1e-8;
    o.max_iter = 300;
    o.preconditioner = "diagonal";
    o.sum_mode = argv[3];
    if (o.sum_mode == "treecode") {
      o.formulation = "jcp2013-treecode";
      o.tree_builder = "topdown-aabb";
      o.tree_leaf_size = 16;
      o.tree_leaf_size_explicit = true;
      o.tree_opening_theta = argc > 5 ? std::stod(argv[5]) : .25;
      o.tree_opening_theta_explicit = true;
    }
    const int threads = argc > 7 ? std::stoi(argv[7]) : 4;
    if (threads < 1 || threads > 256)
      return 2;
    if (argc > 8)
      o.progress = argv[8];
    TitanExecutor executor(static_cast<std::size_t>(threads));
    const auto begin = std::chrono::steady_clock::now();
    auto r = solve_linear_pb(atoms, o, executor);
    std::cout << std::setprecision(17)
              << "RESULT {\"converged\":" << (r.converged ? "true" : "false")
              << ",\"energy\":" << r.energy_kcal_mol
              << ",\"iterations\":" << r.iterations
              << ",\"threads\":" << threads << ",\"mesh\":\""
              << o.sphere_mesh_mode << "\""
              << ",\"discretization\":\"" << r.discretization << "\""
              << ",\"residual\":" << r.true_residual_rel
              << ",\"panels\":" << r.panel_count
              << ",\"area\":" << r.surface_area_a2
              << ",\"mesh_ms\":" << r.surface_mesh_ms
              << ",\"min_d_over_h\":" << r.charge_resolution.min_ratio
              << ",\"inside_atom_count\":"
              << r.charge_resolution.surface_inside_atom_count
              << ",\"far_interactions\":" << r.far_interactions
              << ",\"seconds\":"
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - begin)
                     .count()
              << "}\n";
    return r.converged ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
