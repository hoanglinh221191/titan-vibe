#include "tabi/tabi_types.h"
#include "titan_executor.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

// Atom input: count, then x y z radius charge. Fixed solver settings isolate
// mesh sensitivity; convergence of GMRES is distinct from discretization error.
int main(int argc, char **argv) {
  try {
    if (argc != 4) return 2;
    std::ifstream input(argv[1]);
    std::size_t count{};
    if (!(input >> count) || count == 0) return 2;
    std::vector<TitanPBE::Tabi::Atom> atoms(count);
    for (auto &a : atoms)
      if (!(input >> a.x >> a.y >> a.z >> a.radius >> a.q)) return 2;
    TitanPBE::Tabi::Options o;
    o.primary_surface = "sphere-patch";
    o.surface_mode = "vdw";
    o.sphere_mesh_mode = argv[2];
    o.surface_point_count = std::stoi(argv[3]);
    o.eps_in = 1;
    o.eps_out = 78.5;
    o.salt_m = 0;
    o.sum_mode = "direct";
    o.residual_tol = 1e-7;
    o.max_iter = 300;
    o.preconditioner = "diagonal";
    o.near_integration = "analytic_singular";
    o.discretization = "centroid";
    TitanExecutor executor(4);
    const auto begin = std::chrono::steady_clock::now();
    auto r = TitanPBE::Tabi::solve_linear_pb(atoms, o, executor);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    std::cout << std::setprecision(17) << "RESULT " << r.converged << ' '
              << r.energy_kcal_mol << ' ' << r.iterations << ' '
              << r.true_residual_rel << ' ' << seconds << '\n';
    return r.converged ? 0 : 1;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
