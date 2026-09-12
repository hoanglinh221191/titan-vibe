#include "pbe/grid/grid_solver.h"

#include <cassert>
#include <string_view>

int main() {
  using TitanPBE::GridSolveMode;
  using TitanPBE::GridSolveRequest;

  const GridSolveRequest linear =
      GridSolveRequest::linear(2.0, 80.0, nullptr);
  assert(linear.mode == GridSolveMode::Linear);
  assert(linear.eps_in == 2.0);
  assert(linear.eps_out == 80.0);
  assert(linear.kappa == 0.0);
  assert(linear.map_result == nullptr);

  const GridSolveRequest physical =
      GridSolveRequest::physical(4.0, 78.5, 0.125, nullptr);
  assert(physical.mode == GridSolveMode::Physical);
  assert(physical.eps_in == 4.0);
  assert(physical.eps_out == 78.5);
  assert(physical.kappa == 0.125);

  const GridSolveRequest nonlinear =
      GridSolveRequest::nonlinear(1.0, 60.0, 0.25, nullptr);
  assert(nonlinear.mode == GridSolveMode::Nonlinear);
  assert(nonlinear.kappa == 0.25);

  assert(TitanPBE::grid_solve_mode_name(GridSolveMode::Linear) ==
         std::string_view("linear"));
  assert(TitanPBE::grid_solve_mode_name(GridSolveMode::Physical) ==
         std::string_view("physical"));
  assert(TitanPBE::grid_solve_mode_name(GridSolveMode::Nonlinear) ==
         std::string_view("nonlinear"));

  return 0;
}
