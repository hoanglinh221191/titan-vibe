#include "pbe/grid/grid_workspace.h"

#include <cassert>
#include <cmath>
#include <cstddef>

int main() {
  TitanPBE::GridOptions options;
  options.ext_diel = 78.5;

  TitanPBE::GridWorkspace workspace;
  assert(!workspace.sized(3));
  assert(workspace.memory_bytes() == 0);

  workspace.reset(3, options);
  assert(workspace.sized(3));
  assert(workspace.phi.size() == 3);
  assert(workspace.rho.size() == 3);
  assert(workspace.eps_x.size() == 3);
  assert(workspace.eps_y.size() == 3);
  assert(workspace.eps_z.size() == 3);
  assert(workspace.denom.size() == 3);
  assert(workspace.screen.size() == 3);
  assert(workspace.inside.size() == 3);
  assert(workspace.ion_accessible.size() == 3);

  for (std::size_t i = 0; i < 3; ++i) {
    assert(workspace.phi[i] == 0.0);
    assert(workspace.rho[i] == 0.0);
    assert(workspace.screen[i] == 0.0);
    assert(workspace.denom[i] == 1.0);
    assert(std::abs(workspace.eps_x[i] - options.ext_diel) < 1.0e-12);
    assert(std::abs(workspace.eps_y[i] - options.ext_diel) < 1.0e-12);
    assert(std::abs(workspace.eps_z[i] - options.ext_diel) < 1.0e-12);
    assert(workspace.inside[i] == 0U);
    assert(workspace.ion_accessible[i] == 1U);
  }

  const std::size_t bytes3 = workspace.memory_bytes();
  assert(bytes3 > 0);

  workspace.reset(5, options);
  assert(workspace.sized(5));
  assert(workspace.memory_bytes() > bytes3);

  workspace.clear();
  assert(!workspace.sized(1));
  assert(workspace.phi.empty());
  assert(workspace.rho.empty());
  assert(workspace.eps_x.empty());
  assert(workspace.eps_y.empty());
  assert(workspace.eps_z.empty());
  assert(workspace.denom.empty());
  assert(workspace.screen.empty());
  assert(workspace.inside.empty());
  assert(workspace.ion_accessible.empty());
  assert(workspace.memory_bytes() == 0);

  return 0;
}
