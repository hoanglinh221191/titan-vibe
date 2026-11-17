// PATH-01 oracle: -de-pathway-radius is shorthand for the pathway graph edge
// cutoff and the route-recovery band. It used to set the recovery enter and
// exit radii both to r (a band of zero width: Advance/Recovery could
// alternate every step at the boundary); the exit radius must stay below the
// enter radius, and explicit -pathway-recovery-* / -pathway-graph-edge-cutoff
// options must win in any order. Linked against the engine objects
// (test/engine_linked_oracles.py, case path01).
#include "arguments.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {

struct Case {
  const char *name;
  std::vector<std::string> options;
  double edge, enter, exit;
};

bool close(double a, double b) { return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b)); }

} // namespace

int main() {
  {
    std::ofstream pdb("path01.pdb");
    pdb << "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
           "ATOM      2  CA  ALA B   2       5.000   0.000   0.000  1.00  0.00           C\n"
           "END\n";
  }
  const std::vector<Case> cases = {
      {"defaults", {}, 4.0, 4.0, 3.0},
      {"alias", {"-de-pathway-radius", "6"}, 6.0, 6.0, 4.5},
      {"alias_then_exit", {"-de-pathway-radius", "6", "-pathway-recovery-exit", "2"}, 6.0, 6.0, 2.0},
      {"exit_then_alias", {"-pathway-recovery-exit", "2", "-de-pathway-radius", "6"}, 6.0, 6.0, 2.0},
      {"enter_then_alias", {"-pathway-recovery-enter", "8", "-de-pathway-radius", "6"}, 6.0, 8.0, 6.0},
      {"edge_then_alias", {"-pathway-graph-edge-cutoff", "5", "-de-pathway-radius", "6"}, 5.0, 6.0, 4.5},
  };
  int failures = 0;
  for (const Case &c : cases) {
    std::vector<std::string> tokens = {"titan", "-f", "path01.pdb", "-1pdb"};
    tokens.insert(tokens.end(), c.options.begin(), c.options.end());
    std::vector<char *> argv;
    for (auto &t : tokens)
      argv.push_back(t.data());
    arguments args{};
    // As in a run: parse_arguments, then post_process (engine/run_setup.cpp),
    // where the alias is resolved after every explicit option has been read.
    int status = args.parse_arguments(static_cast<int>(argv.size()), argv.data());
    if (status == 0)
      status = args.post_process();
    const auto &de = args.de;
    const bool ok = status == 0 && close(de.path_graph_edge_cutoff, c.edge) &&
                    close(de.path_recovery_enter_radius, c.enter) &&
                    close(de.path_recovery_exit_radius, c.exit) &&
                    de.path_recovery_exit_radius < de.path_recovery_enter_radius;
    std::printf("%s %-18s status %d edge %.3f enter %.3f exit %.3f (want %.3f %.3f %.3f)\n",
                ok ? "PASS" : "FAIL", c.name, status, de.path_graph_edge_cutoff,
                de.path_recovery_enter_radius, de.path_recovery_exit_radius, c.edge, c.enter,
                c.exit);
    failures += ok ? 0 : 1;
  }
  std::printf("%s PATH-01 pathway radius alias\n", failures ? "FAIL" : "PASS");
  return failures ? 1 : 0;
}
