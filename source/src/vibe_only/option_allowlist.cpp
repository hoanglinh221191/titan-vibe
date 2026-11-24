#include "vibe_only/option_allowlist.h"

#include "titan_ui.h"

#include <algorithm>
#include <iterator>
#include <string>

namespace titan_vibe_only {

namespace {

// Kept sorted; option_allowed() uses binary search.
constexpr std::string_view kAllowed[] = {
    "--help",                      // the VIBE reference (-vibe-help)
    "--version",
    "--vibe-help",
    "-1pdb",                       // one PDB file for all frames
    "-a",                          // argument file
    "-b",                          // first frame time
    "-common_topo",                // GROMACS topology for both molecules
    "-dt",                         // frame stride
    "-e",                          // last frame time
    "-eu",                         // energy unit (no energies here; kept for argument files)
    "-f",                          // structure or frame list
    "-ff",                         // force-field library
    "-fflib-auto",                 // atom typing from the library when no topology is given
    "-fflib-hydrogens",
    "-fflib-missing",
    "-fu",                         // force unit (as -eu)
    "-g",                          // index groups of the two molecules
    "-h",
    "-n",                          // index file
    "-no-dssp",
    "-nocal",                      // coordinate run; TITAN-VIBE requires it
    "-noprint-side_chain_contact",
    "-nt",                         // threads
    "-opdb",                       // output PDB
    "-opdbsep",                    // one PDB file per frame
    "-output-manifest",
    "-path-a",                     // force-field library paths
    "-path-n",
    "-s",                          // reference structure
    "-topo",                       // GROMACS topology of molecule 1
    "-topo2",                      // GROMACS topology of molecule 2
    "-traj",                       // trajectory (xtc, trr)
    "-tu",                         // time unit of -b, -e, -dt
    "-v",
    "-vibe-help",
};
static_assert(std::is_sorted(std::begin(kAllowed), std::end(kAllowed)));

// Every -pbctopo option (VIBE itself), in either spelling.
constexpr std::string_view kAllowedPrefixes[] = {"-pbctopo", "--pbctopo"};

} // namespace

std::string_view option_name(std::string_view token) noexcept {
  const std::size_t eq = token.find('=');
  return eq == std::string_view::npos ? token : token.substr(0, eq);
}

bool option_allowed(std::string_view token) noexcept {
  if (token.empty() || token.front() != '-')
    return true;
  const std::string_view name = option_name(token);
  if (std::binary_search(std::begin(kAllowed), std::end(kAllowed), name))
    return true;
  return std::any_of(std::begin(kAllowedPrefixes), std::end(kAllowedPrefixes),
                     [name](std::string_view prefix) { return name.starts_with(prefix); });
}

void report_refused_option(std::string_view token) {
  TitanUI::print_error("TITAN-VIBE: option " + std::string(option_name(token)) +
                       " belongs to a TITAN module that is not part of this build. TITAN-VIBE "
                       "reads structures, index files, topologies and trajectories, runs VIBE "
                       "(-pbctopo) and writes PDB files; see -vibe-help.");
}

} // namespace titan_vibe_only
