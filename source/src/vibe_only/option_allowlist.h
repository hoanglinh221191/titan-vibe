#pragma once
// TITAN-VIBE (VIBE_ONLY=1): the command-line options this build accepts.
//
// TITAN-VIBE keeps the parts of TITAN that a VIBE run uses: the structure,
// index, topology and trajectory readers, the force-field library that types
// atoms when no topology file is given, VIBE itself (-pbctopo*) and the PDB
// writer. Every other option belongs to a module that is not compiled in, and
// parse_arguments() refuses it before any option family sees it.

#include <string_view>

namespace titan_vibe_only {

// The option name of a command-line token: "-opdb=placed.pdb" -> "-opdb".
std::string_view option_name(std::string_view token) noexcept;

// True if this build accepts the token. Tokens that do not start with '-' are
// left to the parser, which treats them as before.
bool option_allowed(std::string_view token) noexcept;

// Prints why the token is refused.
void report_refused_option(std::string_view token);

} // namespace titan_vibe_only
