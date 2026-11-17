// KIN-23 (a, b) oracle: the input readers of -mbar-demux-map, -tpt-demux-map
// and -dham-temps / -mbar-temps on small files.
//
//   temps   TitanDHAM::read_temps
//     header   "# ...", "@ ...", blank, then 300 / 310.5 / "320 330":
//              all four temperatures (old: stopped at "#", none)
//     commas   "300, 310,320": three (old: one)
//     text     300, 310, "T=320": refused naming file:3 (old: 300, 310)
//     negative 300, -5: refused naming file:2 (old: 300, an error, went on)
//   tpt     postproc_detail::read_tpt_demux_map
//     ints     "@ xvg", "# c", "0 1 2", "1 2 0": two rows of three ids
//     time     GROMACS replica_index.xvg row "0.000 0 1 2": refused naming
//              file:1 (old: the row collapsed to the single id 0)
//     token    "0 1 x": refused naming file:2 (old: row [0 1], column error)
//   mbar    titan_fe::DemuxStateMap::read
//     token    "0 1 x" in a 3-state map: refused naming file:2
// Every refusal must name the file and the line. Linked against the engine
// objects (test/engine_linked_oracles.py).
#include "free_energy/mbar_demux_map.h"
#include "postproc/postproc_common.h"
#include "titan_dham.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures = 0;
void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok   " : "  FAIL ") << what << '\n';
  if (!ok)
    ++failures;
}

void write(const std::string &path, const std::string &text) {
  std::ofstream(path, std::ios::binary) << text;
}

// read_temps with std::cerr captured
std::vector<double> temps_of(const std::string &path, std::string &stderr_text) {
  std::ostringstream capture;
  auto *old = std::cerr.rdbuf(capture.rdbuf());
  auto temps = TitanDHAM::read_temps(path);
  std::cerr.rdbuf(old);
  stderr_text = capture.str();
  return temps;
}

bool same(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.size() != b.size())
    return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::abs(a[i] - b[i]) > 1e-12)
      return false;
  return true;
}

std::string show(const std::vector<double> &v) {
  std::ostringstream o;
  o << '[';
  for (std::size_t i = 0; i < v.size(); ++i)
    o << (i ? " " : "") << v[i];
  o << ']';
  return o.str();
}

bool names(const std::string &message, const std::string &file, int line) {
  return message.find(file + ":" + std::to_string(line)) != std::string::npos;
}

// the last error line of captured stderr, without the colour codes
std::string message(const std::string &captured) {
  const auto at = captured.rfind("[ERROR]");
  if (at == std::string::npos)
    return "no error message";
  std::string text = captured.substr(at + 7);
  const auto reset = text.find("\033[0m");
  if (reset != std::string::npos)
    text.erase(reset, 4);
  while (!text.empty() && (text.back() == '\n' || text.back() == ' '))
    text.pop_back();
  while (!text.empty() && text.front() == ' ')
    text.erase(text.begin());
  return text;
}

} // namespace

int main() {
  std::string err;
  std::cout << "case temps\n";
  write("temps_header.dat", "# replica ladder\n@ title \"T\"\n\n300\n310.5\n320 330\n");
  auto t = temps_of("temps_header.dat", err);
  check(same(t, {300, 310.5, 320, 330}),
        "header, xvg directive and blank lines skipped: " + show(t));
  write("temps_commas.dat", "300, 310,320\n");
  t = temps_of("temps_commas.dat", err);
  check(same(t, {300, 310, 320}), "comma-separated temperatures: " + show(t));
  write("temps_text.dat", "300\n310\nT=320\n");
  t = temps_of("temps_text.dat", err);
  check(t.empty() && names(err, "temps_text.dat", 3),
        "non-numeric line refused with file:line: " + show(t) + " / " + message(err));
  write("temps_negative.dat", "300\n-5\n");
  t = temps_of("temps_negative.dat", err);
  check(t.empty() && names(err, "temps_negative.dat", 2),
        "temperature <= 0 refused with file:line: " + show(t) + " / " + message(err));

  std::cout << "case tpt demux map\n";
  std::vector<std::vector<int>> rows;
  std::string error;
  write("tpt_ints.dat", "@ title \"replica index\"\n# comment\n0 1 2\n1 2 0\n");
  bool ok = postproc_detail::read_tpt_demux_map("tpt_ints.dat", rows, error);
  check(ok && rows.size() == 2 && rows[0] == std::vector<int>{0, 1, 2} &&
            rows[1] == std::vector<int>{1, 2, 0},
        "integer rows read, '@' and '#' lines skipped" + (ok ? std::string() : ": " + error));
  write("tpt_time.dat", "0.000 0 1 2\n2.000 1 0 2\n");
  ok = postproc_detail::read_tpt_demux_map("tpt_time.dat", rows, error);
  check(!ok && names(error, "tpt_time.dat", 1),
        "replica_index.xvg row with a time column refused: " +
            (ok ? "read " + std::to_string(rows.size()) + " row(s) of " +
                      std::to_string(rows.empty() ? 0 : rows[0].size()) + " id(s)"
                : error));
  write("tpt_token.dat", "0 1 2\n0 1 x\n");
  ok = postproc_detail::read_tpt_demux_map("tpt_token.dat", rows, error);
  check(!ok && names(error, "tpt_token.dat", 2) && error.find("'x'") != std::string::npos,
        "unreadable token refused with file:line: " + (ok ? std::string("accepted") : error));

  std::cout << "case mbar demux map\n";
  titan_fe::DemuxStateMap map;
  write("mbar_token.dat", "0 1 2\n0 1 x\n");
  ok = map.read("mbar_token.dat", 3, 3, error);
  check(!ok && names(error, "mbar_token.dat", 2),
        "unreadable token refused with file:line: " + (ok ? std::string("accepted") : error));

  if (failures) {
    std::cout << "FAIL KIN-23 input parsers: " << failures << " check(s) failed\n";
    return 1;
  }
  std::cout << "PASS KIN-23 input parsers: comments skipped, other text refused "
               "with file and line\n";
  return 0;
}
