// Unit gate for the -cscd distance module (src/impl_cscd_class.cpp).
//
// Two chains, A as molecule 1 and B as molecule 2, against brute-force
// distances: nearest and -act layouts, parallel = serial bit for bit, frames
// spread over several engine chunks keep their own rows, replicas average at
// equal times while segments append, the self pair is skipped for -g=x-x,
// and a residue list that changes between frames is refused.
#include "cscd.h"
#include "titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

int failures = 0;

void check(bool ok, const std::string &what) {
  if (!ok) {
    ++failures;
    std::fprintf(stderr, "FAIL %s\n", what.c_str());
  }
}

std::vector<com> chain(const std::string &id, std::size_t residues,
                       std::mt19937_64 &rng, double shift) {
  std::uniform_real_distribution<double> coord(-12.0, 12.0);
  std::vector<com> out;
  for (std::size_t r = 0; r < residues; ++r)
    out.emplace_back(coord(rng) + shift, coord(rng), coord(rng), id, "LYS",
                     r + 1, 0, 0, r + 1);
  return out;
}

double distance(const com &a, const com &b) {
  const double dx = a.get_x() - b.get_x();
  const double dy = a.get_y() - b.get_y();
  const double dz = a.get_z() - b.get_z();
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

struct Table {
  std::vector<std::string> comments;
  std::vector<std::vector<double>> rows; // time, then the distances
};

Table read_table(const fs::path &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("missing " + path.string());
  Table table;
  for (std::string line; std::getline(in, line);) {
    if (!line.empty() && line[0] == '#') {
      table.comments.push_back(line);
      continue;
    }
    std::istringstream fields(line);
    std::vector<double> row;
    for (std::string token; fields >> token;)
      row.push_back(std::strtod(token.c_str(), nullptr));
    if (!row.empty())
      table.rows.push_back(row);
  }
  return table;
}

bool near(double a, double b) {
  if (std::isnan(a) || std::isnan(b))
    return std::isnan(a) && std::isnan(b);
  return std::abs(a - b) <= 6e-5; // the files keep 4 decimals
}

all_cscd configured(bool nearest, bool skip_self) {
  all_cscd cscd;
  all_cscd::Options options;
  options.nearest_only = nearest;
  options.skip_same_residue = skip_self;
  cscd.configure(options);
  return cscd;
}

fs::path fresh(const fs::path &dir) {
  fs::remove_all(dir);
  fs::create_directories(dir);
  return dir;
}

// Nearest and -act against brute force, over frames that the engine would
// have delivered in three chunks (local index 0 again at each chunk start).
void nearest_and_all(const fs::path &work) {
  std::mt19937_64 rng(20261012);
  const std::size_t frames = 7;
  std::vector<std::vector<com>> a, b;
  for (std::size_t f = 0; f < frames; ++f) {
    a.push_back(chain("A", 5, rng, 0.0));
    b.push_back(chain("B", 4, rng, 6.0));
  }
  for (const bool nearest : {true, false}) {
    all_cscd cscd = configured(nearest, false);
    for (std::size_t f = 0; f < frames; ++f)
      cscd.add_frame(f, 100.0 + 10.0 * f, a[f], b[f], nullptr, 1);
    const fs::path dir = fresh(work / (nearest ? "nearest" : "all"));
    check(cscd.write(dir.string()) == 5, "five residue files");
    check(cscd.columns() == (nearest ? 5u : 20u), "column count");
    for (std::size_t r = 0; r < 5; ++r) {
      const Table t = read_table(dir / ("out_cal_dis_com_A" +
                                        std::to_string(r + 1) + ".dat"));
      check(t.rows.size() == frames, "one row per frame, chunks included");
      for (std::size_t f = 0; f < t.rows.size(); ++f) {
        const auto &row = t.rows[f];
        check(near(row[0], 100.0 + 10.0 * f), "row time");
        if (nearest) {
          double best = std::numeric_limits<double>::infinity();
          for (const com &partner : b[f])
            best = std::min(best, distance(a[f][r], partner));
          check(row.size() == 2 && near(row[1], best), "nearest distance");
        } else {
          check(row.size() == 5, "-act: one column per molecule-2 residue");
          for (std::size_t j = 0; j < 4 && j + 1 < row.size(); ++j)
            check(near(row[j + 1], distance(a[f][r], b[f][j])),
                  "-act distance");
        }
      }
      if (!nearest)
        check(!t.comments.empty() &&
                  t.comments.back().find("B4") != std::string::npos,
              "-act header names the partners");
    }
  }
}

// The parallel split must not change a single bit.
void parallel_equals_serial(const fs::path &work) {
  std::mt19937_64 rng(7);
  const std::vector<com> a = chain("A", 3000, rng, 0.0);
  const std::vector<com> b = chain("B", 60, rng, 3.0);
  TitanExecutor executor(4);
  for (const bool nearest : {true, false}) {
    all_cscd serial = configured(nearest, false);
    all_cscd parallel = configured(nearest, false);
    serial.add_frame(0, 0.0, a, b, nullptr, 1);
    parallel.add_frame(0, 0.0, a, b, &executor, 4);
    const fs::path s = fresh(work / "serial");
    const fs::path p = fresh(work / "parallel");
    serial.write(s.string());
    parallel.write(p.string());
    for (const std::size_t r : {std::size_t{0}, std::size_t{1499},
                                std::size_t{2999}}) {
      const std::string name =
          "out_cal_dis_com_A" + std::to_string(r + 1) + ".dat";
      std::ifstream fs_serial(s / name), fs_parallel(p / name);
      std::stringstream ts, tp;
      ts << fs_serial.rdbuf();
      tp << fs_parallel.rdbuf();
      check(!ts.str().empty() && ts.str() == tp.str(),
            "parallel output equals serial output");
    }
  }
}

// Replicas (same times) average; a later segment (new times) appends.
void merge_by_time(const fs::path &work) {
  auto place = [](double x) {
    return std::vector<com>{com(x, 0.0, 0.0, "A", "ASP", 1, 0, 0, 1)};
  };
  const std::vector<com> partner{com(0.0, 0.0, 0.0, "B", "ARG", 1, 0, 0, 1)};

  all_cscd cscd = configured(true, false);
  const fs::path traj1 = fresh(work / "traj1");
  // Replica 1 at 0, 10, 20 ps; distances 1, 2, 3.
  for (std::size_t f = 0; f < 3; ++f)
    cscd.add_frame(f, 10.0 * f, place(1.0 + f), partner, nullptr, 1);
  cscd.end_trajectory(traj1.string());
  // Replica 2 at the same times; distances 3, 4, 5.
  for (std::size_t f = 0; f < 3; ++f)
    cscd.add_frame(f, 10.0 * f, place(3.0 + f), partner, nullptr, 1);
  cscd.end_trajectory({});
  // A segment at 30 and 40 ps; distances 7, 8.
  for (std::size_t f = 0; f < 2; ++f)
    cscd.add_frame(f, 30.0 + 10.0 * f, place(7.0 + f), partner, nullptr, 1);
  cscd.end_trajectory({});

  const fs::path merged = fresh(work / "merged");
  cscd.write(merged.string());
  const Table own = read_table(traj1 / "out_cal_dis_com_A1.dat");
  check(own.rows.size() == 3 && near(own.rows[2][1], 3.0),
        "per-trajectory file keeps its own rows");
  const Table t = read_table(merged / "out_cal_dis_com_A1.dat");
  const double expect_time[] = {0, 10, 20, 30, 40};
  const double expect_dist[] = {2, 3, 4, 7, 8};
  check(t.rows.size() == 5, "replicas merge, the segment appends");
  for (std::size_t k = 0; k < t.rows.size() && k < 5; ++k)
    check(near(t.rows[k][0], expect_time[k]) &&
              near(t.rows[k][1], expect_dist[k]),
          "merged row " + std::to_string(k));
  check(cscd.trajectories() == 3 && cscd.averaged_rows() == 3,
        "merge statistics");
}

// -g=x-x: a residue is not its own partner.
void self_pair_skipped(const fs::path &work) {
  std::mt19937_64 rng(5);
  const std::vector<com> a = chain("A", 4, rng, 0.0);
  for (const bool nearest : {true, false}) {
    all_cscd cscd = configured(nearest, true);
    cscd.add_frame(0, 0.0, a, a, nullptr, 1);
    const fs::path dir = fresh(work / (nearest ? "self_nearest" : "self_all"));
    cscd.write(dir.string());
    for (std::size_t r = 0; r < 4; ++r) {
      const Table t = read_table(dir / ("out_cal_dis_com_A" +
                                        std::to_string(r + 1) + ".dat"));
      std::vector<double> expected;
      for (std::size_t j = 0; j < 4; ++j)
        if (j != r)
          expected.push_back(distance(a[r], a[j]));
      if (nearest)
        expected = {*std::min_element(expected.begin(), expected.end())};
      check(t.rows.size() == 1 && t.rows[0].size() == expected.size() + 1,
            "self pair omitted");
      for (std::size_t c = 0; c < expected.size() && c + 1 < t.rows[0].size();
           ++c)
        check(near(t.rows[0][c + 1], expected[c]), "self-skip distance");
    }
  }
  // A residue whose only partner is itself has no nearest: NaN, not 0.
  const std::vector<com> one{com(1.0, 2.0, 3.0, "A", "GLY", 1, 0, 0, 1)};
  all_cscd lonely = configured(true, true);
  lonely.add_frame(0, 0.0, one, one, nullptr, 1);
  const fs::path dir = fresh(work / "lonely");
  lonely.write(dir.string());
  const Table t = read_table(dir / "out_cal_dis_com_A1.dat");
  check(t.rows.size() == 1 && t.rows[0].size() == 2 && std::isnan(t.rows[0][1]),
        "no partner gives NaN");
}

void guards(const fs::path &work) {
  std::mt19937_64 rng(3);
  const std::vector<com> a = chain("A", 3, rng, 0.0);
  const std::vector<com> b = chain("B", 3, rng, 0.0);
  all_cscd cscd = configured(true, false);
  cscd.add_frame(0, 0.0, a, b, nullptr, 1);
  bool threw = false;
  try {
    cscd.add_frame(1, 1.0, a, chain("B", 2, rng, 0.0), nullptr, 1);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  check(threw, "a changed residue list is refused");
  threw = false;
  try {
    cscd.add_frame(0, 2.0, a, b, nullptr, 1);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  check(threw, "a repeated frame number is refused");

  // Two residues with the same blank chain and number (atoms 1 and 12; the
  // engine's internal residue numbers 6 and 9) still get a file each: named
  // by internal residue number, or by first atom when that is unknown.
  std::vector<com> twins{com(0, 0, 0, " ", "ALA", 7, 0, 0, 1),
                         com(5, 0, 0, " ", "ALA", 7, 0, 1, 12)};
  const std::vector<std::size_t> twin_internal{6, 9};
  for (const bool numbered : {true, false}) {
    all_cscd dup = configured(true, false);
    if (numbered) {
      dup = all_cscd{};
      dup.configure(all_cscd::Options{}, twin_internal, {});
    }
    dup.add_frame(0, 0.0, twins, b, nullptr, 1);
    const fs::path dir = fresh(work / (numbered ? "dup_internal" : "dup_atom"));
    const char *first = numbered ? "out_cal_dis_com_7_r7.dat"
                                 : "out_cal_dis_com_7_a1.dat";
    const char *second = numbered ? "out_cal_dis_com_7_r10.dat"
                                  : "out_cal_dis_com_7_a12.dat";
    check(dup.write(dir.string()) == 2 && fs::exists(dir / first) &&
              fs::exists(dir / second),
          std::string("repeated chain+number gets unique names: ") + first);
  }

  // -g=x-x with the twins: each is skipped only against itself, so its
  // partner is the other twin (5 A), and the -act columns are told apart.
  for (const bool nearest : {true, false}) {
    all_cscd self;
    all_cscd::Options options;
    options.nearest_only = nearest;
    options.skip_same_residue = true;
    self.configure(options, twin_internal, twin_internal);
    self.add_frame(0, 0.0, twins, twins, nullptr, 1);
    const fs::path out = fresh(work / (nearest ? "twins_nearest" : "twins_all"));
    self.write(out.string());
    for (const char *name :
         {"out_cal_dis_com_7_r7.dat", "out_cal_dis_com_7_r10.dat"}) {
      const Table t = read_table(out / name);
      check(t.rows.size() == 1 && t.rows[0].size() == 2 &&
                near(t.rows[0][1], 5.0),
            std::string("twin skips only itself: ") + name);
    }
  }
  all_cscd act;
  all_cscd::Options options;
  options.nearest_only = false;
  act.configure(options, {}, twin_internal);
  act.add_frame(0, 0.0, b, twins, nullptr, 1);
  const fs::path cols = fresh(work / "twin_columns");
  act.write(cols.string());
  const Table t = read_table(cols / "out_cal_dis_com_B1.dat");
  check(!t.comments.empty() &&
            t.comments.back().find(" 7_r7") != std::string::npos &&
            t.comments.back().find(" 7_r10") != std::string::npos,
        "-act header tells repeated partner labels apart");

  // A centre that moves to another first atom is another residue.
  all_cscd moved = configured(true, false);
  moved.add_frame(0, 0.0, twins, b, nullptr, 1);
  std::vector<com> shifted = twins;
  shifted[1].set_respos(13);
  threw = false;
  try {
    moved.add_frame(1, 1.0, shifted, b, nullptr, 1);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  check(threw, "a centre with another first atom is refused");
}

} // namespace

int main(int argc, char **argv) {
  const fs::path work =
      argc > 1 ? fs::path(argv[1]) : fs::temp_directory_path() / "cscd_smoke";
  try {
    fs::create_directories(work);
    nearest_and_all(work);
    parallel_equals_serial(work);
    merge_by_time(work);
    self_pair_skipped(work);
    guards(work);
  } catch (const std::exception &e) {
    std::fprintf(stderr, "FAIL exception: %s\n", e.what());
    return 1;
  }
  if (failures != 0) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("PASS cscd_module_smoke\n");
  return 0;
}
