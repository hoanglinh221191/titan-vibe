#include "titan_path/cartesian_audit.h"
#include <iostream>
#include <sstream>
#include <vector>
struct Info {
  std::size_t id;
  bool valid = true;
  bool has_original_atom_id() const { return valid; }
  std::size_t get_original_atom_id() const { return id; }
};
struct Positions {
  std::vector<double> x, y, z;
};
struct Mol {
  Positions pos;
};
struct Params {
  std::vector<double> mass;
};
struct Multi {
  std::vector<Info> info;
  std::vector<Mol> m_mol;
  Params parameters;
};
int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::filesystem::path dir(argv[1]);
  if (std::filesystem::exists(dir))
    return 3;
  std::filesystem::create_directories(dir);
  Multi a{{{0}}, {{{{1}, {2}, {3}}}}, {{12}}},
      b{{{1}}, {{{{4}, {5}, {6}}}}, {{1}}};
  int n = 0;
  auto check = [&](bool ok) {
    ++n;
    if (!ok)
      throw std::runtime_error("case " + std::to_string(n));
  };
  auto rejected = [&](const Multi &x, const std::filesystem::path &p) {
    try {
      titan_path::write_cartesian_audit(p, x, b);
      return false;
    } catch (const std::exception &) {
      return true;
    }
  };
  try {
    auto path = dir / "audit.tsv";
    titan_path::write_cartesian_audit(path, a, b);
    std::ifstream in(path);
    std::stringstream data;
    data << in.rdbuf();
    const auto original = data.str();
    check(original.find("mol1\t0\t0\t1\t2\t3\t12") != std::string::npos);
    check(original.find("mol2\t0\t1\t4\t5\t6\t1") != std::string::npos);
    check(rejected(a, path));
    Multi x = a;
    x.parameters.mass.clear();
    check(rejected(x, dir / "missing.tsv"));
    x = a;
    x.parameters.mass[0] = 0;
    check(rejected(x, dir / "zero.tsv"));
    x = a;
    x.m_mol[0].pos.x[0] = std::numeric_limits<double>::infinity();
    check(rejected(x, dir / "infinite.tsv"));
    x = a;
    x.info[0].valid = false;
    check(rejected(x, dir / "identity.tsv"));
    check(!std::filesystem::exists(dir / "missing.tsv") &&
          !std::filesystem::exists(dir / "zero.tsv"));
    std::cout << "TITAN_CARTESIAN_AUDIT_SMOKE_PASS cases=" << n << "\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
