// Equivalence harness for ENE-05: the grid-accelerated nearest-heavy-atom
// lookup in src/gbsa_intrin.cpp must return exactly what the previous linear
// scan returned, including the strict cutoff test and the lowest-index tie
// rule. The grid code below is a verbatim copy of GbGeometryIndex's search and
// build logic with the MOL/mol_info accessors replaced by plain vectors.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

namespace {

constexpr double kMaximumHydrogenBondLengthA = 1.55;

struct Mol {
  std::vector<double> x, y, z;
  std::vector<char> element; // 'H' or a heavy-atom symbol
  size_t size() const { return x.size(); }
};

// ---- reference: the implementation that shipped before this change ---------
char reference_nearest_heavy_element(const Mol &mol, size_t atom_index) {
  const double cutoff2 =
      kMaximumHydrogenBondLengthA * kMaximumHydrogenBondLengthA;
  double best2 = cutoff2;
  char best = '\0';
  for (size_t other = 0; other < mol.size(); ++other) {
    if (other == atom_index || mol.element[other] == 'H')
      continue;
    const double dx = mol.x[atom_index] - mol.x[other];
    const double dy = mol.y[atom_index] - mol.y[other];
    const double dz = mol.z[atom_index] - mol.z[other];
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (std::isfinite(r2) && r2 < best2) {
      best2 = r2;
      best = mol.element[other];
    }
  }
  return best;
}

// ---- candidate: GbGeometryIndex search/build, verbatim --------------------
class GbGeometryIndex {
public:
  explicit GbGeometryIndex(const Mol &mol) {
    count_ = mol.size();
    element_.assign(mol.element.begin(), mol.element.end());
    build_grid(mol);
  }

  char nearest_heavy_element(const Mol &mol, size_t atom_index) const {
    if (atom_index >= count_ || heavy_.empty())
      return '\0';
    const double qx = mol.x[atom_index];
    const double qy = mol.y[atom_index];
    const double qz = mol.z[atom_index];
    if (!std::isfinite(qx) || !std::isfinite(qy) || !std::isfinite(qz))
      return '\0';

    constexpr double kCutoff2 =
        kMaximumHydrogenBondLengthA * kMaximumHydrogenBondLengthA;
    double best2 = kCutoff2;
    size_t best = std::numeric_limits<size_t>::max();
    const auto consider = [&](size_t other) {
      if (other == atom_index)
        return;
      const double dx = qx - mol.x[other];
      const double dy = qy - mol.y[other];
      const double dz = qz - mol.z[other];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (!std::isfinite(r2) || !(r2 < kCutoff2))
        return;
      if (r2 < best2 || (r2 == best2 && other < best)) {
        best2 = r2;
        best = other;
      }
    };

    if (cells_valid_) {
      const size_t cx = cell_of(qx, min_x_, nx_);
      const size_t cy = cell_of(qy, min_y_, ny_);
      const size_t cz = cell_of(qz, min_z_, nz_);
      const size_t x0 = cx > 0 ? cx - 1 : 0;
      const size_t y0 = cy > 0 ? cy - 1 : 0;
      const size_t z0 = cz > 0 ? cz - 1 : 0;
      const size_t x1 = std::min(cx + 1, nx_ - 1);
      const size_t y1 = std::min(cy + 1, ny_ - 1);
      const size_t z1 = std::min(cz + 1, nz_ - 1);
      for (size_t ix = x0; ix <= x1; ++ix)
        for (size_t iy = y0; iy <= y1; ++iy)
          for (size_t iz = z0; iz <= z1; ++iz) {
            const size_t cell = (ix * ny_ + iy) * nz_ + iz;
            for (size_t slot = cell_start_[cell]; slot < cell_start_[cell + 1];
                 ++slot)
              consider(cell_items_[slot]);
          }
    } else {
      for (const size_t other : heavy_)
        consider(other);
    }
    return best == std::numeric_limits<size_t>::max() ? '\0' : element_[best];
  }

  bool gridded() const { return cells_valid_; }

private:
  size_t cell_of(double value, double origin, size_t extent) const {
    const double scaled = std::floor((value - origin) / cell_edge_);
    if (!(scaled > 0.0))
      return 0;
    if (scaled >= static_cast<double>(extent))
      return extent - 1;
    return static_cast<size_t>(scaled);
  }

  static size_t axis_extent(double span, double cell_edge) {
    if (!std::isfinite(span) || span <= 0.0)
      return 1;
    const double cells = std::floor(span / cell_edge) + 1.0;
    return cells >= 1.0 ? static_cast<size_t>(cells) : 1;
  }

  void build_grid(const Mol &mol) {
    heavy_.reserve(count_);
    for (size_t atom = 0; atom < count_; ++atom) {
      if (element_[atom] == 'H')
        continue;
      if (!std::isfinite(mol.x[atom]) || !std::isfinite(mol.y[atom]) ||
          !std::isfinite(mol.z[atom]))
        continue;
      heavy_.push_back(atom);
    }
    if (heavy_.empty())
      return;

    min_x_ = max_x_ = mol.x[heavy_.front()];
    min_y_ = max_y_ = mol.y[heavy_.front()];
    min_z_ = max_z_ = mol.z[heavy_.front()];
    for (const size_t atom : heavy_) {
      min_x_ = std::min(min_x_, mol.x[atom]);
      max_x_ = std::max(max_x_, mol.x[atom]);
      min_y_ = std::min(min_y_, mol.y[atom]);
      max_y_ = std::max(max_y_, mol.y[atom]);
      min_z_ = std::min(min_z_, mol.z[atom]);
      max_z_ = std::max(max_z_, mol.z[atom]);
    }

    constexpr double kMaxCells = 4.0e6;
    cell_edge_ = kMaximumHydrogenBondLengthA;
    for (int attempt = 0; attempt < 64; ++attempt) {
      nx_ = axis_extent(max_x_ - min_x_, cell_edge_);
      ny_ = axis_extent(max_y_ - min_y_, cell_edge_);
      nz_ = axis_extent(max_z_ - min_z_, cell_edge_);
      const double total = static_cast<double>(nx_) *
                           static_cast<double>(ny_) * static_cast<double>(nz_);
      if (total <= kMaxCells) {
        cells_valid_ = true;
        break;
      }
      cell_edge_ *= 2.0;
    }
    if (!cells_valid_)
      return;

    const size_t cells = nx_ * ny_ * nz_;
    cell_start_.assign(cells + 1, 0);
    for (const size_t atom : heavy_)
      ++cell_start_[cell_for(mol, atom) + 1];
    for (size_t cell = 0; cell < cells; ++cell)
      cell_start_[cell + 1] += cell_start_[cell];
    cell_items_.resize(heavy_.size());
    std::vector<size_t> cursor(cell_start_.begin(), cell_start_.end() - 1);
    for (const size_t atom : heavy_)
      cell_items_[cursor[cell_for(mol, atom)]++] = atom;
  }

  size_t cell_for(const Mol &mol, size_t atom) const {
    const size_t ix = cell_of(mol.x[atom], min_x_, nx_);
    const size_t iy = cell_of(mol.y[atom], min_y_, ny_);
    const size_t iz = cell_of(mol.z[atom], min_z_, nz_);
    return (ix * ny_ + iy) * nz_ + iz;
  }

  size_t count_{0};
  std::vector<char> element_;
  std::vector<size_t> heavy_;
  std::vector<size_t> cell_start_;
  std::vector<size_t> cell_items_;
  double min_x_{0.0}, max_x_{0.0};
  double min_y_{0.0}, max_y_{0.0};
  double min_z_{0.0}, max_z_{0.0};
  double cell_edge_{kMaximumHydrogenBondLengthA};
  size_t nx_{1}, ny_{1}, nz_{1};
  bool cells_valid_{false};
};

int failures = 0;
long long compared = 0;

void compare_case(const char *label, const Mol &mol) {
  const GbGeometryIndex index(mol);
  int local = 0;
  for (size_t atom = 0; atom < mol.size(); ++atom) {
    const char expected = reference_nearest_heavy_element(mol, atom);
    const char actual = index.nearest_heavy_element(mol, atom);
    ++compared;
    if (expected != actual) {
      if (local < 3)
        std::printf("  MISMATCH %s atom %zu: linear '%c' grid '%c'\n", label,
                    atom, expected ? expected : '0', actual ? actual : '0');
      ++local;
    }
  }
  if (local != 0) {
    failures += local;
    std::printf("FAIL %-34s %d/%zu atoms differ\n", label, local, mol.size());
  } else {
    std::printf("ok   %-34s %zu atoms, grid=%s\n", label, mol.size(),
                index.gridded() ? "yes" : "fallback");
  }
}

Mol random_mol(std::mt19937 &rng, size_t n, double extent, double h_fraction) {
  std::uniform_real_distribution<double> coord(0.0, extent);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  const char heavy_types[] = {'C', 'N', 'O', 'S'};
  Mol mol;
  for (size_t i = 0; i < n; ++i) {
    mol.x.push_back(coord(rng));
    mol.y.push_back(coord(rng));
    mol.z.push_back(coord(rng));
    mol.element.push_back(unit(rng) < h_fraction
                              ? 'H'
                              : heavy_types[rng() % 4]);
  }
  return mol;
}

// Dense clusters put many candidates inside the 1.55 A shell, which is where
// the tie rule and the 27-cell neighbourhood actually get exercised.
Mol clustered_mol(std::mt19937 &rng, size_t clusters, size_t per_cluster) {
  std::uniform_real_distribution<double> centre(0.0, 30.0);
  std::uniform_real_distribution<double> jitter(-1.6, 1.6);
  const char heavy_types[] = {'C', 'N', 'O', 'S'};
  Mol mol;
  for (size_t c = 0; c < clusters; ++c) {
    const double cx = centre(rng), cy = centre(rng), cz = centre(rng);
    for (size_t k = 0; k < per_cluster; ++k) {
      mol.x.push_back(cx + jitter(rng));
      mol.y.push_back(cy + jitter(rng));
      mol.z.push_back(cz + jitter(rng));
      mol.element.push_back(k % 3 == 0 ? 'H' : heavy_types[rng() % 4]);
    }
  }
  return mol;
}

} // namespace

int main() {
  std::mt19937 rng(20260913u);

  for (int trial = 0; trial < 12; ++trial) {
    char label[64];
    std::snprintf(label, sizeof(label), "random n=%d trial %d",
                  200 + trial * 137, trial);
    compare_case(label, random_mol(rng, 200 + trial * 137, 25.0, 0.5));
  }
  for (int trial = 0; trial < 8; ++trial) {
    char label[64];
    std::snprintf(label, sizeof(label), "clustered trial %d", trial);
    compare_case(label, clustered_mol(rng, 40, 12));
  }

  // Exact tie: two heavy atoms equidistant from the hydrogen. The linear scan
  // keeps the lower index, so the grid must too even though it may visit the
  // higher index first.
  {
    Mol mol;
    mol.x = {0.0, -1.0, 1.0};
    mol.y = {0.0, 0.0, 0.0};
    mol.z = {0.0, 0.0, 0.0};
    mol.element = {'H', 'N', 'O'};
    compare_case("exact tie lower index wins", mol);
    Mol swapped = mol;
    swapped.element = {'H', 'O', 'N'};
    compare_case("exact tie reversed elements", swapped);
  }

  // Distance exactly at the cutoff must be rejected by both.
  {
    Mol mol;
    mol.x = {0.0, kMaximumHydrogenBondLengthA, 1.40};
    mol.y = {0.0, 0.0, 0.0};
    mol.z = {0.0, 0.0, 0.0};
    mol.element = {'H', 'N', 'O'};
    compare_case("candidate exactly at cutoff", mol);
  }

  // Non-finite coordinates on the query and on candidates.
  {
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    const double inf_value = std::numeric_limits<double>::infinity();
    Mol mol;
    mol.x = {0.0, 1.0, nan_value, inf_value, 0.9};
    mol.y = {0.0, 0.0, 0.0, 0.0, 0.0};
    mol.z = {0.0, 0.0, 0.0, 0.0, 0.0};
    mol.element = {'H', 'C', 'N', 'O', 'S'};
    compare_case("non-finite candidates", mol);

    Mol bad_query = mol;
    bad_query.x[0] = nan_value;
    compare_case("non-finite query", bad_query);
  }

  // Hydrogen with no heavy atom in range, and an all-hydrogen molecule.
  {
    Mol mol;
    mol.x = {0.0, 8.0};
    mol.y = {0.0, 0.0};
    mol.z = {0.0, 0.0};
    mol.element = {'H', 'C'};
    compare_case("no heavy atom within cutoff", mol);

    Mol only_h;
    only_h.x = {0.0, 1.0, 2.0};
    only_h.y = {0.0, 0.0, 0.0};
    only_h.z = {0.0, 0.0, 0.0};
    only_h.element = {'H', 'H', 'H'};
    compare_case("all hydrogens", only_h);
  }

  // Query far outside the heavy-atom bounding box (cell clamping path).
  {
    Mol mol = random_mol(rng, 300, 20.0, 0.4);
    mol.x.push_back(5000.0);
    mol.y.push_back(-5000.0);
    mol.z.push_back(5000.0);
    mol.element.push_back('H');
    compare_case("query far outside bbox", mol);
  }

  // Very large extent forces the cell edge to grow past 1.55 A.
  {
    Mol mol = random_mol(rng, 400, 4.0e5, 0.5);
    compare_case("large extent grows cell edge", mol);
  }

  // Single heavy atom and single hydrogen degenerate grids.
  {
    Mol mol;
    mol.x = {0.0, 1.0};
    mol.y = {0.0, 0.2};
    mol.z = {0.0, 0.0};
    mol.element = {'H', 'C'};
    compare_case("two atoms", mol);
  }

  std::printf("\ncompared %lld atom lookups, %d mismatches\n", compared,
              failures);
  return failures == 0 ? 0 : 1;
}
