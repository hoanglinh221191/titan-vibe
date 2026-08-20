#ifndef TITAN_VIBE_MOL_H
#define TITAN_VIBE_MOL_H

#include "aligned_allocator.h"
#include "box_type.h"

#include <cstddef>
#include <string_view>

struct AlignedVector3D {
  AlignedDoubleVector x;
  AlignedDoubleVector y;
  AlignedDoubleVector z;

  void resize(std::size_t size) {
    x.resize(size, 0.0);
    y.resize(size, 0.0);
    z.resize(size, 0.0);
  }

  void clear() noexcept {
    x.clear();
    y.clear();
    z.clear();
  }

  [[nodiscard]] std::size_t size() const noexcept { return x.size(); }
};

// Minimal frame representation required by VIBE.
class MOL {
public:
  AlignedVector3D pos;
  box Box;
  double time_of_snap{0.0};

  MOL() noexcept = default;
  MOL(const box &box_value, double time, std::string_view,
      bool = false) noexcept
      : Box(box_value), time_of_snap(time) {}

  [[nodiscard]] bool coordinates_consistent() const noexcept {
    return pos.x.size() == pos.y.size() && pos.x.size() == pos.z.size();
  }

  [[nodiscard]] bool topology_consistent(std::size_t atom_count) const
      noexcept {
    return coordinates_consistent() && pos.size() == atom_count;
  }

  [[nodiscard]] std::size_t size() const noexcept { return pos.size(); }
};

#endif
