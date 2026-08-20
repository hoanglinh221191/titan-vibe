#ifndef TITAN_VIBE_MULTI_MOL_H
#define TITAN_VIBE_MULTI_MOL_H

#include "MOL.h"
#include "MOL_info.h"
#include "MOL_para.h"
#include "arguments.h"
#include "index_v2.h"
#include "mol_topology.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

class TitanExecutor;

// Minimal system container required by the published VIBE adapter.
class Multi_MOL {
private:
  std::vector<unsigned char> pbctopo_frame_certified_;
  std::vector<std::string> pbctopo_frame_status_;

  int find_best_image_envelope(std::size_t frame_index, int chain_ref,
                               arguments &Arguments);
  int find_best_image_of_mol2(Multi_MOL &mol2, std::size_t frame_index,
                              int chain_ref, arguments &Arguments);

public:
  std::vector<MOL> m_mol;
  Multi_MOL_parameters parameters;
  std::vector<class index> Index;
  std::vector<mol_info> info;
  MolTopology topology;

  Multi_MOL() noexcept = default;
  explicit Multi_MOL(std::string_view, bool = false) noexcept {}

  [[nodiscard]] bool parameters_match_info() const noexcept {
    return parameters.size() == info.size();
  }

  void cache_topology_vdw_radii_from_parameters();

  [[nodiscard]] bool has_pbctopo_frame_status(std::size_t frame_index) const
      noexcept {
    return frame_index < pbctopo_frame_status_.size();
  }

  [[nodiscard]] bool pbctopo_frame_certified(std::size_t frame_index) const
      noexcept {
    return frame_index < pbctopo_frame_certified_.size() &&
           pbctopo_frame_certified_[frame_index] != 0;
  }

  [[nodiscard]] std::string_view
  pbctopo_frame_status(std::size_t frame_index) const noexcept {
    static constexpr std::string_view missing{"missing_status"};
    return frame_index < pbctopo_frame_status_.size()
               ? std::string_view(pbctopo_frame_status_[frame_index])
               : missing;
  }

  int solve_pbc_pbctopo(Multi_MOL &mol2, arguments &Arguments,
                        TitanExecutor *executor = nullptr,
                        std::size_t frame_report_offset = 0);
};

#endif
