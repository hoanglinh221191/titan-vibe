#include "titan_tica_estimator.h"

extern "C" {
#include "xdrfile_xtc.h"
}

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_tica_estimator_smoke: " << message << '\n';
    std::exit(1);
  }
}

void check_xtc_header_scan() {
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "titan_tica_frame_scan.xtc";
  constexpr int natoms = 12;
  matrix box{};
  box[0][0] = 4.0f;
  box[1][1] = 4.0f;
  box[2][2] = 4.0f;
  std::vector<rvec> coords(static_cast<std::size_t>(natoms));

  XDRFILE *out = xdrfile_open(path.string().c_str(), "w");
  require(out != nullptr, "could not create XTC frame-scan fixture");
  for (int frame = 0; frame < 8; ++frame) {
    for (int atom = 0; atom < natoms; ++atom) {
      coords[static_cast<std::size_t>(atom)][0] =
          0.01f * static_cast<float>(atom + frame);
      coords[static_cast<std::size_t>(atom)][1] =
          0.02f * static_cast<float>(atom - frame);
      coords[static_cast<std::size_t>(atom)][2] =
          0.03f * static_cast<float>(atom + 2 * frame);
    }
    const int rc = write_xtc(out, natoms, frame * 50,
                             static_cast<float>(frame * 50), box,
                             coords.data(), 1000.0f);
    require(rc == exdrOK, "could not write XTC frame-scan fixture");
  }
  require(xdrfile_close(out) == exdrOK,
          "could not close XTC frame-scan fixture");

  TitanTicaFrameEstimateConfig scan;
  scan.dt_ps = 100.0;
  scan.xtc_paths = {path.string()};
  scan.fallback_estimated_frames = 1.0;
  const auto estimated = estimate_tica_selected_frames(scan);
  require(estimated.frames == 4,
          "XTC header scan should select four 100-ps frames");
  require(estimated.source == "xtc-header-sample",
          "XTC estimate should report its sampled-header source");
  require(estimated.diagnostic.empty(),
          "successful XTC scan should not report a diagnostic");

  std::error_code ec;
  std::filesystem::remove(path, ec);
}

} // namespace

int main() {
  TitanTicaFrameEstimateConfig remd48;
  remd48.begin_time_ps = 100000.0;
  remd48.end_time_ps = 827850.0;
  remd48.dt_ps = 200.0;
  remd48.begin_set = true;
  remd48.end_set = true;
  remd48.trajectory_count = 48;
  remd48.number_of_files = 48;
  remd48.fallback_estimated_frames = 316.0;

  const auto estimated = estimate_tica_selected_frames(remd48);
  require(estimated.frames == 174720,
          "48-replica selected-frame estimate should be 174720");
  require(estimated.source == "time-window",
          "48-replica estimate should report time-window source");

  TitanTicaFrameEstimateConfig fallback;
  fallback.dt_ps = 0.0;
  fallback.end_set = false;
  fallback.trajectory_count = 0;
  fallback.number_of_files = 0;
  fallback.fallback_estimated_frames = -12.0;

  const auto fallback_estimated = estimate_tica_selected_frames(fallback);
  require(fallback_estimated.frames == 1,
          "fallback estimate should never be less than one frame");
  require(fallback_estimated.source == "fallback",
          "fallback estimate should report fallback source");

  TitanTicaFrameEstimateConfig missing_xtc;
  missing_xtc.xtc_paths = {"definitely_missing_titan_tica_probe.xtc"};
  missing_xtc.fallback_estimated_frames = 27.0;
  const auto missing_estimated = estimate_tica_selected_frames(missing_xtc);
  require(missing_estimated.frames == 27,
          "failed XTC scan should use the fallback estimate");
  require(!missing_estimated.diagnostic.empty(),
          "failed XTC scan should report a diagnostic");

  check_xtc_header_scan();

  std::cout << "PASS titan_tica_estimator_smoke\n";
  return 0;
}
