#include "xdrfile_xtc.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "usage: titan_xtc_scan_probe <trajectory.xtc> [max_frames]\n";
    return 2;
  }

  const std::string path = argv[1];
  const long long max_frames = argc == 3 ? std::atoll(argv[2]) : -1;
  int natoms = 0;
  const int natoms_rc = read_xtc_natoms(argv[1], &natoms);
  if (natoms_rc != exdrOK || natoms <= 0) {
    std::cerr << "read_xtc_natoms failed: rc=" << natoms_rc
              << " natoms=" << natoms << '\n';
    return 1;
  }

  XDRFILE *handle = xdrfile_open(path.c_str(), "r");
  if (handle == nullptr) {
    std::cerr << "xdrfile_open failed\n";
    return 1;
  }

  std::vector<rvec> coordinates(static_cast<std::size_t>(natoms));
  matrix box{};
  int step = 0;
  float time_ps = 0.0f;
  float precision = 0.0f;
  long long frames = 0;
  int rc = exdrOK;
  while (max_frames < 0 || frames < max_frames) {
    rc = read_xtc(handle, natoms, &step, &time_ps, box, coordinates.data(),
                  &precision);
    if (rc == exdrENDOFFILE)
      break;
    if (rc != exdrOK) {
      std::cerr << "decode failed at frame=" << frames << " rc=" << rc
                << " step=" << step << " time_ps=" << time_ps << '\n';
      xdrfile_close(handle);
      return 1;
    }
    ++frames;
  }
  const int close_rc = xdrfile_close(handle);
  if (close_rc != exdrOK) {
    std::cerr << "xdrfile_close failed: rc=" << close_rc << '\n';
    return 1;
  }

  std::cout << "OK XTC scan frames=" << frames << " natoms=" << natoms
            << " last_step=" << step << " last_time_ps=" << time_ps << '\n';
  return 0;
}
