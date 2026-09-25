// Synthetic writer fixture: uses the SAME helper as the native status writer.
#include "titan_path/status_numeric_io.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

struct CommaDecimal : std::numpunct<char> {
  char do_decimal_point() const override { return ','; }
};

int main(int argc, char **argv) {
  try {
    if (argc != 4)
      throw std::invalid_argument("expected three finite vector components");
    std::array<double, 3> vector{};
    for (int i = 0; i < 3; ++i) {
      std::size_t end = 0;
      const std::string value(argv[i + 1]);
      vector[i] = std::stod(value, &end);
      if (end != value.size() || !std::isfinite(vector[i]))
        throw std::invalid_argument("invalid component");
    }
    // Deliberately inherit a non-machine locale and low precision: the real
    // writer helper must neutralize both rather than rely on process defaults.
    std::cout.imbue(std::locale(std::locale::classic(), new CommaDecimal));
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "SCHEMA=TITAN_PATH_STATUS_V1\nSTATUS=CONTINUE\n";
    titan_path::write_next_vector_status(
        std::cout, vector, std::hypot(vector[0], vector[1], vector[2]), 0.);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 2;
  }
}
