#include "vibe_only/absent.h"

#include <cstdio>
#include <cstdlib>

namespace titan_vibe_only {

void absent(const char *entry_point) noexcept {
  std::fflush(stdout);
  std::fprintf(stderr,
               "TITAN-VIBE internal error: %s belongs to a TITAN module that is not part of this "
               "build. Please report the command line that led here.\n",
               entry_point);
  std::fflush(stderr);
  std::abort();
}

} // namespace titan_vibe_only
