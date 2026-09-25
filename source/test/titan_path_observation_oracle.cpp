// Read-only oracle: existing C++ geometry, not a candidate search or MD runner.
#include "titan_path/registration.h"
#include "titan_path/geometry.h"
#include "titan_path/continuity.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
int main() {
  std::size_t na, nw;
  while (std::cin >> na >> nw) {
    if (na < 3 || na > 10000 || nw < 2 || nw > 10000) return 2;
    std::vector<titan_path::Vec3> ref(na), cur(na), points(nw);
    for (auto& p:ref) for (auto& x:p) std::cin >> x;
    for (auto& p:cur) for (auto& x:p) std::cin >> x;
    for (auto& p:points) for (auto& x:p) std::cin >> x;
    titan_path::Vec3 com;
    for (auto& x:com) std::cin >> x;
    double previous, window;
    if (!(std::cin >> previous >> window)) return 3;
    const auto fit=titan_path::fit_rigid_transform(ref,cur);
    if (!fit.valid) return 4;
    for (auto& p:points) p=titan_path::apply_transform(fit,p);
    std::vector<titan_path::ArcLengthSegment> segments;
    double arc=0;
    for (std::size_t i=1;i<nw;++i) {
      titan_path::Vec3 d;
      for (int k=0;k<3;++k) d[k]=points[i][k]-points[i-1][k];
      const double length=titan_path::norm(d);
      segments.push_back({arc,length,points[i-1],d}); arc+=length;
    }
    auto global=titan_path::project_point_to_polyline(segments,com,previous);
    std::cout<<std::setprecision(17)<<global.progress_angstrom<<" "
             <<global.distance_angstrom<<" "<<fit.rmsd_angstrom<<" "
             <<global.segment_index<<" ";
    try {
      const auto local=titan_path::project_with_continuity(segments,com,previous,window,arc);
      std::cout<<"VALID "<<local.progress_angstrom<<" "<<local.distance_angstrom<<" "<<local.segment_index<<"\n";
    } catch (const std::invalid_argument& e) {
      if (std::string(e.what()).find("PROJECTION_CONTINUITY_WINDOW_EXHAUSTED")==std::string::npos) return 5;
      std::cout<<"CLIPPED 0 0 0\n";
    }
  }
  return std::cin.eof()?0:6;
}
