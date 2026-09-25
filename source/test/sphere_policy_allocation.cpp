#include "tabi/tabi_near.h"
#include "titan_executor.h"
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>

static std::atomic<bool> count_allocations{false};
static std::atomic<std::size_t> allocations{0};
void *operator new(std::size_t size) {
  if(count_allocations.load(std::memory_order_relaxed))
    allocations.fetch_add(1,std::memory_order_relaxed);
  if(void *p=std::malloc(size?size:1))return p;
  throw std::bad_alloc();
}
void operator delete(void *p) noexcept {std::free(p);}
void operator delete(void *p,std::size_t) noexcept {std::free(p);}

int main() {
  using namespace TitanPBE::Tabi;
  TitanExecutor executor(1);Options o;o.salt_m=0;o.surface_point_count=80;
  auto panels=build_surface({{0,0,0,1,1.7}},o,nullptr,&executor);
  const std::string near="analytic_singular",off="off";
  TreecodeStats stats;
  count_allocations=true;void *probe=::operator new(64);count_allocations=false;
  ::operator delete(probe);
  if(allocations!=1)return 2;
  // Warm the stats string and static quadrature tables before counting.
  for(const std::string &mode:std::array<std::string,2>{near,"centroid"}) {
    evaluate_jcp_near_pair_kernel(panels[0],0,panels[0],0,78.5,0,false,mode,off,&stats);
    // Centroid self is intentionally skipped, so warm a non-self pair too.
    evaluate_jcp_near_pair_kernel(panels[0],0,panels[1],1,78.5,0,false,mode,off,&stats);
    allocations=0;count_allocations=true;
    double sum=0;
    for(int repeat=0;repeat<10000;++repeat) {
      auto k=evaluate_jcp_near_pair_kernel(panels[0],0,panels[repeat%panels.size()],repeat%panels.size(),
                                          78.5,0,false,mode,off,&stats);
      sum+=k.k2;
    }
    count_allocations=false;
    std::cout<<mode<<" allocations="<<allocations<<" pairs=10000 checksum="<<sum<<'\n';
    if(allocations!=0)return 1;
  }
  return 0;
}
