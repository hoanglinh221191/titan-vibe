#include "tabi/tabi_near.h"
#include "titan_executor.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <array>
#include <vector>
#include <numbers>
using namespace TitanPBE::Tabi;

// Independent curved-triangle rule from physical vertices. It does not read
// the prepared production nodes; used for an actual cross-atom shared vertex.
double reference_cross_kernel(const Panel &target,const Panel &source,int order) {
  using V=std::array<double,3>;
  auto add=[](V a,V b){return V{a[0]+b[0],a[1]+b[1],a[2]+b[2]};};
  auto sub=[](V a,V b){return V{a[0]-b[0],a[1]-b[1],a[2]-b[2]};};
  auto scale=[](V a,double s){return V{a[0]*s,a[1]*s,a[2]*s};};
  auto dot=[](V a,V b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];};
  auto cross=[](V a,V b){return V{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};};
  const auto &sphere=sphere_patch_sphere(source);V origin{sphere[0],sphere[1],sphere[2]};
  std::array<V,3> vertices={sub(V{source.v0x,source.v0y,source.v0z},origin),
      sub(V{source.v1x,source.v1y,source.v1z},origin),sub(V{source.v2x,source.v2y,source.v2z},origin)};
  V c=scale(add(add(vertices[0],vertices[1]),vertices[2]),1./3.);
  V t=sub(V{target.x,target.y,target.z},origin);
  std::vector<std::array<double,2>> rule;
  for(int i=0;i<order;++i) {
    double x=std::cos(std::numbers::pi*(i+.75)/(order+.5)),derivative=0;
    for(int it=0;it<40;++it) {
      double p=1,previous=0;
      for(int k=1;k<=order;++k){double next=((2*k-1)*x*p-(k-1)*previous)/k;previous=p;p=next;}
      derivative=order*(x*p-previous)/(x*x-1);double dx=p/derivative;x-=dx;
      if(std::abs(dx)<1e-15)break;
    }
    rule.push_back({(1+x)/2,1/((1-x*x)*derivative*derivative)});
  }
  double integral=0;
  for(int e=0;e<3;++e)for(auto u:rule)for(auto v:rule) {
    V a=vertices[e],b=vertices[(e+1)%3];
    V du=sub(add(scale(a,1-v[0]),scale(b,v[0])),c),p=add(c,scale(du,u[0]));
    double length=std::sqrt(dot(p,p));V normal=scale(p,1/length),y=scale(normal,sphere[3]),d=sub(t,y);
    double distance=std::sqrt(dot(d,d));
    double jac=sphere[3]*sphere[3]*u[0]*std::abs(dot(p,cross(du,sub(b,a))))/(length*length*length);
    integral+=u[1]*v[1]*jac*dot(normal,d)/(4*std::numbers::pi*distance*distance*distance);
  }
  return 77.5*integral;
}

int main() {
  try {
    int checks=0;
    auto require=[&](bool condition,const char *message) {
      ++checks;if(!condition)throw std::runtime_error(message);
    };
    TitanExecutor executor(2);
    Options o;o.salt_m=0;o.surface_mode="vdw";o.sphere_mesh_mode="icosphere";
    o.surface_point_count=320;o.near_integration="auto";
    std::vector<Atom> atoms={{0,0,0,1,1.7},{10,0,0,-1,2.3}};
    auto panels=build_surface(atoms,o,nullptr,&executor);
    require(panels.size()==640,"two disjoint sphere meshes");
    require(resolve_near_integration_mode("auto",panels)=="analytic_singular","auto resolution");
    std::unordered_set<const SpherePatchQuadrature*> unique;
    for(const auto &p:panels) unique.insert(p.sphere_patch_quadrature.get());
    require(unique.size()==320 && !unique.count(nullptr),"unit geometry shared across atoms");
    require(panels[0].sphere_patch_quadrature==panels[320].sphere_patch_quadrature,"scaled atom shares template");
    require(sphere_patch_sphere(panels[0])[3]==1.7 && sphere_patch_sphere(panels[320])[3]==2.3,"atom-index geometry lookup");
    auto self=evaluate_jcp_near_pair_kernel(panels[0],0,panels[0],0,78.5,0,false,"analytic_singular","off");
    require(self.valid && self.k2<0 && self.used_duffy,"curved self is finite and nonzero");
    auto scaled=evaluate_jcp_near_pair_kernel(panels[320],320,panels[320],320,78.5,0,false,"analytic_singular","off");
    require(self.k2==scaled.k2,"self integral is scale invariant");
    auto centroid=evaluate_jcp_near_pair_kernel(panels[0],0,panels[0],0,78.5,0,false,"centroid","off");
    require(is_intentional_zero_self_kernel(centroid,false,"centroid"),"explicit centroid legacy self");
    for(const auto &mode: {"centroid","analytic_singular","duffy_adaptive"}) {
      for(std::size_t j: {std::size_t(0),std::size_t(1),std::size_t(30)}) {
        auto target=panels[0],source=panels[j];
        auto tagged=evaluate_jcp_near_pair_kernel(target,0,source,j,78.5,.1,true,mode,"off");
        target.sphere_patch_quadrature.reset();source.sphere_patch_quadrature.reset();
        auto plain=evaluate_jcp_near_pair_kernel(target,0,source,j,78.5,.1,true,mode,"off");
        require(tagged.valid==plain.valid && tagged.k1==plain.k1 && tagged.k2==plain.k2 &&
                tagged.k3==plain.k3 && tagged.k4==plain.k4,"coupled kernel ignores scalar metadata");
      }
    }
    for(int mode=0;mode<5;++mode) {
      auto excluded=o;
      if(mode==0)excluded.surface_mode="sas";
      if(mode==1)excluded.salt_m=.15;
      if(mode==2)excluded.surface_only=true;
      if(mode==3)excluded.surface_geometry="point-panel";
      if(mode==4)excluded.near_integration="centroid";
      auto mesh=build_surface(atoms,excluded,nullptr,&executor);
      bool clean=true;for(auto &p:mesh)clean=clean&&!p.sphere_patch_quadrature;
      require(clean,"out-of-scope paths have no new quadrature");
    }
    auto touching=build_surface({{0,0,0,0,1.7},{3.4,0,0,0,1.7}},o,nullptr,&executor);
    bool found=false;
    for(std::size_t i=0;i<touching.size()&&!found;++i)
      for(std::size_t j=i+1;j<touching.size()&&!found;++j) {
        if(touching[i].atom_index==touching[j].atom_index)continue;
        auto k=evaluate_jcp_near_pair_kernel(touching[i],i,touching[j],j,78.5,0,false,"analytic_singular","off");
        if(!k.one_ring)continue;
        found=true;
        const double r32=reference_cross_kernel(touching[i],touching[j],32);
        const double r64=reference_cross_kernel(touching[i],touching[j],64);
        require(std::abs(r32-r64)<1e-9,"cross-sphere reference convergence");
        require(k.valid && k.used_duffy && std::abs(k.k2-r64)<1e-5*std::max(1.,std::abs(r64)),"cross-sphere branch agrees with independent quadrature");
        auto collision=touching[i];const auto &source=touching[j];
        const auto &s=sphere_patch_sphere(source);auto node=source.sphere_patch_quadrature->nodes.front();
        collision.x=s[0]+s[3]*node[0];collision.y=s[1]+s[3]*node[1];collision.z=s[2]+s[3]*node[2];
        TreecodeStats stats;
        auto fallback=evaluate_jcp_near_pair_kernel(collision,i,source,j,78.5,0,false,"analytic_singular","off",&stats);
        auto legacy=evaluate_jcp_near_pair_kernel(collision,i,source,j,78.5,0,false,"centroid","off");
        require(fallback.valid && fallback.used_centroid && stats.near_fallbacks==1 && fallback.k2==legacy.k2,"collision falls back with an audit count");
      }
    require(found,"real shared vertex reaches cross-sphere branch");
    // Test the legacy latlong cap's most elongated polar triangle as well as
    // ordinary panels; the adaptive self rule must remain finite there.
    o.sphere_mesh_mode="latlong";o.surface_point_count=3969;
    auto fine=build_surface({atoms[0]},o,nullptr,&executor);
    require(fine.size()==16128,"legacy latlong cap");
    require(std::isfinite(fine[0].sphere_patch_quadrature->self_double_layer),"polar self quadrature converges");
    std::cout<<"PASS "<<checks<<" kernel/compatibility checks; sizeof(Panel)="<<sizeof(Panel)<<"\n";
    return 0;
  } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
