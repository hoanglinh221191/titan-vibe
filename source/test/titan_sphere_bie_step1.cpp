// Step 1: production sphere solver versus an independent Kirkwood oracle.
// No production kernel, diagonal, or energy is modified by this driver.
#include "tabi/tabi_types.h"
#include "tabi/tabi_near.h"
#include "titan_executor.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

using namespace TitanPBE::Tabi;
constexpr double pi = std::numbers::pi;
constexpr double coulomb = 332.063713299;
struct V { double x, y, z; };
V operator+(V a,V b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
V operator-(V a,V b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
V operator*(V a,double b) { return {a.x*b,a.y*b,a.z*b}; }
double dot(V a,V b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
V cross(V a,V b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
double norm(V a) { return std::sqrt(dot(a,a)); }

// Kirkwood boundary trace in the solver's 1/(4*pi) normalization.
double exact_trace(double cosine,double offset,double radius,double q,
                   double ei,double eo) {
  long double sum=0, pprev=1, p=1, power=1;
  for (int l=0;l<400;++l) {
    if(l==1) p=cosine;
    else if(l>1) {long double next=((2*l-1)*cosine*p-(l-1)*pprev)/l; pprev=p;p=next;}
    sum+=(2*l+1)*power*p/(ei*l+eo*(l+1));
    power*=offset;
  }
  return q*static_cast<double>(sum)/(4*pi*radius);
}
double exact_energy(double offset,double radius,double q,double ei,double eo) {
  long double sum=0, power=1;
  for(int l=0;l<400;++l) {
    sum+=(ei-eo)*(l+1)*power/(ei*l+eo*(l+1));
    power*=offset*offset;
  }
  return 0.5*coulomb*q*q*static_cast<double>(sum)/(ei*radius);
}

// Screened (linearised PB) Kirkwood energy, for salt > 0.
// Matching phi and eps*dphi/dr on the sphere gives the exterior response
// gamma_l = kappa*a*k_l'(kappa*a)/k_l(kappa*a), with k_l the modified spherical
// Bessel function of the second kind.  Since k_l'(z) = -k_{l+1}(z)+(l/z)k_l(z),
// gamma_l = l - z*r_l with r_l = k_{l+1}(z)/k_l(z).  That ratio obeys
// r_{l+1} = (2l+3)/z + 1/r_l, r_0 = 1 + 1/z, which stays bounded where the
// functions themselves overflow.  As z -> 0, gamma_l -> -(l+1) and the series
// reduces term by term to the unscreened exact_energy above.
double exact_energy_screened(double offset,double radius,double q,double ei,
                             double eo,double kappa) {
  if(kappa<=0) return exact_energy(offset,radius,q,ei,eo);
  const long double z=kappa*radius;
  long double sum=0, power=1, r=1+1/z;
  for(int l=0;l<400;++l) {
    const long double gamma=l-z*r;
    sum+=(eo*gamma+ei*(l+1))*power/(ei*l-eo*gamma);
    power*=offset*offset;
    r=(2.0L*l+3.0L)/z+1/r;
  }
  return 0.5*coulomb*q*q*static_cast<double>(sum)/(ei*radius);
}

struct Rule {std::vector<double> x,w;};
Rule gauss(int n) {
  Rule r;
  for(int i=0;i<n;++i) {
    double z=std::cos(pi*(i+.75)/(n+.5)), d=0;
    for(int it=0;it<30;++it) {
      double p=1, prev=0;
      for(int k=1;k<=n;++k) {double old=p;p=((2*k-1)*z*p-(k-1)*prev)/k;prev=old;}
      d=n*(z*p-prev)/(z*z-1);
      double next=z-p/d;
      if(std::abs(next-z)<1e-15) {z=next;break;}
      z=next;
    }
    r.x.push_back((z+1)/2);r.w.push_back(1/((1-z*z)*d*d));
  }
  return r;
}

// Independent curved-sphere integration: radial map of a chord triangle,
// three Duffy fans from its barycenter. Jacobian is derived from the radial
// map, not the production flat-triangle near integration implementation.
double curved_k2(const Panel &target,const Panel &source,V origin,double radius,
                 double ratio,const Rule &rule) {
  V vertices[3]={{source.v0x,source.v0y,source.v0z},
                 {source.v1x,source.v1y,source.v1z},
                 {source.v2x,source.v2y,source.v2z}};
  for(auto &v:vertices) v=v-origin;
  V c=(vertices[0]+vertices[1]+vertices[2])*(1./3);
  V t=V{target.x,target.y,target.z}-origin;
  double total=0;
  for(int k=0;k<3;++k) {
    V a=vertices[k], b=vertices[(k+1)%3];
    for(std::size_t iu=0;iu<rule.x.size();++iu)
      for(std::size_t iv=0;iv<rule.x.size();++iv) {
        double u=rule.x[iu],v=rule.x[iv];
        V edge=a*(1-v)+b*v, du=edge-c, dv=(b-a)*u;
        V p=c+du*u;double len=norm(p);V normal=p*(1/len), y=normal*radius;
        double jac=radius*radius*std::abs(dot(p,cross(du,dv)))/(len*len*len);
        // For two points on the same sphere, n_y.(x-y)=-|x-y|^2/(2R).
        // This equivalent expression avoids cancellation near the self point.
        double kernel=-(ratio-1)/(8*pi*radius*norm(t-y));
        total+=rule.w[iu]*rule.w[iv]*jac*kernel;
      }
  }
  return total;
}

bool neighbors(const Panel&a,const Panel&b) {
  std::size_t av[3]={a.v0_index,a.v1_index,a.v2_index};
  std::size_t bv[3]={b.v0_index,b.v1_index,b.v2_index};
  for(auto i:av) for(auto j:bv) if(i==j)return true;
  return false;
}

int main(int argc,char**argv) {
  try {
    if(argc<7) throw std::runtime_error("mesh points offset near diagonal tol [radius q epsin epsout shift direction]");
    Options o;o.primary_surface="sphere-patch";o.surface_mode="vdw";
    o.sphere_mesh_mode=argv[1];o.surface_point_count=std::stoi(argv[2]);
    double offset=std::stod(argv[3]);o.near_integration=argv[4];
    o.diagonal_mode=argv[5];o.residual_tol=std::stod(argv[6]);
    double radius=argc>7?std::stod(argv[7]):1.7, q=argc>8?std::stod(argv[8]):1;
    o.eps_in=argc>9?std::stod(argv[9]):1;o.eps_out=argc>10?std::stod(argv[10]):78.5;
    double shift=argc>11?std::stod(argv[11]):0;int direction=argc>12?std::stoi(argv[12]):0;
    if(offset<0 || offset>.8 || radius<=0 || q==0 || o.eps_in<=0 || o.eps_out<o.eps_in)
      throw std::runtime_error("Unsupported oracle domain");
    V origin{shift,2*shift,3*shift}, axis=direction?V{1/std::sqrt(3.),1/std::sqrt(3.),1/std::sqrt(3.)}:V{1,0,0};
    V charge=origin+axis*(offset*radius);
    // Neutral enclosing sphere fixes the boundary; the tiny charge carrier is
    // fully buried even at offset=.8. Assert that it contributes no panels.
    std::vector<Atom> atoms={{origin.x,origin.y,origin.z,0,radius},
                             {charge.x,charge.y,charge.z,q,.01*radius}};
    // Discretization and salt come from the environment so every existing
    // positional invocation reproduces bit-identically.  argv[13] is already
    // the treecode switch, so extra positional slots are not available.
    const char *env_disc=std::getenv("TITAN_STEP1_DISCRETIZATION");
    const char *env_salt=std::getenv("TITAN_STEP1_SALT");
    const std::string discretization=env_disc?env_disc:"centroid";
    if(discretization!="centroid" && discretization!="node-patch")
      throw std::runtime_error("discretization must be centroid or node-patch");
    o.salt_m=env_salt?std::stod(env_salt):0;
    // Curved node-patch pair-cache budget, so the dense and the bounded-sparse
    // regimes can be compared at one fixed node count.
    if(const char *env_cache=std::getenv("TITAN_STEP1_CACHE_MB"))
      o.curved_pair_cache_mb=std::stoi(env_cache);
    // Surface selection.  The isolated sphere is the one geometry where every
    // surface builder should agree, so it can validate a closed SES mesh
    // against the same Kirkwood oracle used for sphere-patch.  An SES of one
    // atom is its vdW sphere: the probe finds nothing to roll into.
    if(const char *env_surface=std::getenv("TITAN_STEP1_SURFACE"))
      o.primary_surface=env_surface;
    if(const char *env_engine=std::getenv("TITAN_STEP1_SES_ENGINE"))
      o.analytic_ses_engine=env_engine;
    if(const char *env_spacing=std::getenv("TITAN_STEP1_SES_SPACING"))
      o.ses_grid_spacing=std::stod(env_spacing);
    if(const char *env_probe=std::getenv("TITAN_STEP1_PROBE"))
      o.probe_radius=std::stod(env_probe);
    const bool sphere_patch_surface=o.primary_surface=="sphere-patch";
    o.formulation="jcp2013";o.discretization=discretization;
    o.energy_projection="boundary";o.sum_mode="direct";o.max_iter=300;
    o.preconditioner="diagonal";o.quad_buffer="off";o.highres_surface=o.surface_point_count>4096;
    TitanExecutor executor(4);Result geometry;
    const auto panels=build_surface(atoms,o,&geometry,&executor);
    const auto n=panels.size();double area=0, radial_error=0;std::size_t carrier_tagged=0;
    std::vector<double> phi(n);double projection=0;
    const double ratio=o.eps_out/o.eps_in;
    for(std::size_t i=0;i<n;++i) {
      const auto &p=panels[i];
      // sphere-patch tags every panel with its parent atom, so a carrier panel
      // would mean the carrier broke the boundary.  SES builders attribute
      // panels differently and are only reported, not asserted.
      if(p.atom_index!=0) {
        if(sphere_patch_surface) throw std::runtime_error("carrier changed boundary");
        ++carrier_tagged;
      }
      V at=V{p.x,p.y,p.z}-origin;radial_error=std::max(radial_error,std::abs(norm(at)-radius));
      phi[i]=exact_trace(dot(at,axis)/radius,offset,radius,q,o.eps_in,o.eps_out);
      area+=p.area;V dr=charge-V{p.x,p.y,p.z};double d=norm(dr);
      projection+=p.area*(ratio-1)*dot(V{p.nx,p.ny,p.nz},dr)/(4*pi*d*d*d)*phi[i];
    }
    projection*=.5*coulomb*q*4*pi;
    // Sample eight deterministic rows. Report RMS normalized by sampled RHS;
    // this is a diagnostic, not the full-system solver residual.
    Rule r8=gauss(8),r64=gauss(64),r128=gauss(128);
    double err0=0,errself=0,errremoved=0,errcurved=0,bnorm=0,selfsize=0,quad_diff=0;
    double constant_error=0;std::size_t invalid_self=0;
    // The sampler's independent reference integrates a spherical TRIANGLE read
    // from panel corners.  A node-patch barycentric dual cell has no such
    // corners, and its diagonal and coupled unknowns differ as well, so the
    // sampler is left to centroid.  Node-patch is judged end to end against
    // the Kirkwood energy, which needs no panel-shape model at all.
    const double kappa_a=o.salt_m>0?0.3291*std::sqrt(o.salt_m)*std::sqrt(298.15/o.temperature_k):0.0;
    // The sampler also assumes the unscreened trace, so salt > 0 turns it off.
    const bool sampled=discretization=="centroid" && o.salt_m<=0 && sphere_patch_surface;
    for(std::size_t sample=0;sampled && sample<8;++sample) {
      std::size_t i=sample*n/8;const auto &p=panels[i];double prod=0,curved=0,kconstant=0;
      for(std::size_t j=0;j<n;++j) {
        auto kernel=evaluate_jcp_near_pair_kernel(p,i,panels[j],j,ratio,0,false,o.near_integration,"off");
        if(kernel.valid) prod+=kernel.k2*phi[j];
        else if(i==j && is_intentional_zero_self_kernel(kernel,false,o.near_integration)) ++invalid_self;
        else throw std::runtime_error("invalid non-self kernel");
        double k=curved_k2(p,panels[j],origin,radius,ratio,neighbors(p,panels[j])?r64:r8);
        curved+=k*phi[j];kconstant+=k;
      }
      double kself=curved_k2(p,p,origin,radius,ratio,r64);
      double kself128=curved_k2(p,p,origin,radius,ratio,r128);
      quad_diff=std::max(quad_diff,std::abs(kself-kself128));
      double rhs=q/(4*pi*o.eps_in*norm(V{p.x,p.y,p.z}-charge));
      double diagonal=.5*(1+ratio)*phi[i];
      err0+=std::pow(diagonal-prod-rhs,2);
      errself+=std::pow(diagonal-prod-kself*phi[i]-rhs,2);
      errremoved+=std::pow(diagonal-prod+kself*phi[i]-rhs,2);
      errcurved+=std::pow(diagonal-curved-rhs,2);
      bnorm+=rhs*rhs;selfsize+=std::pow(kself*phi[i],2);
      constant_error=std::max(constant_error,std::abs(kconstant+.5*(ratio-1)));
    }
    if(!sampled) bnorm=1; // keep the unmeasured ratios at 0 instead of NaN
    const auto begin=std::chrono::steady_clock::now();
    // The curved node-patch representations support direct and the classic
    // treecode only; the jcp2013-treecode formulation enables FMM M2L by
    // default, which the solver rejects for them, so pin it off here.
    if(argc>13) {o.sum_mode="treecode";o.formulation="jcp2013-treecode";o.fmm_m2l="off";
      o.tree_opening_theta=std::stod(argv[13]);o.tree_opening_theta_explicit=true;}
    auto result=solve_linear_pb(atoms,o,executor);
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    if(!result.converged || !std::isfinite(result.energy_kcal_mol)) throw std::runtime_error("solver failed");
    std::cout<<std::setprecision(17)<<"RESULT {"
      <<"\"discretization\":\""<<discretization<<"\",\"salt_m\":"<<o.salt_m
      <<",\"surface\":\""<<o.primary_surface<<"\",\"ses_engine\":\""<<o.analytic_ses_engine
      <<"\",\"ses_spacing\":"<<o.ses_grid_spacing<<",\"probe\":"<<o.probe_radius
      <<",\"carrier_tagged_panels\":"<<carrier_tagged
      <<",\"node_patch_representation\":\""<<result.node_patch_representation<<"\""
      <<",\"sampled\":"<<(sampled?"true":"false")<<",\"kappa_inv_a\":"<<kappa_a<<",\"cache_mb\":"<<o.curved_pair_cache_mb
      <<",\"panels\":"<<n<<",\"area\":"<<area<<",\"radial_error\":"<<radial_error
      <<",\"energy\":"<<result.energy_kcal_mol<<",\"reference\":"<<exact_energy_screened(offset,radius,q,o.eps_in,o.eps_out,kappa_a)
      <<",\"exact_trace_projection_energy\":"<<projection
      <<",\"sample_operator_rel\":"<<std::sqrt(err0/bnorm)
      <<",\"sample_self_added_rel\":"<<std::sqrt(errself/bnorm)
      <<",\"sample_self_removed_rel\":"<<std::sqrt(errremoved/bnorm)
      <<",\"sample_curved_quadrature_rel\":"<<std::sqrt(errcurved/bnorm)
      <<",\"sample_self_size_rel\":"<<std::sqrt(selfsize/bnorm)
      <<",\"self_quadrature_64_128_max_abs\":"<<quad_diff
      <<",\"curved_constant_identity_max_abs\":"<<constant_error
      <<",\"zero_self_samples\":"<<invalid_self
      <<",\"residual\":"<<result.true_residual_rel<<",\"iterations\":"<<result.iterations
      <<",\"solver_seconds\":"<<seconds<<",\"near_duffy_pairs\":"<<result.near_duffy_pairs
      <<",\"near_centroid_pairs\":"<<result.near_centroid_pairs
      <<",\"reported_near_policy\":\""<<result.near_kernel_policy<<"\""
      <<",\"memory_bytes\":"<<result.memory_bytes
      <<",\"quadrature_bytes\":"<<sphere_patch_quadrature_storage_bytes(panels)
      <<",\"matvec_mode\":\""<<result.matvec_mode<<"\""
      <<",\"projection_ms\":"<<result.energy_projection_ms
      <<",\"projection_kernel\":\""<<result.energy_projection_kernel<<"\""
      <<",\"charge_resolution_status\":\""<<result.charge_resolution.status<<"\""
      <<",\"min_d_over_h\":"<<result.charge_resolution.min_ratio
      <<",\"charge_resolution_ms\":"<<result.charge_resolution.milliseconds
      <<",\"energy_source\":\""<<result.reaction_energy_source<<"\"}\n";
    return 0;
  } catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
}
