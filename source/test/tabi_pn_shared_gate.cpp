// Regression tests for allowed vertex/edge contact after the curved PN lift.
#include "../src/tabi/tabi_node_patch_curved.cpp"
#include <iomanip>
using namespace TitanPBE::Tabi;
namespace {
int checks = 0;
void require(bool ok, const char *name) {
  ++checks;
  if (!ok) throw std::runtime_error(name);
}
PnPatch plane(V a, V b, V c) {
  const V normal{0, 0, 1};
  return pn_patch(a, b, c, normal, normal, normal, 7);
}
void direct_pairs() {
  auto a = plane({0,0,0}, {3,0,0}, {0,3,0});
  auto b = plane({3,0,0}, {0,0,0}, {3,-3,0});
  require(pn_adjacent_patches_apart(a, {0,1,2}, b, {1,0,3}, true),
          "valid common edge rejected");
  auto same_side = plane({3,0,0}, {0,0,0}, {1,1,0});
  require(!pn_adjacent_patches_apart(a, {0,1,2}, same_side, {1,0,3}, true),
          "common-edge overlap accepted");
  auto v = plane({0,0,0}, {-3,0,0}, {0,-3,0});
  require(pn_adjacent_patches_apart(a, {0,1,2}, v, {0,3,4}, true),
          "valid vertex-only contact rejected");
  auto contained = plane({0,0,0}, {1,0,0}, {0,1,0});
  require(!pn_adjacent_patches_apart(a, {0,1,2}, contained, {0,3,4}, true),
          "shared-vertex containment accepted");
  auto broken = b;
  broken.b[4][2] = 1e-8;
  require(!pn_adjacent_patches_apart(a, {0,1,2}, broken, {1,0,3}, false),
          "different common cubic accepted");
  require(!pn_adjacent_patches_apart(a, {0,0,2}, b, {1,0,3}, true),
          "duplicate logical vertex accepted");
  for (double shift : {0.,10000.}) for (double scale : {.001,10.}) {
    auto x=a, y=b;
    for (auto *p : {&x,&y}) for (auto &v : p->b) {
      const auto old=v;
      v={shift+scale*(.6*old[0]-.8*old[1]),
         -shift+scale*(.8*old[0]+.6*old[1]), shift};
    }
    require(pn_adjacent_patches_apart(x,{0,1,2},y,{1,0,3},true),
            "translation/rotation/scale changed valid contact");
  }
}
void shared_vertex_repair() {
  NodePatchGeometry g;
  std::vector<V> normals;
  for (int k=0;k<3;++k) {
    const double a=2*pi*k/3;
    g.vertex_xyz.push_back({std::cos(a),std::sin(a),0});
    normals.push_back({.6*std::cos(a),.6*std::sin(a),.8});
  }
  for (int k=1;k<3;++k) {
    const double a=2*pi*k/3;
    g.vertex_xyz.push_back({std::cos(a),std::sin(a),.1});
    normals.push_back({.6*std::cos(a),.6*std::sin(a),-.8});
  }
  NodePatchTriangle a{}, b{};
  a.vertices={0,1,2}; b.vertices={0,4,3}; a.nz=1;
  V n=cross(sub(g.vertex_xyz[4],g.vertex_xyz[0]),sub(g.vertex_xyz[3],g.vertex_xyz[0]));
  n=mul(n,1/norm(n));b.nx=n[0];b.ny=n[1];b.nz=n[2];g.triangles={a,b};
  require(!pn_adjacent_patches_apart(pn_triangle_patch(normals,g,0,0),a.vertices,
                                     pn_triangle_patch(normals,g,1,7),b.vertices,false),
          "review shared-vertex crossing accepted");
  // The lift mends the pair through its vertex normals (one per vertex, so
  // the patches keep a common tangent plane at vertex 0) or refuses it; it
  // never flattens.  The star normal of a vertex is its area-weighted face
  // normal, as build_node_patch_panels makes it.
  std::vector<V> star(g.vertex_xyz.size(), V{});
  for (const auto &t : g.triangles) {
    const double area = norm(cross(sub(g.vertex_xyz[t.vertices[1]], g.vertex_xyz[t.vertices[0]]),
                                   sub(g.vertex_xyz[t.vertices[2]], g.vertex_xyz[t.vertices[0]])));
    for (auto v : t.vertices)
      star[v] = add(star[v], mul(V{t.nx, t.ny, t.nz}, area));
  }
  for (auto &s : star)
    s = mul(s, 1 / norm(s));
  PnRepairReport report;
  const auto residual = pn_repair_normals(normals, star, g, report);
  std::vector<std::uint8_t> mask(2, 0U);
  if (!residual.empty()) {
    PnLiftReport lift;
    mask = pn_flatten_residual(normals, g, residual, lift);
  }
  require(pn_adjacent_patches_apart(pn_triangle_patch(normals,g,0,mask[0]),a.vertices,
                                    pn_triangle_patch(normals,g,1,mask[1]),b.vertices,
                                    mask[0]==7 && mask[1]==7),
          "lifted pair has unqualified contact");
  std::cout<<"SHARED_VERTEX "<<(residual.empty()?"mended":"flattened")<<" star_normals="<<report.star
           <<" turned="<<report.turned<<" masks="<<int(mask[0])<<","<<int(mask[1])<<"\n";
}
void unresolved_paths() {
  NodePatchGeometry g;g.vertex_xyz={{0,0,0},{1,0,0},{0,1,0}};
  NodePatchTriangle t{};t.vertices={0,1,2};t.nz=1;g.triangles={t};
  g.incident_offsets={0,1,2,3};g.incident_triangles={0,0,0};
  CurvedNodePatchGeometry c;c.pn_node_normals={{0,0,1},{0,0,1},{0,0,1}};
  c.pn_straight_edges={0};c.cells.resize(3);c.cells[0].radius=.5;
  c.cells[0].samples=pn_samples(c.pn_node_normals,c.pn_straight_edges,g,0,2);
  Panel p;p.x=.2;p.y=.1;p.z=1e-12;p.nz=1;
  for(double zone:{0.,6.}) {
    c.one_point_beyond_radii=zone;bool refused=false;
    try { (void)geometric_pair(p,g,c,0,80,0); }
    catch(const std::runtime_error&e) { refused=std::string(e.what())=="curved_near_quadrature_unresolved"; }
    require(refused,"unresolved pair returned on a solver route");
  }
}
}
int main(){try{direct_pairs();shared_vertex_repair();unresolved_paths();
 std::cout<<"PASS "<<checks<<" shared-contact and unresolved-route checks\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<" after "<<checks<<" checks\n";return 1;}}
