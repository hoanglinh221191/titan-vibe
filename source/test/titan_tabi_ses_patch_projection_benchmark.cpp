#include "tabi/tabi_analytic_ses.h"
#include "tabi/tabi_ses_patch_intersections.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::analytic_ses;

std::vector<Atom> read_pqr(const std::string &path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open PQR: " + path);
  std::vector<Atom> atoms;
  std::string line;
  while (std::getline(input, line)) {
    if (!(line.starts_with("ATOM") || line.starts_with("HETATM")))
      continue;
    std::istringstream stream(line);
    std::vector<std::string> fields;
    std::string field;
    while (stream >> field)
      fields.push_back(field);
    if (fields.size() < 9U)
      continue;
    const std::size_t count = fields.size();
    Atom atom;
    atom.x = std::stod(fields[count - 5U]);
    atom.y = std::stod(fields[count - 4U]);
    atom.z = std::stod(fields[count - 3U]);
    atom.q = std::stod(fields[count - 2U]);
    atom.radius = std::stod(fields[count - 1U]);
    atoms.push_back(atom);
  }
  if (atoms.empty())
    throw std::runtime_error("PQR contains no atoms: " + path);
  return atoms;
}

struct Query {
  std::size_t atom{0U};
  Vec3 point;
  std::size_t contact_patch{0U};
};

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2 || argc > 3 ||
        (argc == 3 && std::string(argv[2]) != "--ab-only")) {
      std::cerr << "usage: patch_projection_benchmark <input.pqr> "
                   "[--ab-only]\n";
      return EXIT_FAILURE;
    }
    const bool ab_only = argc == 3;
    const std::vector<Atom> atoms = read_pqr(argv[1]);
    const auto patch_begin = std::chrono::steady_clock::now();
    const PatchComplex complex =
        build_patch_complex(atoms, 1.4, nullptr, false);
    const auto patch_end = std::chrono::steady_clock::now();
    if (complex.status != "pass")
      throw std::runtime_error("patch complex status=" + complex.status);
    const LocalPatchProjectionIndex index =
        build_local_patch_projection_index(complex);
    const auto index_end = std::chrono::steady_clock::now();
    if (!index.bounded_domains_certified ||
        !index.bvh_invariants_certified) {
      const auto constraint_diagnostics =
          [](const std::vector<std::size_t> &offsets,
             const std::vector<SpherePatchConstraint> &constraints) {
            std::size_t empty_constraints = 0U;
            std::size_t near_parallel_pairs = 0U;
            double maximum_minimum_dot =
                -std::numeric_limits<double>::infinity();
            for (const SpherePatchConstraint &constraint : constraints) {
              maximum_minimum_dot =
                  std::max(maximum_minimum_dot,
                           constraint.minimum_dot);
              if (constraint.minimum_dot > 1.0 + 4.0e-10)
                ++empty_constraints;
            }
            if (!offsets.empty()) {
              for (std::size_t patch = 0;
                   patch + 1U < offsets.size(); ++patch) {
                for (std::size_t first = offsets[patch];
                     first < offsets[patch + 1U]; ++first) {
                  for (std::size_t second = first + 1U;
                       second < offsets[patch + 1U]; ++second) {
                    const Vec3 &a = constraints[first].normal;
                    const Vec3 &b = constraints[second].normal;
                    const Vec3 cross{
                        a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x};
                    const double cross2 =
                        cross.x * cross.x + cross.y * cross.y +
                        cross.z * cross.z;
                    if (cross2 > 0.0 && cross2 <= 1.0e-14)
                      ++near_parallel_pairs;
                  }
                }
              }
            }
            std::ostringstream output;
            output << "empty_constraints=" << empty_constraints
                   << ",near_parallel_pairs=" << near_parallel_pairs
                   << ",max_minimum_dot=" << maximum_minimum_dot;
            return output.str();
          };
      const ActivePairCacheFailure &failure =
          index.active_pair_cache_failure;
      std::ostringstream failure_diagnostics;
      failure_diagnostics
          << "reason=" << static_cast<int>(failure.reason)
          << ",kind=" << static_cast<int>(failure.kind)
          << ",patch=" << failure.patch
          << ",first=" << failure.first
          << ",second=" << failure.second
          << ",normal2=(" << failure.first_normal2 << ','
          << failure.second_normal2 << ')'
          << ",minimum_dot=(" << failure.first_minimum_dot << ','
          << failure.second_minimum_dot << ')'
          << ",cross2=" << failure.cross_squared
          << ",base2=" << failure.base_squared
          << ",active_residual=("
          << failure.first_active_residual << ','
          << failure.second_active_residual << ')'
          << ",feasibility_margin="
          << failure.minimum_feasibility_margin
          << ",candidate_feasible=" << failure.candidate_feasible;
      throw std::runtime_error(
          "prepared patch index is not closest-point certified; "
          "bvh=" +
          std::to_string(index.bvh_invariants_certified) +
          "; contact_cache_offsets=" +
          std::to_string(index.contact_two_active_offsets.size()) +
          "/" + std::to_string(complex.contacts.size() + 1U) +
          "; reentrant_cache_offsets=" +
          std::to_string(index.reentrant_two_active_offsets.size()) +
          "/" + std::to_string(complex.reentrant.size() + 1U) +
          "; contact{" +
          constraint_diagnostics(index.contact_constraint_offsets,
                                 index.contact_constraints) +
          "}; reentrant{" +
          constraint_diagnostics(index.reentrant_constraint_offsets,
                                 index.reentrant_constraints) +
          "}; cache_failure{" + failure_diagnostics.str() + "}");
    }

    std::size_t maximum_constraints = 0U;
    std::size_t total_constraints = 0U;
    std::vector<std::size_t> constraint_counts;
    constraint_counts.reserve(complex.contacts.size());
    for (std::size_t patch = 0; patch < complex.contacts.size(); ++patch) {
      const std::size_t count =
          index.contact_constraint_offsets[patch + 1U] -
          index.contact_constraint_offsets[patch];
      constraint_counts.push_back(count);
      maximum_constraints = std::max(maximum_constraints, count);
      total_constraints += count;
    }
    std::sort(constraint_counts.begin(), constraint_counts.end());
    const auto percentile = [&](double fraction) {
      if (constraint_counts.empty())
        return std::size_t{0U};
      const double rank =
          std::ceil(fraction *
                    static_cast<double>(constraint_counts.size()));
      const std::size_t position = static_cast<std::size_t>(
          std::max(1.0, rank) - 1.0);
      return constraint_counts[
          std::min(position, constraint_counts.size() - 1U)];
    };
    const std::size_t constraints_p50 = percentile(0.50);
    const std::size_t constraints_p95 = percentile(0.95);
    const std::size_t constraints_p99 = percentile(0.99);
    const double mean_constraints =
        complex.contacts.empty()
            ? 0.0
            : static_cast<double>(total_constraints) /
                  static_cast<double>(complex.contacts.size());
    struct CountDistribution {
      std::size_t total{0U};
      std::size_t p50{0U};
      std::size_t p95{0U};
      std::size_t p99{0U};
      std::size_t maximum{0U};
    };
    const auto count_distribution =
        [](const std::vector<std::size_t> &offsets,
           std::size_t patch_count) {
          if (offsets.size() != patch_count + 1U)
            throw std::runtime_error(
                "active-pair cache offsets are malformed");
          std::vector<std::size_t> counts;
          counts.reserve(patch_count);
          for (std::size_t patch = 0; patch < patch_count; ++patch) {
            if (offsets[patch] > offsets[patch + 1U])
              throw std::runtime_error(
                  "active-pair cache offsets are not monotone");
            counts.push_back(
                offsets[patch + 1U] - offsets[patch]);
          }
          std::sort(counts.begin(), counts.end());
          const auto at = [&](double fraction) {
            if (counts.empty())
              return std::size_t{0U};
            const std::size_t rank = static_cast<std::size_t>(
                std::max(
                    1.0,
                    std::ceil(fraction *
                              static_cast<double>(counts.size()))) -
                1.0);
            return counts[std::min(rank, counts.size() - 1U)];
          };
          CountDistribution result;
          result.total = offsets.empty() ? 0U : offsets.back();
          result.p50 = at(0.50);
          result.p95 = at(0.95);
          result.p99 = at(0.99);
          result.maximum = counts.empty() ? 0U : counts.back();
          return result;
        };
    const CountDistribution contact_cache =
        count_distribution(index.contact_two_active_offsets,
                           complex.contacts.size());
    const CountDistribution reentrant_cache =
        count_distribution(index.reentrant_two_active_offsets,
                           complex.reentrant.size());
    if (contact_cache.total !=
            index.contact_two_active_directions.size() ||
        reentrant_cache.total !=
            index.reentrant_two_active_directions.size())
      throw std::runtime_error(
          "active-pair cache totals do not match storage");
    const std::size_t empty_contact_domains =
        static_cast<std::size_t>(std::count_if(
            index.contact_domain_empty.begin(),
            index.contact_domain_empty.end(),
            [](std::uint8_t value) { return value != 0U; }));
    const std::size_t empty_reentrant_domains =
        static_cast<std::size_t>(std::count_if(
            index.reentrant_domain_empty.begin(),
            index.reentrant_domain_empty.end(),
            [](std::uint8_t value) { return value != 0U; }));
    const double patch_build_ms =
        std::chrono::duration<double, std::milli>(
            patch_end - patch_begin)
            .count();
    const double index_build_ms =
        std::chrono::duration<double, std::milli>(
            index_end - patch_end)
            .count();

    constexpr std::size_t kDirectionsPerAtom = 4U;
    const double kGoldenAngle =
        std::numbers::pi * (3.0 - std::sqrt(5.0));
    std::vector<Query> queries;
    queries.reserve(complex.contacts.size() * kDirectionsPerAtom);
    for (std::size_t patch = 0; patch < complex.contacts.size(); ++patch) {
      const ContactPatch &contact = complex.contacts[patch];
      for (std::size_t sample = 0; sample < kDirectionsPerAtom; ++sample) {
        const std::size_t serial =
            patch * kDirectionsPerAtom + sample;
        const double z =
            1.0 - 2.0 * (static_cast<double>(sample) + 0.5) /
                      static_cast<double>(kDirectionsPerAtom);
        const double radial =
            std::sqrt(std::max(0.0, 1.0 - z * z));
        const double angle =
            kGoldenAngle * static_cast<double>(serial);
        const Vec3 direction{radial * std::cos(angle),
                             radial * std::sin(angle), z};
        const double query_radius = contact.atom_radius + 0.05;
        queries.push_back(
            {contact.atom,
             {contact.center.x + query_radius * direction.x,
              contact.center.y + query_radius * direction.y,
              contact.center.z + query_radius * direction.z},
             patch});
      }
    }

    struct ContactRun {
      double elapsed_us{0.0};
      std::size_t hits{0U};
      std::size_t certified{0U};
      double distance_sum{0.0};
    };
    const auto run_contact = [&]() {
      ContactRun run;
      const auto begin = std::chrono::steady_clock::now();
      for (const Query &query : queries) {
        const LocalProjectionHit hit =
            project_patch_bounded_closest(
                complex, index, PatchKind::contact,
                query.contact_patch, query.point, 1.0);
        run.hits += hit.found ? 1U : 0U;
        run.certified +=
            hit.bounded_closest_certified ? 1U : 0U;
        if (hit.found)
          run.distance_sum += hit.distance;
      }
      run.elapsed_us =
          std::chrono::duration<double, std::micro>(
              std::chrono::steady_clock::now() - begin)
              .count();
      return run;
    };

    struct GlobalRun {
      double elapsed_us{0.0};
      std::size_t hits{0U};
      std::size_t certified{0U};
      std::size_t owner_hint_hits{0U};
      std::size_t owner_hint_candidates{0U};
      std::size_t spatial_fallbacks{0U};
      std::size_t axis_fallbacks{0U};
      std::size_t candidates{0U};
      std::array<std::size_t, 3> candidates_by_kind{};
      std::array<std::size_t, 3> empty_short_circuits_by_kind{};
      std::array<std::uint64_t, 3> projector_nanoseconds_by_kind{};
      std::size_t pruned{0U};
      std::size_t nodes_visited{0U};
      std::size_t actual_stamp_bytes{0U};
      double distance_sum{0.0};
      double selection_checksum{0.0};
    };
    struct SelectedProjection {
      Vec3 point;
      double distance{0.0};
      std::size_t patch_index{0U};
      PatchKind kind{PatchKind::contact};
      bool found{false};
      bool certified{false};
    };
    const auto run_global =
        [&](const PatchComplex &run_complex,
            const LocalPatchProjectionIndex &run_index,
            std::span<const Query> run_queries,
            std::vector<SelectedProjection> *selections,
            ProjectionHintPolicy hint_policy,
            bool collect_kind_timing) {
      GlobalRun run;
      if (selections != nullptr) {
        selections->clear();
        selections->reserve(run_queries.size());
      }
      const auto begin = std::chrono::steady_clock::now();
      for (const Query &query : run_queries) {
        const LocalProjectionHit hit =
            project_local_patch_closest_with_hint_policy_for_audit(
                run_complex, run_index, query.atom, query.point, 1.0,
                hint_policy, collect_kind_timing);
        run.hits += hit.found ? 1U : 0U;
        run.certified +=
            hit.bounded_closest_certified ? 1U : 0U;
        run.owner_hint_hits += hit.owner_hint_hit ? 1U : 0U;
        run.owner_hint_candidates += hit.owner_hint_candidates;
        run.spatial_fallbacks += hit.spatial_fallback ? 1U : 0U;
        run.axis_fallbacks += hit.axis_fallback ? 1U : 0U;
        run.candidates += hit.candidate_patches;
        for (std::size_t kind = 0; kind < 3U; ++kind) {
          run.candidates_by_kind[kind] +=
              hit.candidate_patches_by_kind[kind];
          run.empty_short_circuits_by_kind[kind] +=
              hit.empty_domain_short_circuits_by_kind[kind];
          run.projector_nanoseconds_by_kind[kind] +=
              hit.projector_nanoseconds_by_kind[kind];
        }
        run.pruned += hit.bound_pruned_patches;
        run.nodes_visited += hit.spatial_nodes_visited;
        run.actual_stamp_bytes =
            std::max(run.actual_stamp_bytes,
                     hit.stamp_bytes_per_worker);
        if (hit.found)
          run.distance_sum += hit.distance;
        if (hit.found) {
          run.selection_checksum +=
              hit.distance + 1.0e-7 * hit.point.x +
              3.0e-8 * hit.point.y + 7.0e-8 * hit.point.z +
              1.0e-11 * static_cast<double>(hit.patch_index) +
              1.0e-13 * static_cast<double>(
                              static_cast<int>(hit.kind));
        }
        if (selections != nullptr) {
          selections->push_back(
              {hit.point, hit.distance, hit.patch_index, hit.kind,
               hit.found, hit.bounded_closest_certified});
        }
      }
      run.elapsed_us =
          std::chrono::duration<double, std::micro>(
              std::chrono::steady_clock::now() - begin)
              .count();
      return run;
    };
    const auto same_global_counters =
        [](const GlobalRun &a, const GlobalRun &b) {
          return a.hits == b.hits &&
                 a.certified == b.certified &&
                 a.owner_hint_hits == b.owner_hint_hits &&
                 a.owner_hint_candidates ==
                     b.owner_hint_candidates &&
                 a.spatial_fallbacks == b.spatial_fallbacks &&
                 a.axis_fallbacks == b.axis_fallbacks &&
                 a.candidates == b.candidates &&
                 a.candidates_by_kind == b.candidates_by_kind &&
                 a.empty_short_circuits_by_kind ==
                     b.empty_short_circuits_by_kind &&
                 a.pruned == b.pruned &&
                 a.nodes_visited == b.nodes_visited &&
                 a.distance_sum == b.distance_sum &&
                 a.selection_checksum == b.selection_checksum;
        };
    const auto same_selections =
        [](std::span<const SelectedProjection> a,
           std::span<const SelectedProjection> b) {
          if (a.size() != b.size())
            return false;
          for (std::size_t i = 0; i < a.size(); ++i) {
            if (a[i].found != b[i].found ||
                a[i].certified != b[i].certified ||
                a[i].kind != b[i].kind ||
                a[i].patch_index != b[i].patch_index ||
                a[i].distance != b[i].distance ||
                a[i].point.x != b[i].point.x ||
                a[i].point.y != b[i].point.y ||
                a[i].point.z != b[i].point.z)
              return false;
          }
          return true;
        };

    if (ab_only) {
      constexpr std::size_t kMaximumAbQueries = 512U;
      const std::size_t ab_count =
          std::min(kMaximumAbQueries, queries.size());
      std::vector<Query> ab_queries;
      ab_queries.reserve(ab_count);
      for (std::size_t sample = 0; sample < ab_count; ++sample)
        ab_queries.push_back(
            queries[(sample * queries.size()) / ab_count]);

      struct PolicyResult {
        const char *name;
        GlobalRun timed;
        std::vector<SelectedProjection> selections;
      };
      const auto run_policy =
          [&](const char *name, ProjectionHintPolicy hint_policy) {
            const GlobalRun warm =
                run_global(complex, index, ab_queries, nullptr,
                           hint_policy, false);
            PolicyResult result{name, {}, {}};
            result.timed = run_global(
                complex, index, ab_queries, &result.selections,
                hint_policy, false);
            if (!same_global_counters(warm, result.timed))
              throw std::runtime_error(
                  std::string("A/B counters changed after warm-up: ") +
                  name);
            return result;
          };
      const PolicyResult current =
          run_policy("owner_plus_neighbors_audit",
                     ProjectionHintPolicy::owner_and_neighbors);
      const PolicyResult owner_only =
          run_policy("owner_only_production",
                     ProjectionHintPolicy::owner_only);
      const PolicyResult bvh_only =
          run_policy("bvh_only",
                     ProjectionHintPolicy::bvh_only);
      const auto report_first_mismatch =
          [&](const char *policy,
              std::span<const SelectedProjection> reference,
              std::span<const SelectedProjection> candidate) {
            if (same_selections(reference, candidate))
              return true;
            for (std::size_t i = 0;
                 i < reference.size() && i < candidate.size(); ++i) {
              const SelectedProjection &a = reference[i];
              const SelectedProjection &b = candidate[i];
              if (a.found == b.found &&
                  a.certified == b.certified &&
                  a.kind == b.kind &&
                  a.patch_index == b.patch_index &&
                  a.distance == b.distance &&
                  a.point.x == b.point.x &&
                  a.point.y == b.point.y &&
                  a.point.z == b.point.z)
                continue;
              const Query &query = ab_queries[i];
              std::cerr
                  << std::setprecision(17)
                  << "ab_first_mismatch policy=" << policy
                  << " sample=" << i << " owner=" << query.atom
                  << " query=(" << query.point.x << ','
                  << query.point.y << ',' << query.point.z << ')'
                  << " reference=(found=" << a.found
                  << ",cert=" << a.certified
                  << ",kind=" << static_cast<int>(a.kind)
                  << ",patch=" << a.patch_index
                  << ",distance=" << a.distance << ",point="
                  << a.point.x << ':' << a.point.y << ':'
                  << a.point.z << ") candidate=(found="
                  << b.found << ",cert=" << b.certified
                  << ",kind=" << static_cast<int>(b.kind)
                  << ",patch=" << b.patch_index
                  << ",distance=" << b.distance << ",point="
                  << b.point.x << ':' << b.point.y << ':'
                  << b.point.z << ")\n";
              break;
            }
            return false;
          };
      const bool owner_equal = report_first_mismatch(
          "owner_only", current.selections, owner_only.selections);
      const bool bvh_equal = report_first_mismatch(
          "bvh_only", current.selections, bvh_only.selections);
      if (!owner_equal || !bvh_equal)
        throw std::runtime_error(
            "A/B hint policies changed selected projection");
      std::vector<SelectedProjection> profiled_selections;
      const GlobalRun kind_profile =
          run_global(complex, index, ab_queries, &profiled_selections,
                     ProjectionHintPolicy::owner_only, true);
      if (!same_global_counters(owner_only.timed, kind_profile) ||
          !same_selections(owner_only.selections,
                           profiled_selections))
        throw std::runtime_error(
            "kind profiling changed owner-only projection results");

      const double ab_query_count =
          static_cast<double>(ab_queries.size());
      std::cout
          << std::setprecision(12)
          << "ab_policy,query_distribution,queries,warmups,"
             "timed_repeats,patch_build_ms,index_build_ms,"
             "active_pair_cache_build_ms,active_pair_cache_bytes,"
             "empty_contact_domains,empty_reentrant_domains,"
             "contact_cache_vertices_total,contact_cache_p50,"
             "contact_cache_p95,contact_cache_p99,contact_cache_max,"
             "reentrant_cache_vertices_total,reentrant_cache_p50,"
             "reentrant_cache_p95,reentrant_cache_p99,"
             "reentrant_cache_max,us_per_query,hits,certified,"
             "candidates_per_query,owner_hint_candidates_per_query,"
             "spatial_nodes_visited_per_query,"
             "owner_hint_hit_fraction,spatial_fallback_fraction,"
             "axis_fallback_fraction,checksum,"
             "selection_equal_reference\n";
      const auto print_policy =
          [&](const PolicyResult &policy, bool equal_reference) {
            const GlobalRun &run = policy.timed;
            std::cout
                << policy.name
                << ",deterministic_stratified_512,"
                << ab_queries.size() << ",1,1,"
                << patch_build_ms << ',' << index_build_ms << ','
                << index.active_pair_cache_build_ms << ','
                << index.active_pair_cache_bytes() << ','
                << empty_contact_domains << ','
                << empty_reentrant_domains << ','
                << contact_cache.total << ',' << contact_cache.p50
                << ',' << contact_cache.p95 << ','
                << contact_cache.p99 << ',' << contact_cache.maximum
                << ',' << reentrant_cache.total << ','
                << reentrant_cache.p50 << ','
                << reentrant_cache.p95 << ','
                << reentrant_cache.p99 << ','
                << reentrant_cache.maximum << ','
                << run.elapsed_us / ab_query_count << ','
                << run.hits << ',' << run.certified << ','
                << static_cast<double>(run.candidates) /
                       ab_query_count
                << ','
                << static_cast<double>(
                       run.owner_hint_candidates) /
                       ab_query_count
                << ','
                << static_cast<double>(run.nodes_visited) /
                       ab_query_count
                << ','
                << static_cast<double>(run.owner_hint_hits) /
                       ab_query_count
                << ','
                << static_cast<double>(run.spatial_fallbacks) /
                       ab_query_count
                << ','
                << static_cast<double>(run.axis_fallbacks) /
                       ab_query_count
                << ',' << run.selection_checksum << ','
                << equal_reference << '\n';
          };
      print_policy(current, true);
      print_policy(
          owner_only,
          same_selections(current.selections,
                          owner_only.selections));
      print_policy(
          bvh_only,
          same_selections(current.selections, bvh_only.selections));
      const auto projector_us_per_candidate =
          [&](std::size_t kind) {
            if (kind_profile.candidates_by_kind[kind] == 0U)
              return 0.0;
            return static_cast<double>(
                       kind_profile.projector_nanoseconds_by_kind[kind]) /
                   1000.0 /
                   static_cast<double>(
                       kind_profile.candidates_by_kind[kind]);
          };
      std::cout
          << "kind_profile,policy,query_distribution,queries,"
             "profiled_wall_us_per_query,"
             "contact_candidates_total,contact_candidates_per_query,"
             "contact_empty_short_circuits,contact_projector_us_total,"
             "contact_projector_us_per_candidate,"
             "toroidal_candidates_total,toroidal_candidates_per_query,"
             "toroidal_empty_short_circuits,toroidal_projector_us_total,"
             "toroidal_projector_us_per_candidate,"
             "reentrant_candidates_total,reentrant_candidates_per_query,"
             "reentrant_empty_short_circuits,"
             "reentrant_projector_us_total,"
             "reentrant_projector_us_per_candidate\n"
          << "kind_profile,owner_only_production,"
             "deterministic_stratified_512,"
          << ab_queries.size() << ','
          << kind_profile.elapsed_us / ab_query_count;
      for (std::size_t kind = 0; kind < 3U; ++kind) {
        std::cout
            << ',' << kind_profile.candidates_by_kind[kind] << ','
            << static_cast<double>(
                   kind_profile.candidates_by_kind[kind]) /
                   ab_query_count
            << ',' << kind_profile.empty_short_circuits_by_kind[kind]
            << ','
            << static_cast<double>(
                   kind_profile.projector_nanoseconds_by_kind[kind]) /
                   1000.0
            << ',' << projector_us_per_candidate(kind);
      }
      std::cout << '\n';
      return EXIT_SUCCESS;
    }

    // Warm both projector paths once so allocation and instruction-cache
    // effects are outside the three reported timed repeats.
    const ContactRun contact_warm = run_contact();
    const GlobalRun global_warm =
        run_global(complex, index, queries, nullptr,
                   ProjectionHintPolicy::owner_only, false);
    constexpr std::size_t kTimedRepeats = 3U;
    std::array<ContactRun, kTimedRepeats> contact_runs;
    std::array<GlobalRun, kTimedRepeats> global_runs;
    for (std::size_t repeat = 0; repeat < kTimedRepeats; ++repeat) {
      contact_runs[repeat] = run_contact();
      global_runs[repeat] =
          run_global(complex, index, queries, nullptr,
                     ProjectionHintPolicy::owner_only, false);
    }
    const auto same_contact_counters =
        [](const ContactRun &a, const ContactRun &b) {
          return a.hits == b.hits &&
                 a.certified == b.certified &&
                 a.distance_sum == b.distance_sum;
        };
    if (!same_contact_counters(contact_warm, contact_runs[0]) ||
        !same_global_counters(global_warm, global_runs[0]))
      throw std::runtime_error(
          "projection counters changed after warm-up");
    for (std::size_t repeat = 1U; repeat < kTimedRepeats; ++repeat) {
      if (!same_contact_counters(contact_runs[0], contact_runs[repeat]) ||
          !same_global_counters(global_runs[0], global_runs[repeat]))
        throw std::runtime_error(
            "projection counters changed across timed repeats");
    }
    const auto median_and_best_us = [](const auto &runs) {
      std::array<double, kTimedRepeats> timings{};
      for (std::size_t repeat = 0; repeat < kTimedRepeats; ++repeat)
        timings[repeat] = runs[repeat].elapsed_us;
      std::sort(timings.begin(), timings.end());
      return std::pair{timings[kTimedRepeats / 2U], timings.front()};
    };
    const auto [contact_median_us, contact_best_us] =
        median_and_best_us(contact_runs);
    const auto [global_median_us, global_best_us] =
        median_and_best_us(global_runs);
    const ContactRun &contact = contact_runs[0];
    const GlobalRun &global = global_runs[0];

    const double query_count = static_cast<double>(queries.size());
    const double contact_median_us_per_query =
        contact_median_us / query_count;
    const double contact_best_us_per_query =
        contact_best_us / query_count;
    const double global_median_us_per_query =
        global_median_us / query_count;
    const double global_best_us_per_query =
        global_best_us / query_count;

    std::cout << std::setprecision(12)
              << "atoms,contacts,tori,reentrant,index_cert,bvh_cert,"
                 "contact_constraints_mean,contact_constraints_p50,"
                 "contact_constraints_p95,contact_constraints_p99,"
                 "contact_constraints_max,"
                 "patch_build_ms,index_build_ms,"
                 "active_pair_cache_build_ms,active_pair_cache_bytes,"
                 "empty_contact_domains,empty_reentrant_domains,"
                 "contact_cache_vertices_total,contact_cache_p50,"
                 "contact_cache_p95,contact_cache_p99,"
                 "contact_cache_max,reentrant_cache_vertices_total,"
                 "reentrant_cache_p50,reentrant_cache_p95,"
                 "reentrant_cache_p99,reentrant_cache_max,"
                 "torus_reentrant_audit_pairs,torus_reentrant_audit_ms,"
                 "query_distribution,global_hint_policy,queries,"
                 "timed_repeats,"
                 "contact_us_per_query_median,"
                 "contact_us_per_query_best,contact_hits,"
                 "contact_certified,global_us_per_query_median,"
                 "global_us_per_query_best,global_hits,"
                 "global_certified,candidates_per_query,"
                 "bound_pruned_raw_per_query,"
                 "owner_hint_candidates_per_query,"
                 "spatial_nodes_visited_per_query,"
                 "owner_hint_hit_fraction,spatial_fallback_fraction,"
                 "axis_fallback_fraction,index_bytes,"
                 "stamp_required_bytes_per_worker,"
                 "stamp_actual_bytes_per_worker,checksum\n"
              << atoms.size() << ',' << complex.contacts.size() << ','
              << complex.tori.size() << ',' << complex.reentrant.size()
              << ',' << index.bounded_domains_certified << ','
              << index.bvh_invariants_certified << ','
              << mean_constraints << ',' << constraints_p50 << ','
              << constraints_p95 << ',' << constraints_p99 << ','
              << maximum_constraints << ','
              << patch_build_ms << ',' << index_build_ms << ','
              << index.active_pair_cache_build_ms << ','
              << index.active_pair_cache_bytes() << ','
              << empty_contact_domains << ','
              << empty_reentrant_domains << ','
              << contact_cache.total << ',' << contact_cache.p50 << ','
              << contact_cache.p95 << ',' << contact_cache.p99 << ','
              << contact_cache.maximum << ','
              << reentrant_cache.total << ',' << reentrant_cache.p50
              << ',' << reentrant_cache.p95 << ','
              << reentrant_cache.p99 << ','
              << reentrant_cache.maximum << ','
              << index.torus_reentrant_pairs_audited << ','
              << index.torus_reentrant_audit_ms << ','
              << "four_synthetic_contact_directions_per_atom,"
              << "owner_only_production,"
              << queries.size() << ',' << kTimedRepeats << ','
              << contact_median_us_per_query << ','
              << contact_best_us_per_query << ','
              << contact.hits << ',' << contact.certified << ','
              << global_median_us_per_query << ','
              << global_best_us_per_query << ',' << global.hits << ','
              << global.certified << ','
              << static_cast<double>(global.candidates) / query_count
              << ','
              << static_cast<double>(global.pruned) / query_count << ','
              << static_cast<double>(global.owner_hint_candidates) /
                     query_count
              << ','
              << static_cast<double>(global.nodes_visited) / query_count
              << ','
              << static_cast<double>(global.owner_hint_hits) /
                     query_count
              << ','
              << static_cast<double>(global.spatial_fallbacks) /
                     query_count
              << ','
              << static_cast<double>(global.axis_fallbacks) /
                     query_count
              << ',' << index.retained_bytes() << ','
              << index.stamp_bytes_per_worker() << ','
              << global.actual_stamp_bytes << ','
              << contact.distance_sum + global.distance_sum << '\n';
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "patch_projection_benchmark=fail reason="
              << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
