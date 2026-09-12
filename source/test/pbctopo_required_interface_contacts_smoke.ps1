param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_required_interface_contacts_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
if (Test-Path -LiteralPath $WorkRoot) {
    Remove-Item -LiteralPath $WorkRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$contacts = Join-Path $WorkRoot 'required_contacts.dat'
$invalidContacts = Join-Path $WorkRoot 'invalid_contacts.dat'
$src = Join-Path $WorkRoot 'smoke.cpp'
$exe = Join-Path $WorkRoot 'smoke.exe'

@'
# selected atom indices are 1-based
1 2 4.0 2.0
2 1 5.0 1.0 ; duplicate: stricter distance and larger weight win
'@ | Set-Content -LiteralPath $contacts -Encoding ASCII
'1 3 4.0' | Set-Content -LiteralPath $invalidContacts -Encoding ASCII

@'
#include "pbctopo/certification.h"
#include "pbctopo/required_contacts.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace titan_pbctopo;

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL pbctopo_required_interface_contacts_smoke: "
              << message << "\n";
    std::exit(1);
  }
}

int main(int argc, char **argv) {
  expect(argc == 3, "fixture paths were not supplied");
  const auto loaded = load_pbctopo_required_contacts(argv[1], 2);
  expect(loaded.has_value() && loaded->size() == 1,
         "required-contact parser did not deduplicate an unordered pair");
  expect((*loaded)[0].a == 0 && (*loaded)[0].b == 1 &&
             (*loaded)[0].max_distance_A == 4.0 &&
             (*loaded)[0].weight == 2.0,
         "required-contact parser did not preserve strict deterministic data");
  expect(!load_pbctopo_required_contacts(argv[2], 2).has_value(),
         "out-of-range required-contact atom index was accepted");

  std::vector<PbctopoAtom> atoms(2);
  atoms[1].x = 3.0;
  const std::vector<std::vector<std::size_t>> components{{0}, {1}};
  const box Box(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
  const auto constraint_context =
      build_pbctopo_required_contact_constraint_context(
          atoms, components, Box, *loaded);
  expect(constraint_context.valid && !constraint_context.infeasible &&
             constraint_context.groups.size() == 1,
         "required-contact component constraint context is invalid");
  expect(pbctopo_required_contact_relative_shift_allowed(
             constraint_context, atoms, 0, 1, {0, 0, 0}),
         "intact relative shift was removed from the component domain");
  expect(!pbctopo_required_contact_relative_shift_allowed(
             constraint_context, atoms, 0, 1, {1, 0, 0}),
         "split relative shift survived the component-domain constraint");
  expect(!pbctopo_required_contact_relative_shift_allowed(
             constraint_context, atoms, 1, 0,
             {std::numeric_limits<std::int64_t>::min(), 0, 0}),
         "required-contact reverse shift did not fail closed on overflow");

  std::vector<PbctopoAtom> near_boundary_atoms(2);
  near_boundary_atoms[1].x = std::sqrt(16.0 + 0.5e-6);
  const auto near_boundary_context =
      build_pbctopo_required_contact_constraint_context(
          near_boundary_atoms, components, Box, *loaded);
  expect(pbctopo_required_contact_relative_shift_allowed(
             near_boundary_context, near_boundary_atoms, 0, 1, {0, 0, 0}),
         "search rejected a required contact inside the shared tolerance");
  const std::vector<Vec3> near_boundary_points{
      {0.0, 0.0, 0.0}, {near_boundary_atoms[1].x, 0.0, 0.0}};
  const auto near_boundary_residuals = evaluate_pbctopo_required_contacts(
      near_boundary_atoms, near_boundary_points, {}, *loaded);
  expect(near_boundary_residuals.input_valid &&
             near_boundary_residuals.explicit_contacts_lost == 0,
         "certificate rejected a required contact accepted by search");

  std::vector<PbctopoAtom> outside_boundary_atoms(2);
  outside_boundary_atoms[1].x = std::sqrt(16.0 + 2.0e-6);
  const auto outside_boundary_context =
      build_pbctopo_required_contact_constraint_context(
          outside_boundary_atoms, components, Box, *loaded);
  expect(!pbctopo_required_contact_relative_shift_allowed(
             outside_boundary_context, outside_boundary_atoms, 0, 1,
             {0, 0, 0}),
         "search accepted a required contact outside the shared tolerance");
  const std::vector<Vec3> outside_boundary_points{
      {0.0, 0.0, 0.0}, {outside_boundary_atoms[1].x, 0.0, 0.0}};
  const auto outside_boundary_residuals = evaluate_pbctopo_required_contacts(
      outside_boundary_atoms, outside_boundary_points, {}, *loaded);
  expect(outside_boundary_residuals.input_valid &&
             outside_boundary_residuals.explicit_contacts_lost == 1,
         "certificate accepted a required contact rejected by search");

  std::vector<Vec3> intact{{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
  std::vector<Vec3> split{{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}};
  const auto intact_residuals = evaluate_pbctopo_required_contacts(
      atoms, intact, {}, *loaded);
  expect(intact_residuals.input_valid &&
             intact_residuals.explicit_contacts_checked == 1 &&
             intact_residuals.explicit_contacts_lost == 0 &&
             intact_residuals.explicit_weight_checked == 2.0,
         "intact explicit interface contact did not pass");
  const auto split_residuals = evaluate_pbctopo_required_contacts(
      atoms, split, {}, *loaded);
  expect(split_residuals.input_valid &&
             split_residuals.explicit_contacts_lost == 1 &&
             split_residuals.explicit_weight_lost == 2.0,
         "split explicit interface contact was not rejected");

  PbctopoCutCandidate candidate;
  candidate.valid = true;
  candidate.assignment.valid = true;
  PbctopoImageEquationResiduals image_equations;
  image_equations.input_valid = true;
  PbctopoAssignmentPointResiduals assignment_points;
  assignment_points.input_valid = true;
  PbctopoStericAudit sterics;
  sterics.input_valid = true;
  sterics.scan_complete = true;
  expect(certify_pbctopo_candidate(candidate, image_equations,
                                    assignment_points, sterics,
                                    intact_residuals),
         "intact explicit interface certificate failed");
  candidate.valid = true;
  candidate.assignment.valid = true;
  sterics.scan_complete = false;
  expect(!certify_pbctopo_candidate(candidate, image_equations,
                                     assignment_points, sterics,
                                     intact_residuals),
         "incomplete steric scan was certified");
  candidate.valid = true;
  candidate.assignment.valid = true;
  sterics.scan_complete = true;
  expect(!certify_pbctopo_candidate(candidate, image_equations,
                                     assignment_points, sterics,
                                     split_residuals),
         "split explicit interface was certified");

  std::cout << "PASS pbctopo_required_interface_contacts_smoke\n";
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

$compilerPath = (Get-Command $Compiler -ErrorAction Stop).Source
Push-Location $repo
try {
    & $compilerPath -std=c++23 -O2 -Isrc $src `
        src\pbctopo\required_contacts.cpp `
        src\pbctopo\certification.cpp `
        src\pbctopo\steric_clash.cpp `
        src\pbctopo\lattice_math.cpp -o $exe
    if ($LASTEXITCODE -ne 0) {
        throw 'Failed to compile PBCTOPO required-interface contact smoke'
    }
    & $exe $contacts $invalidContacts
    if ($LASTEXITCODE -ne 0) {
        throw "PBCTOPO required-interface contact smoke failed: $LASTEXITCODE"
    }
}
finally {
    Pop-Location
}
