param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_tile_kernel_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_tile_kernel_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_tile_kernel_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI leaf-tile root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $Compiler)) {
    throw "Missing clang++ compiler: $Compiler"
}

New-CleanDir $WorkRoot

$source = @'
#include "titan_executor.h"
#include "tabi/tabi_types.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

double rms(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values)
    sum += value * value;
  return values.empty() ? 0.0
                        : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  require_true("same vector size", reference.size() == candidate.size());
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

std::vector<TitanPBE::Tabi::Panel> build_panels() {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  opt.surface_point_count = 24;
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Atom> atoms;
  for (int i = 0; i < 10; ++i) {
    const double x = -6.5 + 1.45 * static_cast<double>(i);
    const double y = (i % 2 == 0) ? -0.28 : 0.31;
    const double z = (i % 3 == 0) ? 0.25 : -0.21;
    const double q = (i % 2 == 0) ? 0.35 : -0.30;
    atoms.push_back({x, y, z, q, 1.50 + 0.04 * static_cast<double>(i % 4)});
  }
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, opt, &result);
  require_true("triangle panels", panels.size() > 150);
  return panels;
}

} // namespace

int main() {
  std::vector<TitanPBE::Tabi::Panel> panels = build_panels();
  const std::size_t n = panels.size();
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.06 * std::sin(0.013 * static_cast<double>(i));
    unknowns[n + i] = 0.03 * std::cos(0.019 * static_cast<double>(i));
  }

  TitanExecutor executor1(1);
  TitanPBE::Tabi::TreeWorkspace scalar_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor1, 32, 0.8,
                                           "exact", "scalar");
  TitanPBE::Tabi::TreeWorkspace tile_scalar_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor1, 32, 0.8,
                                           "exact", "tile-scalar");
  TitanPBE::Tabi::TreeWorkspace tile_avx2_tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor1, 32, 0.8,
                                           "exact", "tile-avx2");
  require_true("scalar kernel", scalar_tree.leaf_panel_kernel == "scalar");
  require_true("tile scalar kernel",
               tile_scalar_tree.leaf_panel_kernel == "tile-scalar");
  require_true("tile avx2 kernel",
               tile_avx2_tree.leaf_panel_kernel == "tile-avx2");
  require_true("tile nodes emitted",
               tile_scalar_tree.interactions.leaf_tile_node_count > 0);
  require_true("tile panel count emitted",
               tile_scalar_tree.interactions.leaf_tile_panel_count > 0);
  require_true("tile leaves reduce flat list",
               tile_scalar_tree.interactions.leaf_indices.size() <
                   scalar_tree.interactions.leaf_indices.size());
  require_true("tile total preserves leaf panels",
               tile_scalar_tree.interactions.leaf_pair_count ==
                   tile_scalar_tree.interactions.leaf_indices.size() +
                       tile_scalar_tree.interactions.leaf_tile_panel_count);

  std::vector<double> scalar_out;
  std::vector<double> tile_scalar_out;
  std::vector<double> tile_avx2_out;
  TitanPBE::Tabi::TreecodeStats scalar_stats;
  TitanPBE::Tabi::TreecodeStats tile_scalar_stats;
  TitanPBE::Tabi::TreecodeStats tile_avx2_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, scalar_tree,
      "centroid", "off", executor1, scalar_out, &scalar_stats);
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true,
      tile_scalar_tree, "centroid", "off", executor1, tile_scalar_out,
      &tile_scalar_stats);
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, tile_avx2_tree,
      "centroid", "off", executor1, tile_avx2_out, &tile_avx2_stats);
  require_true("tile scalar stats kernel",
               tile_scalar_stats.leaf_panel_kernel == "tile-scalar");
  require_true("tile avx2 stats kernel",
               tile_avx2_stats.leaf_panel_kernel == "tile-avx2");
  require_true("tile stats active",
               tile_scalar_stats.leaf_tile_status == "active");
  require_true("tile stats node count",
               tile_scalar_stats.leaf_tile_interactions > 0);
  require_true("tile stats panel count",
               tile_scalar_stats.leaf_tile_panel_interactions > 0);
  require_true("tile scalar vs flat scalar exact",
               relative_rms(scalar_out, tile_scalar_out) <= 1.0e-12);
  require_true("tile avx2 vs tile scalar exact",
               relative_rms(tile_scalar_out, tile_avx2_out) <= 1.0e-12);
  require_true("tile avx2 counts same",
               tile_avx2_stats.leaf_tile_interactions ==
                   tile_scalar_stats.leaf_tile_interactions);

  TitanExecutor executor4(4);
  TitanPBE::Tabi::TreeWorkspace tile_tree4 =
      TitanPBE::Tabi::build_tree_workspace(panels, executor4, 32, 0.8,
                                           "exact", "tile-avx2");
  std::vector<double> tile_out4;
  TitanPBE::Tabi::TreecodeStats tile_stats4;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, tile_tree4,
      "centroid", "off", executor4, tile_out4, &tile_stats4);
  require_true("tile nt deterministic",
               relative_rms(tile_avx2_out, tile_out4) < 1.0e-14);

  std::cout << "PBE TABI leaf-tile kernel smoke PASS panels=" << n
            << " flat_residual="
            << tile_scalar_stats.leaf_residual_panel_interactions
            << " tile_nodes=" << tile_scalar_stats.leaf_tile_interactions
            << " tile_panels="
            << tile_scalar_stats.leaf_tile_panel_interactions
            << "\n";
  return 0;
}
'@

Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII

$compileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',
    ("-I" + (Join-Path $RepoRoot "src")),
    (Join-Path $RepoRoot 'src\tabi\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_aux_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_calderon_trace_space.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_analytic_ses.cpp'),
    (Get-ChildItem -LiteralPath (Join-Path $RepoRoot 'src\tabi') -Filter 'tabi_ses_*.cpp' |
        Sort-Object Name | ForEach-Object { $_.FullName }),
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI leaf-tile kernel fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-tile kernel fixture failed."
}
