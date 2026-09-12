param()

$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot

function Assert-PathExists {
  param(
    [Parameter(Mandatory = $true)][string]$Path,
    [Parameter(Mandatory = $true)][string]$Description
  )
  if (-not (Test-Path -LiteralPath $Path)) {
    throw "Missing ${Description}: ${Path}"
  }
}

function Assert-ContainsLiteral {
  param(
    [Parameter(Mandatory = $true)][string]$Path,
    [Parameter(Mandatory = $true)][string]$Needle
  )
  $text = Get-Content -LiteralPath $Path -Raw
  if (-not $text.Contains($Needle)) {
    throw "Expected ${Path} to contain literal: ${Needle}"
  }
}

function Assert-NotMatches {
  param(
    [Parameter(Mandatory = $true)][string]$Path,
    [Parameter(Mandatory = $true)][string]$Pattern,
    [Parameter(Mandatory = $true)][string]$Description
  )
  $text = Get-Content -LiteralPath $Path -Raw
  if ($text -match $Pattern) {
    throw "${Description}: ${Path}"
  }
}

$gridDir = Join-Path $RepoRoot 'src\pbe\grid'
Assert-PathExists -Path $gridDir -Description 'PBE grid solver directory'
Assert-PathExists -Path (Join-Path $gridDir 'grid_solver.h') -Description 'PBE grid solver API header'
Assert-PathExists -Path (Join-Path $gridDir 'grid_workspace.h') -Description 'PBE grid workspace header'
Assert-PathExists -Path (Join-Path $gridDir 'ses_mask.h') -Description 'PBE grid SES mask header'
Assert-PathExists -Path (Join-Path $gridDir 'ses_mask.cpp') -Description 'PBE grid SES mask implementation'

$gridFiles = @(
  'grid_surface.cpp',
  'charge_operator.cpp',
  'linear_solver.cpp',
  'multigrid.cpp',
  'focus_probe.cpp',
  'gaussian_rpb.cpp'
)

foreach ($file in $gridFiles) {
  Assert-PathExists -Path (Join-Path $gridDir $file) -Description "grid solver file ${file}"
  $oldPath = Join-Path (Join-Path $RepoRoot 'src\tabi') $file
  if (Test-Path -LiteralPath $oldPath) {
    throw "Grid solver file still lives under src\tabi: ${oldPath}"
  }
}

$gridCommon = Join-Path $gridDir 'grid_common.h'
Assert-PathExists -Path $gridCommon -Description 'PBE grid common helper header'

$analyzerImpl = Join-Path $RepoRoot 'src\tabi\analyzer_impl.h'
Assert-ContainsLiteral -Path $analyzerImpl -Needle '#include "pbe/grid/grid_solver.h"'
Assert-ContainsLiteral -Path $analyzerImpl -Needle '#include "pbe/grid/grid_workspace.h"'
Assert-ContainsLiteral -Path $analyzerImpl -Needle 'GridWorkspace workspace;'
Assert-ContainsLiteral -Path $analyzerImpl -Needle 'GridWorkspace grid_;'
Assert-ContainsLiteral -Path $analyzerImpl -Needle 'FieldSolve solve_grid_field'
Assert-NotMatches -Path $analyzerImpl -Pattern 'AlignedDoubleVector\s+(?:phi_|rho_|eps_x_|eps_y_|eps_z_|denom_|screen_)\s*(?:[;=])' -Description 'AnalyzerImpl still declares raw grid double buffer member'
Assert-NotMatches -Path $analyzerImpl -Pattern 'std::vector<unsigned char>\s+(?:inside_|ion_accessible_)\s*(?:[;=])' -Description 'AnalyzerImpl still declares raw grid mask member'

$oldRawTokenPattern = '\b(phi_|rho_|eps_x_|eps_y_|eps_z_|denom_|screen_|inside_|ion_accessible_)\b'
$sourceFilesToCheck = @(
  (Join-Path $RepoRoot 'src\tabi\analyzer_impl.cpp'),
  (Join-Path $RepoRoot 'src\tabi\audit_io.cpp')
)
foreach ($file in $gridFiles) {
  $sourceFilesToCheck += Join-Path $gridDir $file
}
foreach ($path in $sourceFilesToCheck) {
  Assert-NotMatches -Path $path -Pattern $oldRawTokenPattern -Description 'Grid source still references direct AnalyzerImpl raw grid token'
}

foreach ($file in $gridFiles) {
  $path = Join-Path $gridDir $file
  Assert-ContainsLiteral -Path $path -Needle '#include "pbe/grid/grid_common.h"'
  $text = Get-Content -LiteralPath $path -Raw
  if ($text -match 'std::size_t\s+idx3\s*\(') {
    throw "Grid source still defines local idx3 helper: ${path}"
  }
  if ($text -match 'std::size_t\s+effective_workers\s*\(') {
    throw "Grid source still defines local effective_workers helper: ${path}"
  }
}

Assert-PathExists -Path (Join-Path $RepoRoot 'src\tabi\tabi_solver.cpp') -Description 'TABI solver remains in src\tabi'
Assert-PathExists -Path (Join-Path $RepoRoot 'src\tabi\pbsa.cpp') -Description 'PBSA remains outside grid solver refactor'

$makefile = Join-Path $RepoRoot 'Makefile'
Assert-ContainsLiteral -Path $makefile -Needle 'PBE_GRID_SRCS'
Assert-ContainsLiteral -Path $makefile -Needle '$(SRC_ROOT)/pbe/grid/*.cpp'
Assert-ContainsLiteral -Path $makefile -Needle 'PBE_SRCS  = $(wildcard $(SRC_ROOT)/tabi/*.cpp) $(PBE_GRID_SRCS)'

$gridSurface = Join-Path $gridDir 'grid_surface.cpp'
Assert-ContainsLiteral -Path $gridSurface -Needle 'point_near_inflated_atom_surface'
$gridSurfaceText = Get-Content -LiteralPath $gridSurface -Raw
$hysteresisStart = $gridSurfaceText.IndexOf('void AnalyzerImpl::apply_stability_hysteresis')
$hysteresisEnd = $gridSurfaceText.IndexOf('void AnalyzerImpl::store_stability_medium_state')
if ($hysteresisStart -lt 0 -or $hysteresisEnd -le $hysteresisStart) {
  throw 'Unable to locate apply_stability_hysteresis block in grid_surface.cpp.'
}
$hysteresisText = $gridSurfaceText.Substring($hysteresisStart, $hysteresisEnd - $hysteresisStart)
if ($hysteresisText.Contains('for (const Atom &atom : surface_atoms)')) {
  throw 'Grid hysteresis still scans all surface atoms per changed voxel.'
}

$multigrid = Join-Path $gridDir 'multigrid.cpp'
$multigridText = Get-Content -LiteralPath $multigrid -Raw
if ($multigridText.Contains('const AlignedDoubleVector before = levels.front().phi;')) {
  throw 'Grid multigrid still allocates a fresh fine-grid copy every V-cycle.'
}
Assert-ContainsLiteral -Path $multigrid -Needle 'previous_phi'

Write-Host 'PBE grid layout smoke passed.'
