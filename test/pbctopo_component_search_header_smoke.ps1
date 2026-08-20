param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbctopo_component_search_header_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'pbctopo_component_search_header_smoke.cpp'
$obj = Join-Path $WorkRoot 'pbctopo_component_search_header_smoke.o'

@'
#include "pbctopo/component_search.h"

#include <cstddef>
#include <type_traits>

int main() {
  using titan_pbctopo::PbctopoComponentSearchCallbacks;
  static_assert(std::is_default_constructible_v<PbctopoComponentSearchCallbacks>);
  if (titan_pbctopo::pbctopo_component_assignment_upper_bound(4, 27) != 19683)
    return 2;
  if (titan_pbctopo::pbctopo_component_assignment_upper_bound(5, 27, 100001) !=
      100001)
    return 3;
  PbctopoComponentSearchCallbacks callbacks;
  return callbacks.build_seed == nullptr && callbacks.score_assignment == nullptr &&
                 callbacks.better_candidate == nullptr
             ? 0
             : 1;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc -c $src -o $obj
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile pbctopo component search header smoke'
}

Write-Output 'PASS pbctopo_component_search_header_smoke'
