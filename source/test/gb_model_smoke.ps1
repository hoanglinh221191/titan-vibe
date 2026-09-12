param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\gb_model_smoke")
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$src = Join-Path $repo "test\titan_gb_model_smoke.cpp"
$exe = Join-Path $WorkRoot "titan_gb_model_smoke.exe"
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

& $Compiler -std=c++23 -O2 -mavx2 -mfma "-I$(Join-Path $repo 'src')" $src -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "gb_model_smoke: failed to compile"
}
& $exe
if ($LASTEXITCODE -ne 0) {
  throw "gb_model_smoke: executable failed"
}
Write-Host "PASS gb_model_smoke"
