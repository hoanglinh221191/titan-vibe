param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\cscd_module_smoke")
)

$ErrorActionPreference = "Stop"
$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$srcDir = Join-Path $repo "src"
$exe = Join-Path $WorkRoot "cscd_module_smoke.exe"
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

& $Compiler -std=c++23 -O2 "-I$srcDir" `
  (Join-Path $repo "test\cscd_module_smoke.cpp") `
  (Join-Path $srcDir "impl_cscd_class.cpp") `
  (Join-Path $srcDir "titan_executor.cpp") -o $exe
if ($LASTEXITCODE -ne 0) {
  throw "cscd_module_smoke: failed to compile"
}
& $exe (Join-Path $WorkRoot "out")
if ($LASTEXITCODE -ne 0) {
  throw "cscd_module_smoke: executable failed"
}
Write-Host "PASS cscd_module_smoke"
