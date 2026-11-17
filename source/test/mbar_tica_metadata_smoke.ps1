param(
  [string]$Compiler = "C:\msys64\ucrt64\bin\clang++.exe",
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\mbar_tica_metadata_smoke")
)

$ErrorActionPreference = "Stop"

$repo = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

if (-not (Test-Path -LiteralPath $Compiler)) {
  throw "Missing clang++ compiler: $Compiler"
}

$parser = (Get-ChildItem -LiteralPath (Join-Path $repo "src\cli") -File -Filter "parse_*.cpp" |
  Sort-Object Name | ForEach-Object { Get-Content -Raw -LiteralPath $_.FullName }) -join "`n"
$mbarEnableBlock = [regex]::Match(
  $parser,
  'if \(string\(argv\[i\]\) == "-mbar"\) \{(?<body>.*?)\n\s*\}',
  [Text.RegularExpressions.RegexOptions]::Singleline)
if (-not $mbarEnableBlock.Success) {
  throw "mbar_tica_metadata_smoke: could not find the -mbar parser block"
}
if ($mbarEnableBlock.Groups['body'].Value -match 'msm\.flag_(tica|msm)') {
  throw "mbar_tica_metadata_smoke: -mbar must defer TICA/MSM activation to the resolved MBAR domain"
}
# the option families return option_consumed() where the single loop had `continue;`
if ($parser -notmatch '(?s)const int arg_file_status = parse_arg_file\(\);.*?if \(arg_file_status != titan_cli::kParseSuccess\).*?return (?:option_exit\()?arg_file_status\)?;.*?\+\+i;\s*(?:continue;|return option_consumed\(\);)') {
  throw "mbar_tica_metadata_smoke: command-line overrides after -a must remain reachable"
}

$ticaExe = Join-Path $WorkRoot "titan_tica_estimator_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc -Isrc\gmxtraj\include `
  test\titan_tica_estimator_smoke.cpp src\titan_tica_estimator.cpp `
  src\gmxtraj\xdrfile.cpp src\gmxtraj\xdrfile_xtc.cpp `
  -o $ticaExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: failed to compile TICA estimator smoke"
}
& $ticaExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: TICA estimator smoke failed"
}

$mbarExe = Join-Path $WorkRoot "titan_mbar_weights_smoke.exe"
& $Compiler -std=c++23 -O2 -Isrc `
  test\titan_mbar_weights_smoke.cpp src\titan_mbar_weights.cpp `
  -o $mbarExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: failed to compile MBAR weights smoke"
}
& $mbarExe
if ($LASTEXITCODE -ne 0) {
  throw "mbar_tica_metadata_smoke: MBAR weights smoke failed"
}

Write-Host "PASS mbar_tica_metadata_smoke"
