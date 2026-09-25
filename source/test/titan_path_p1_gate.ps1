$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $PSScriptRoot
Push-Location $Root
try {
  $steps = @(
    @('powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', 'test\titan_path_route_smoke.ps1'),
    @('powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', 'test\titan_path_registration_smoke.ps1'),
    @('powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', 'test\titan_path_geometry_smoke.ps1'),
    @('powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', 'test\titan_path_control_smoke.ps1'),
    @('python', '-X', 'utf8', 'test\titan_path_handoff_smoke.py'),
    @('python', '-X', 'utf8', 'test\titan_direction_adapter_smoke.py'),
    @('python', '-X', 'utf8', 'test\titan_path_run_snapshot_smoke.py'),
    @('python', '-X', 'utf8', 'test\titan_path_controller_smoke.py')
  )
  foreach ($step in $steps) {
    & $step[0] $step[1..($step.Count - 1)]
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
  }
  & mingw32-make PROFILE=native all -j1
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
  & mingw32-make PROFILE=native cli-validation-test -j1
  if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
  Write-Output 'PASS titan_path_p1_gate'
} finally {
  Pop-Location
}
