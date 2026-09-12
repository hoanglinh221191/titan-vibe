param(
    [Parameter(Mandatory = $true)][string]$Protocol,
    [Parameter(Mandatory = $true)][string]$OutputRoot
)

$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$protocolPath = (Resolve-Path $Protocol).Path
$protocolData = Get-Content -LiteralPath $protocolPath -Raw | ConvertFrom-Json
$titanExe = Join-Path $repo 'titan.exe'

if (Test-Path -LiteralPath $OutputRoot) {
    throw "Output already exists: $OutputRoot"
}
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null

$ledger = [Collections.Generic.List[object]]::new()
Push-Location $repo
try {
    foreach ($item in $protocolData.matrix) {
        $command = [string]$item.command
        $command = $command.Replace('${TITAN_EXE}', $titanExe.Replace("'", "''"))
        $log = Join-Path $OutputRoot "$($item.name).log"
        $started = Get-Date
        $savedPreference = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        & powershell.exe -NoProfile -ExecutionPolicy Bypass -Command $command *> $log
        $exitCode = $LASTEXITCODE
        $ErrorActionPreference = $savedPreference
        $binaryHash = if (Test-Path -LiteralPath $titanExe) {
            (Get-FileHash -Algorithm SHA256 -LiteralPath $titanExe).Hash.ToLowerInvariant()
        } else {
            ''
        }
        $ledger.Add([pscustomobject]@{
            name = [string]$item.name
            command = $command
            exit_code = $exitCode
            elapsed_seconds = ((Get-Date) - $started).TotalSeconds
            log = (Resolve-Path $log).Path
            log_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $log).Hash.ToLowerInvariant()
            titan_sha256 = $binaryHash
        })
        if ($exitCode -ne 0) {
            throw "$($item.name) failed with exit code $exitCode"
        }
    }
} finally {
    Pop-Location
    $ledger | Export-Csv (Join-Path $OutputRoot 'test_ledger.csv') -NoTypeInformation -Encoding utf8
}

$ledger | Format-Table name, exit_code, elapsed_seconds -AutoSize
