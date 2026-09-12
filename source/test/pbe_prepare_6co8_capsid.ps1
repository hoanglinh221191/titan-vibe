[CmdletBinding()]
param(
    [string]$OutputRoot = 'test\pbe_inputs_capsid',
    [string]$AssemblyUrl = 'https://files.rcsb.org/download/6CO8.pdb1',
    [string]$PythonExe = 'C:\Users\ADMIN\AppData\Local\Programs\Python\Python314\python.exe',
    [string]$Pdb2PqrTarget = '.codex_tmp\pdb2pqr-3.7.1-py314',
    [ValidateRange(1,100)][int]$ExpectedModels = 60,
    [ValidateRange(0,100)][int]$MaximumModels = 0,
    [ValidateRange(1,10000000)][int64]$MinimumAssemblyAtoms = 700000,
    [ValidateRange(1,10000000)][int64]$MinimumFullPqrAtoms = 1000000,
    [switch]$InstallPdb2Pqr,
    [switch]$ForceDownload,
    [switch]$ForceSplit,
    [switch]$ForcePrepare
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$RepoRoot = Split-Path -Parent $PSScriptRoot
$OutputRoot = if ([IO.Path]::IsPathRooted($OutputRoot)) {
    $OutputRoot
} else { Join-Path $RepoRoot $OutputRoot }
$Pdb2PqrTarget = if ([IO.Path]::IsPathRooted($Pdb2PqrTarget)) {
    $Pdb2PqrTarget
} else { Join-Path $RepoRoot $Pdb2PqrTarget }
$PythonExe = [IO.Path]::GetFullPath($PythonExe)
if (-not (Test-Path -LiteralPath $PythonExe)) {
    throw "Python executable not found: $PythonExe"
}
New-Item -ItemType Directory -Path $OutputRoot -Force | Out-Null
$modelPdbRoot = Join-Path $OutputRoot 'models_pdb'
$modelPqrRoot = Join-Path $OutputRoot 'models_pqr'
New-Item -ItemType Directory -Path $modelPdbRoot,$modelPqrRoot -Force | Out-Null
$assemblyPath = Join-Path $OutputRoot '6co8_assembly1.pdb'
$combinedPqrPath = Join-Path $OutputRoot '6co8_assembly1_charmm.pqr'
$manifestPath = Join-Path $OutputRoot '6co8_capsid_fixture_manifest.csv'

function Count-Records([string]$Path, [string[]]$Prefixes) {
    $count = 0L
    $reader = [IO.File]::OpenText($Path)
    try {
        while ($null -ne ($line = $reader.ReadLine())) {
            foreach ($prefix in $Prefixes) {
                if ($line.StartsWith($prefix)) { ++$count; break }
            }
        }
    } finally { $reader.Dispose() }
    return $count
}
if ($ForceDownload -or -not (Test-Path -LiteralPath $assemblyPath)) {
    Write-Host "download: $AssemblyUrl"
    $partial = "$assemblyPath.partial"
    if (Test-Path -LiteralPath $partial) { Remove-Item $partial -Force }
    Invoke-WebRequest -Uri $AssemblyUrl -OutFile $partial -UseBasicParsing `
        -TimeoutSec 1800
    Move-Item -LiteralPath $partial -Destination $assemblyPath -Force
}
$assemblyAtoms = Count-Records $assemblyPath @('ATOM  ','HETATM')
$assemblyModels = Count-Records $assemblyPath @('MODEL ')
if ($assemblyAtoms -lt $MinimumAssemblyAtoms -or
    $assemblyModels -ne $ExpectedModels) {
    throw "6CO8 assembly contract failed: atoms=$assemblyAtoms models=$assemblyModels; expected at least $MinimumAssemblyAtoms atoms and exactly $ExpectedModels models."
}
$assemblySha = (Get-FileHash $assemblyPath -Algorithm SHA256).Hash.ToLowerInvariant()

$splitComplete = Join-Path $modelPdbRoot 'split_complete.dat'
if ($ForceSplit -or -not (Test-Path -LiteralPath $splitComplete)) {
    Get-ChildItem $modelPdbRoot -Filter 'model_*.pdb' -ErrorAction SilentlyContinue |
        Remove-Item -Force
    $reader = [IO.File]::OpenText($assemblyPath)
    $writer = $null
    $modelIndex = 0
    try {
        while ($null -ne ($line = $reader.ReadLine())) {
            if ($line.StartsWith('MODEL ')) {
                if ($null -ne $writer) { throw 'Nested MODEL records in 6CO8 assembly.' }
                ++$modelIndex
                $modelPath = Join-Path $modelPdbRoot ('model_{0:D3}.pdb' -f $modelIndex)
                $writer = [IO.StreamWriter]::new($modelPath, $false,
                    [Text.Encoding]::ASCII)
                $writer.WriteLine('HEADER    6CO8 BIOLOGICAL ASSEMBLY 1 MODEL {0}' -f $modelIndex)
            } elseif ($line.StartsWith('ENDMDL')) {
                if ($null -eq $writer) { throw 'ENDMDL without MODEL in 6CO8 assembly.' }
                $writer.WriteLine('END')
                $writer.Dispose()
                $writer = $null
            } elseif ($null -ne $writer -and
                ($line.StartsWith('ATOM  ') -or $line.StartsWith('HETATM') -or
                 $line.StartsWith('TER'))) {
                $writer.WriteLine($line)
            }
        }
    } finally {
        if ($null -ne $writer) { $writer.Dispose() }
        $reader.Dispose()
    }
    if ($modelIndex -ne $ExpectedModels) {
        throw "Split produced $modelIndex models; expected $ExpectedModels."
    }
    Set-Content -LiteralPath $splitComplete -Value @(
        "assembly_sha256=$assemblySha", "models=$modelIndex") -Encoding ASCII
}

$pdb2pqrMarker = Join-Path $Pdb2PqrTarget 'pdb2pqr\__init__.py'
if ($InstallPdb2Pqr -and -not (Test-Path -LiteralPath $pdb2pqrMarker)) {
    New-Item -ItemType Directory -Path $Pdb2PqrTarget -Force | Out-Null
    & $PythonExe -m pip install --disable-pip-version-check `
        --target $Pdb2PqrTarget 'pdb2pqr==3.7.1'
    if ($LASTEXITCODE -ne 0) { throw 'Unable to install pdb2pqr==3.7.1.' }
}
if (-not (Test-Path -LiteralPath $pdb2pqrMarker)) {
    throw "Pinned PDB2PQR runtime missing: $Pdb2PqrTarget. Rerun with -InstallPdb2Pqr."
}
$env:PYTHONPATH = $Pdb2PqrTarget

$modelLimit = if ($MaximumModels -gt 0) {
    [Math]::Min($MaximumModels, $ExpectedModels)
} else { $ExpectedModels }
$modelRows = [Collections.Generic.List[object]]::new()
for ($model = 1; $model -le $modelLimit; ++$model) {
    $modelPdb = Join-Path $modelPdbRoot ('model_{0:D3}.pdb' -f $model)
    $modelPqr = Join-Path $modelPqrRoot ('model_{0:D3}.pqr' -f $model)
    $modelLog = Join-Path $modelPqrRoot ('model_{0:D3}.log' -f $model)
    $modelErr = Join-Path $modelPqrRoot ('model_{0:D3}.err.log' -f $model)
    $reuse = -not $ForcePrepare -and (Test-Path -LiteralPath $modelPqr) -and
        (Count-Records $modelPqr @('ATOM  ','HETATM')) -gt 0
    if ($reuse) {
        Write-Host ('resume: model {0:D3}/{1}' -f $model,$modelLimit)
    } else {
        Write-Host ('prepare: model {0:D3}/{1}' -f $model,$modelLimit)
        if (Test-Path -LiteralPath $modelPqr) { Remove-Item $modelPqr -Force }
        $process = Start-Process -FilePath $PythonExe -ArgumentList @(
            '-m','pdb2pqr','--ff=CHARMM','--ffout=CHARMM','--keep-chain',
            '--nodebump','--noopt','--drop-water',$modelPdb,$modelPqr
        ) -WorkingDirectory $RepoRoot -NoNewWindow -Wait -PassThru `
            -RedirectStandardOutput $modelLog -RedirectStandardError $modelErr
        if ($process.ExitCode -ne 0 -or
            -not (Test-Path -LiteralPath $modelPqr)) {
            throw "PDB2PQR failed for model $model. See $modelLog and $modelErr"
        }
    }
    $inputCount = Count-Records $modelPdb @('ATOM  ','HETATM')
    $outputCount = Count-Records $modelPqr @('ATOM  ','HETATM')
    if ($outputCount -lt $inputCount) {
        throw "PDB2PQR model $model lost atoms: input=$inputCount output=$outputCount"
    }
    $modelRows.Add([pscustomobject]@{
        model=$model; input_atoms=$inputCount; output_atoms=$outputCount
        pqr_sha256=(Get-FileHash $modelPqr -Algorithm SHA256).Hash.ToLowerInvariant()
        status='ok'
    })
}

if ($modelLimit -eq $ExpectedModels) {
    Write-Host "combine: $ExpectedModels model PQR files"
    $combiner = Join-Path $PSScriptRoot 'pbe_prepare_6co8_capsid.py'
    & $PythonExe $combiner --model-root $modelPqrRoot `
        --output $combinedPqrPath --models $ExpectedModels
    if ($LASTEXITCODE -ne 0) { throw '6CO8 PQR combination failed.' }
    $pqrAtoms = Count-Records $combinedPqrPath @('ATOM  ','HETATM')
    if ($pqrAtoms -lt $MinimumFullPqrAtoms) {
        throw "Combined PQR has only $pqrAtoms atoms; expected at least $MinimumFullPqrAtoms."
    }
    @([pscustomobject]@{
        fixture='6co8'; pdb_id='6CO8'; assembly='biological assembly 1'
        source_url='https://www.rcsb.org/structure/6CO8'
        download_url=$AssemblyUrl; assembly_models=$assemblyModels
        assembly_atoms=$assemblyAtoms; pqr_atoms=$pqrAtoms
        pdb2pqr_version='3.7.1'; force_field='CHARMM'
        preparation='per-model --ff=CHARMM --ffout=CHARMM --keep-chain --nodebump --noopt --drop-water; concatenate 60 model PQRs'
        pdb_sha256=$assemblySha
        pqr_sha256=(Get-FileHash $combinedPqrPath -Algorithm SHA256).Hash.ToLowerInvariant()
        pqr_path=$combinedPqrPath; status='ok'
    }) | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
    Write-Host "6CO8 fixture ready: atoms=$pqrAtoms path=$combinedPqrPath"
} else {
    $diagnosticPath = Join-Path $OutputRoot '6co8_partial_model_manifest.csv'
    $modelRows.ToArray() | Export-Csv -LiteralPath $diagnosticPath `
        -NoTypeInformation -Encoding UTF8
    Write-Host "Partial preparation complete ($modelLimit/$ExpectedModels models): $diagnosticPath"
}
