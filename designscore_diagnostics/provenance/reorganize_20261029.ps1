$ErrorActionPreference = 'Stop'
$taskRepo = [IO.Path]::GetFullPath((Get-Location).Path).TrimEnd('\')
if ((Split-Path $taskRepo -Leaf) -ne 'hchs_mixedmodels_simulations') { throw 'Run from the submodule root.' }
$taskDestination = Join-Path $taskRepo 'designscore_diagnostics'
$taskManifest = Join-Path $taskDestination 'provenance\relocation_manifest.json'
if (Test-Path -LiteralPath $taskManifest) { throw 'Relocation already recorded; refusing to repeat.' }
$taskMoves = [Collections.Generic.List[object]]::new()
function Assert-WithinRepo([string]$Path) {
    $absolute = [IO.Path]::GetFullPath($Path)
    if (-not $absolute.StartsWith($taskRepo + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Outside repository: $absolute" }
    return $absolute
}
function Add-Relocation([string]$Source, [string]$Destination) {
    $sourcePath = Assert-WithinRepo (Join-Path $taskRepo $Source)
    $destinationPath = Assert-WithinRepo (Join-Path $taskDestination $Destination)
    if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) { throw "Missing source: $Source" }
    if (Test-Path -LiteralPath $destinationPath) { throw "Destination exists: $Destination" }
    $taskMoves.Add([ordered]@{old_path=$Source.Replace('\','/'); new_path=$Destination.Replace('\','/');
        sha256_before=(Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash.ToLower();
        bytes=(Get-Item -LiteralPath $sourcePath).Length})
}
$taskSources = @(
 'scripts/stata/scen3_designscore.do', 'scripts/stata/variance_designscore.mata',
 'scripts/stata/run_designscore_local.do', 'scripts/evaluation/summarize_designscore.py',
 'scripts/execution/extend_designscore_1000.py', 'scripts/preparation/prepare_designscore_frame.py',
 'scripts/preparation/recover_designscore_counts.R', 'scripts/tests/test_designscore_algebra.py',
 'scripts/tests/test_designscore_summary.py')
foreach ($source in $taskSources) { Add-Relocation $source ($source -replace '^scripts/', 'code/') }
Add-Relocation 'scripts/stata/SCEN3_DESIGNSCORE.md' 'README.md'
Get-ChildItem -LiteralPath (Join-Path $taskRepo 'data/design_frames') -File | ForEach-Object {
    Add-Relocation ('data/design_frames/' + $_.Name) ('data/design_frames/' + $_.Name)
}
Get-ChildItem -LiteralPath (Join-Path $taskRepo 'data/estimates/oct2026/scen3_designscore/exact') -File | ForEach-Object {
    $target = if ($_.Extension -eq '.log') { 'logs/oct2026/exact/' }
              elseif ($_.Name -like 'local_run_*') { 'logs/local_batches/' }
              else { 'data/estimates/oct2026/exact/' }
    Add-Relocation ('data/estimates/oct2026/scen3_designscore/exact/' + $_.Name) ($target + $_.Name)
}
foreach ($name in @('coverage_scen3_designscore_exact.csv', 'audit_scen3_designscore_exact.csv', 'provenance_scen3_designscore_exact.json')) {
    Add-Relocation ('data/summaries/oct2026/' + $name) ('summaries/oct2026/' + $name)
    Add-Relocation ('data/summaries/oct2026/archive_200/' + $name) ('summaries/oct2026/archive_200/' + $name)
}
Get-ChildItem -LiteralPath $taskRepo -File -Filter '*.log' | ForEach-Object {
    Add-Relocation $_.Name ('logs/local_batches/' + $_.Name)
}
foreach ($run in @('designscore_local', 'designscore_1000')) {
    $sourceRoot = Join-Path $taskRepo ('.statatmp/' + $run)
    Get-ChildItem -LiteralPath $sourceRoot -Recurse -File | ForEach-Object {
        $relative = $_.FullName.Substring($sourceRoot.Length + 1).Replace('\','/')
        $target = if ($_.Extension -eq '.log') { 'logs/local_runs/' }
                  elseif ($_.Extension -in @('.do','.py','.R')) { 'code/archive/local_runs/' }
                  elseif ($_.Extension -eq '.tmp') { '.statatmp/archive/' }
                  elseif ($_.Extension -eq '.csv') { 'data/run_artifacts/' }
                  else { 'provenance/local_runs/' }
        Add-Relocation ('.statatmp/' + $run + '/' + $relative) ($target + $run + '/' + $relative)
    }
}
if (($taskMoves.new_path | Sort-Object -Unique).Count -ne $taskMoves.Count) { throw 'Duplicate destinations' }
$taskProtected = [ordered]@{}
foreach ($subfolder in @('scripts','data/summaries','data/population_coefficients')) {
    Get-ChildItem -LiteralPath (Join-Path $taskRepo $subfolder) -Recurse -File | ForEach-Object {
        $relative = $_.FullName.Substring($taskRepo.Length + 1).Replace('\','/')
        if ($relative -notin $taskMoves.old_path) {
            $taskProtected[$relative] = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLower()
        }
    }
}
# Preserve the exact pre-relocation sources, including the two files fingerprinted by the run.
foreach ($source in ($taskSources + @('scripts/stata/SCEN3_DESIGNSCORE.md','README.md','.gitignore'))) {
    $snapshot = Assert-WithinRepo (Join-Path $taskDestination ('provenance/pre_reorganization/' + $source))
    New-Item -ItemType Directory -Path (Split-Path $snapshot -Parent) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $taskRepo $source) -Destination $snapshot
}
$taskCopies = @()
foreach ($population in @('oct2026','dec2024')) {
    $relative = 'data/population_coefficients/population_coef_' + $population + '_hhsub.csv'
    $target = Assert-WithinRepo (Join-Path $taskDestination $relative)
    New-Item -ItemType Directory -Path (Split-Path $target -Parent) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $taskRepo $relative) -Destination $target
    $taskCopies += [ordered]@{old_path=$relative; new_path=$relative; sha256=(Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLower()}
}
$taskRecord = [ordered]@{date='2026-10-29'; original_root=$taskRepo; destination_root=$taskDestination;
    old_paths_relative_to='original_root'; new_paths_relative_to='destination_root';
    status='planned'; moves=@($taskMoves.ToArray()); copies=$taskCopies; protected_files_sha256=$taskProtected}
$taskRecord | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $taskManifest -Encoding utf8
foreach ($entry in $taskMoves) {
    $source = Assert-WithinRepo (Join-Path $taskRepo $entry.old_path)
    $target = Assert-WithinRepo (Join-Path $taskDestination $entry.new_path)
    New-Item -ItemType Directory -Path (Split-Path $target -Parent) -Force | Out-Null
    Move-Item -LiteralPath $source -Destination $target
    if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLower() -ne $entry.sha256_before) { throw "Hash mismatch: $target" }
}
$taskRecord.status = 'moved_and_hash_verified'
$taskRecord | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $taskManifest -Encoding utf8
# Remove only now-empty source directories; never recursively delete any content.
foreach ($subfolder in @('scripts','data',' .statatmp'.Trim())) {
    $parent = Assert-WithinRepo (Join-Path $taskRepo $subfolder)
    Get-ChildItem -LiteralPath $parent -Directory -Recurse | Sort-Object { $_.FullName.Length } -Descending | ForEach-Object {
        $directory = Assert-WithinRepo $_.FullName
        if (-not (Get-ChildItem -LiteralPath $directory -Force | Select-Object -First 1)) { Remove-Item -LiteralPath $directory }
    }
}
Write-Output ("Moved and SHA-256 verified {0} files; preserved {1} original files." -f $taskMoves.Count, $taskProtected.Count)
