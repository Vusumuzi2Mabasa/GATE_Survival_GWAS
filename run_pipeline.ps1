param([switch]$DryRun)
$ErrorActionPreference = 'Stop'
$projectWindowsPath = (Get-Item -LiteralPath $PSScriptRoot).FullName.Replace('\', '/')
if ($projectWindowsPath -notmatch '^([A-Za-z]):/(.*)$') {
    throw 'This launcher expects a normal Windows drive path.'
}
$projectLinuxPath = '/mnt/' + $Matches[1].ToLowerInvariant() + '/' + $Matches[2]
if ($DryRun) {
    wsl -d Ubuntu-24.04 --cd $projectLinuxPath -- bash run_snakemake.sh --dry-run
} else {
    wsl -d Ubuntu-24.04 --cd $projectLinuxPath -- bash run_snakemake.sh
}
if ($LASTEXITCODE -ne 0) { throw 'Pipeline failed. Inspect project logs.' }
