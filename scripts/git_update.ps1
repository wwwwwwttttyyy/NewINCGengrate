param(
    [string]$Remote = "origin",
    [string]$Branch = "main"
)

$ErrorActionPreference = "Stop"

function Run-Git {
    & git @args
    if ($LASTEXITCODE -ne 0) {
        throw "git $($args -join ' ') failed with exit code $LASTEXITCODE"
    }
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir
Set-Location $repoRoot

Run-Git rev-parse --is-inside-work-tree | Out-Null

Write-Host "Repository: $repoRoot"
Write-Host "Updating from $Remote/$Branch ..."

Run-Git fetch $Remote
Run-Git pull --rebase --autostash $Remote $Branch

Write-Host ""
Write-Host "Update complete."
Run-Git status -sb
