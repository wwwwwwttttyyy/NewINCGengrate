param(
    [string]$Message = "",
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

if ([string]::IsNullOrWhiteSpace($Message)) {
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $Message = "sync update $stamp"
}

Write-Host "Repository: $repoRoot"
Write-Host "Staging tracked and unignored changes ..."
Run-Git add -A

& git diff --cached --quiet
if ($LASTEXITCODE -eq 1) {
    Write-Host "Creating commit: $Message"
    Run-Git commit -m $Message
} elseif ($LASTEXITCODE -eq 0) {
    Write-Host "No staged changes to commit."
} else {
    throw "git diff --cached --quiet failed with exit code $LASTEXITCODE"
}

Write-Host "Rebasing on latest $Remote/$Branch before push ..."
Run-Git fetch $Remote
Run-Git pull --rebase --autostash $Remote $Branch

Write-Host "Pushing to $Remote/$Branch ..."
Run-Git push $Remote $Branch

Write-Host ""
Write-Host "Upload complete."
Run-Git status -sb
