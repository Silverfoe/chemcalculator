$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$classDir = Join-Path $repoRoot "build\manual\classes"

if (Test-Path $classDir) {
    Remove-Item -Recurse -Force $classDir
}

New-Item -ItemType Directory -Force -Path $classDir | Out-Null

$sourceRoots = @(
    (Join-Path $repoRoot "chem-domain\src\main\java"),
    (Join-Path $repoRoot "chem-data\src\main\java"),
    (Join-Path $repoRoot "chem-engine\src\main\java"),
    (Join-Path $repoRoot "chem-engine\src\test\java")
)

$sources = Get-ChildItem -Recurse -Filter *.java $sourceRoots | ForEach-Object { $_.FullName }

& javac -d $classDir $sources

$classpath = @(
    $classDir,
    (Join-Path $repoRoot "chem-data\src\main\resources")
) -join ";"

& java -cp $classpath com.chemcalc.engine.EngineVerificationMain
