# chemcalculator

This repository is being rebuilt from a single-file Swing prototype into a modular chemistry desktop app.

## Modules

- `chem-domain`: immutable public model types
- `chem-data`: bundled element data, repository implementations, topology helpers
- `chem-engine`: parser, balancer, nomenclature, stoichiometry, verification harness
- `chem-structure-cdk`: CDK-backed structure normalization, formula resolution, and Molfile/SVG viewer payload generation
- `chem-import`: offline import pipeline and SQLite schema/bootstrap entry point
- `chem-ui-desktop`: JavaFX desktop shell with tool tabs, a WebView-backed periodic table explorer, rich balancing traces, and an embedded Chromium 3Dmol.js structure viewer

## Current Status

- The legacy prototype remains in `ChemistryTool.java`.
- The new foundation is the source of truth going forward.
- Core services are designed to avoid silent guessing. Unsupported cases must return warnings.
- Formula-only structure viewing is offline-first and best-effort. The desktop viewer resolves formulas through a bundled compound index that includes all elements plus a larger catalog of common compounds, shows explicit candidate choices for ambiguous formulas such as `C2H6O`, and supports direct SMILES, Molfile/SDF, and InChI input for arbitrary explicit structures.

## Verification

The desktop app targets **JDK 21 LTS**. Gradle toolchains are configured through Foojay so `gradlew` can provision JDK 21 automatically when it is not already installed.

- Preferred when Gradle is available: `.\gradlew.bat verifyCore`
- Desktop app: `.\gradlew.bat :chem-ui-desktop:run`
- Manual fallback: compile the core modules and run `com.chemcalc.engine.EngineVerificationMain`
- Scripted fallback in this repo: `powershell -ExecutionPolicy Bypass -File .\scripts\verify-core.ps1`
