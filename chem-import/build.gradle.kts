plugins {
    application
}

dependencies {
    implementation(project(":chem-data"))
    runtimeOnly("org.xerial:sqlite-jdbc:3.46.1.3")
}

tasks.withType<JavaExec>().configureEach {
    javaLauncher.set(javaToolchains.launcherFor(java.toolchain))
}

application {
    mainClass = "com.chemcalc.importer.ChemDataImportTool"
}
