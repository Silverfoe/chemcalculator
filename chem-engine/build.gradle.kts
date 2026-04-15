plugins {
    `java-library`
}

dependencies {
    api(project(":chem-domain"))
    implementation(project(":chem-data"))
}

tasks.withType<JavaExec>().configureEach {
    javaLauncher.set(javaToolchains.launcherFor(java.toolchain))
}

tasks.register<JavaExec>("runManualVerification") {
    group = "verification"
    description = "Runs the core manual verification harness without external test dependencies."
    classpath = sourceSets["test"].runtimeClasspath
    mainClass.set("com.chemcalc.engine.EngineVerificationMain")
}
