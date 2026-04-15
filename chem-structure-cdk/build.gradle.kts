plugins {
    `java-library`
}

dependencies {
    implementation(project(":chem-engine"))
    implementation(project(":chem-data"))

    implementation("org.openscience.cdk:cdk-bundle:2.12")
    implementation("org.openscience.cdk:cdk-inchi:2.12")
}

tasks.withType<JavaExec>().configureEach {
    javaLauncher.set(javaToolchains.launcherFor(java.toolchain))
}

tasks.register<JavaExec>("runStructureVerification") {
    group = "verification"
    description = "Runs the CDK structure verification harness."
    classpath = sourceSets["test"].runtimeClasspath
    mainClass.set("com.chemcalc.structurecdk.StructureVerificationMain")
}
