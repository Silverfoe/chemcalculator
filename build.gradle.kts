plugins {
    base
}

group = "com.chemcalc"
version = "0.1.0-SNAPSHOT"

subprojects {
    group = rootProject.group
    version = rootProject.version

    plugins.withId("java") {
        extensions.configure<org.gradle.api.plugins.JavaPluginExtension>("java") {
            toolchain {
                languageVersion = org.gradle.jvm.toolchain.JavaLanguageVersion.of(21)
            }
        }
    }
}

tasks.register("verifyCore") {
    group = "verification"
    description = "Compiles and runs the chemistry verification suites."
    dependsOn(":chem-engine:runManualVerification", ":chem-structure-cdk:runStructureVerification")
}
