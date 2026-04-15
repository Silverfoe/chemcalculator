plugins {
    application
    id("org.openjfx.javafxplugin") version "0.1.0"
}

javafx {
    version = "21.0.5"
    modules = listOf("javafx.controls", "javafx.swing", "javafx.web")
}

dependencies {
    implementation(project(":chem-domain"))
    implementation(project(":chem-data"))
    implementation(project(":chem-engine"))
    implementation(project(":chem-structure-cdk"))
    implementation("me.friwi:jcefmaven:143.0.14")
}

tasks.withType<JavaExec>().configureEach {
    javaLauncher.set(javaToolchains.launcherFor(java.toolchain))
}

application {
    mainClass = "com.chemcalc.ui.ChemistrySuperApp"
}
