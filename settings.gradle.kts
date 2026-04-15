plugins {
    id("org.gradle.toolchains.foojay-resolver-convention") version "1.0.0"
}

rootProject.name = "chemcalculator"

include(
    "chem-domain",
    "chem-data",
    "chem-engine",
    "chem-structure-cdk",
    "chem-import",
    "chem-ui-desktop"
)

dependencyResolutionManagement {
    repositoriesMode.set(RepositoriesMode.FAIL_ON_PROJECT_REPOS)
    repositories {
        mavenCentral()
    }
}
