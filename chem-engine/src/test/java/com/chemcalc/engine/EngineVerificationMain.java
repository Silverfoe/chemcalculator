package com.chemcalc.engine;

import com.chemcalc.data.BundledPeriodicTableRepository;
import com.chemcalc.domain.BalancedReaction;
import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.NomenclatureResult;
import com.chemcalc.domain.Reaction;
import com.chemcalc.domain.Species;
import java.util.List;
import java.util.Map;

public final class EngineVerificationMain {
    private EngineVerificationMain() {
    }

    public static void main(String[] args) {
        BundledPeriodicTableRepository repository = new BundledPeriodicTableRepository();
        FormulaParser formulaParser = new DefaultFormulaParser();
        ReactionParser reactionParser = new DefaultReactionParser(formulaParser);
        BalancerService balancerService = new LinearEquationBalancer();
        StoichiometryService stoichiometryService = new DefaultStoichiometryService(repository);
        NomenclatureService nomenclatureService = new RuleBasedNomenclatureService(repository);

        verifyData(repository.findAll());
        verifyParser(formulaParser);
        verifyBalancing(reactionParser, balancerService);
        verifyStoichiometry(formulaParser, stoichiometryService);
        verifyNomenclature(formulaParser, nomenclatureService);

        System.out.println("Core verification passed.");
    }

    private static void verifyData(List<ElementRecord> elements) {
        require(elements.size() == 118, "Expected 118 bundled element records");
        require(elements.stream().allMatch(element -> !element.symbol().isBlank()), "All elements need symbols");
        require(elements.stream().allMatch(element -> !element.name().isBlank()), "All elements need names");
        require(elements.stream().allMatch(element -> element.group() != null), "All elements need groups");
        require(elements.stream().allMatch(element -> element.period() >= 1), "All elements need periods");
    }

    private static void verifyParser(FormulaParser parser) {
        assertSpecies(parser.parseSpecies("Fe2(SO4)3"), Map.of("Fe", 2, "S", 3, "O", 12), 0, "Fe2(SO4)3");
        assertSpecies(parser.parseSpecies("CuSO4·5H2O"), Map.of("Cu", 1, "S", 1, "O", 9, "H", 10), 0, "CuSO4·5H2O");
        assertSpecies(parser.parseSpecies("NH4+"), Map.of("N", 1, "H", 4), 1, "NH4");
        Species leadNitrate = parser.parseSpecies("Pb(NO3)2(aq)");
        assertSpecies(leadNitrate, Map.of("Pb", 1, "N", 2, "O", 6), 0, "Pb(NO3)2");
        require(leadNitrate.phase().token().equals("aq"), "Expected aqueous phase");
        assertSpecies(parser.parseSpecies("[Fe(CN)6]4-"), Map.of("Fe", 1, "C", 6, "N", 6), -4, "[Fe(CN)6]");
    }

    private static void verifyBalancing(ReactionParser parser, BalancerService balancerService) {
        verifyBalanced(parser, balancerService, "H2 + O2 -> H2O", "2 H2 + O2 -> 2 H2O");
        verifyBalanced(parser, balancerService, "Fe + O2 -> Fe2O3", "4 Fe + 3 O2 -> 2 Fe2O3");
        verifyBalanced(parser, balancerService, "C2H6 + O2 -> CO2 + H2O", "2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O");
        verifyBalanced(parser, balancerService, "KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O", "2 KMnO4 + 16 HCl -> 2 KCl + 2 MnCl2 + 5 Cl2 + 8 H2O");
        verifyBalanced(parser, balancerService, "Cr2O7^2- + Fe^2+ + H+ -> Cr^3+ + Fe^3+ + H2O", "Cr2O7^2- + 6 Fe2+ + 14 H+ -> 2 Cr3+ + 6 Fe3+ + 7 H2O");
    }

    private static void verifyStoichiometry(FormulaParser parser, StoichiometryService stoichiometryService) {
        double water = stoichiometryService.molarMass(parser.parseSpecies("H2O")).totalMass();
        require(Math.abs(water - 18.015) < 0.01, "Water molar mass should be about 18.015 g/mol");
        double oxygenPercent = stoichiometryService.percentComposition(parser.parseSpecies("H2O")).entries().stream()
                .filter(entry -> entry.symbol().equals("O"))
                .findFirst()
                .orElseThrow()
                .percent();
        require(oxygenPercent > 88 && oxygenPercent < 89, "Water oxygen mass percent should be about 88.8%");
    }

    private static void verifyNomenclature(FormulaParser parser, NomenclatureService nomenclatureService) {
        NomenclatureResult sodiumChloride = nomenclatureService.name(parser.parseSpecies("NaCl"));
        require(sodiumChloride.supported(), "NaCl should be supported");
        require("Sodium chloride".equals(sodiumChloride.displayName()), "Unexpected NaCl name: " + sodiumChloride.displayName());

        NomenclatureResult dinitrogenPentoxide = nomenclatureService.name(parser.parseSpecies("N2O5"));
        require(dinitrogenPentoxide.supported(), "N2O5 should be supported");
        require("Dinitrogen pentoxide".equals(dinitrogenPentoxide.displayName()), "Unexpected N2O5 name");

        NomenclatureResult unsupported = nomenclatureService.name(parser.parseSpecies("C6H12O6"));
        require(!unsupported.supported(), "Glucose should not be named confidently by the bootstrap nomenclature engine");
    }

    private static void verifyBalanced(ReactionParser parser, BalancerService balancerService, String input, String expected) {
        Reaction reaction = parser.parse(input);
        BalancedReaction balanced = balancerService.balance(reaction);
        require(expected.equals(balanced.format()), "Expected '" + expected + "' but got '" + balanced.format() + "'");
        require(!balanced.trace().quantities().isEmpty(), "Expected conserved quantities trace for " + input);
        require(!balanced.trace().initialMatrix().values().isEmpty(), "Expected initial matrix trace for " + input);
        require(!balanced.trace().rowOperations().isEmpty(), "Expected row operations trace for " + input);
        require(!balanced.trace().solution().nullspaceVector().isEmpty(), "Expected solution trace for " + input);
        require(!balanced.trace().normalization().coefficients().isEmpty(), "Expected normalization trace for " + input);
        require(!balanced.trace().conservationChecks().isEmpty(), "Expected conservation checks for " + input);
        require(balanced.trace().conservationChecks().stream()
                        .allMatch(check -> check.reactantTotal() == check.productTotal()),
                "Expected balanced conservation checks for " + input);
    }

    private static void assertSpecies(Species species, Map<String, Integer> expectedComposition, int expectedCharge, String expectedCanonicalFormula) {
        require(species.composition().equals(expectedComposition), "Unexpected composition for " + expectedCanonicalFormula + ": " + species.composition());
        require(species.charge() == expectedCharge, "Unexpected charge for " + expectedCanonicalFormula);
        require(species.canonicalFormula().equals(expectedCanonicalFormula), "Unexpected canonical formula for " + expectedCanonicalFormula);
    }

    private static void require(boolean condition, String message) {
        if (!condition) {
            throw new IllegalStateException(message);
        }
    }
}
