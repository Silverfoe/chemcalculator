package com.chemcalc.engine;

import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.ExplanationStep;
import com.chemcalc.domain.NomenclatureResult;
import com.chemcalc.domain.PeriodicTableRepository;
import com.chemcalc.domain.Species;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

public final class RuleBasedNomenclatureService implements NomenclatureService {
    private static final Map<String, String> ANION_NAMES = Map.ofEntries(
            Map.entry("H", "hydride"),
            Map.entry("F", "fluoride"),
            Map.entry("Cl", "chloride"),
            Map.entry("Br", "bromide"),
            Map.entry("I", "iodide"),
            Map.entry("O", "oxide"),
            Map.entry("S", "sulfide"),
            Map.entry("N", "nitride"),
            Map.entry("P", "phosphide"),
            Map.entry("C", "carbide"));

    private static final Map<String, String> PREFIXES = Map.ofEntries(
            Map.entry("1", "mono"),
            Map.entry("2", "di"),
            Map.entry("3", "tri"),
            Map.entry("4", "tetra"),
            Map.entry("5", "penta"),
            Map.entry("6", "hexa"),
            Map.entry("7", "hepta"),
            Map.entry("8", "octa"),
            Map.entry("9", "nona"),
            Map.entry("10", "deca"));

    private static final Map<String, String> EXACT_ACIDS = Map.of(
            "HCl", "Hydrochloric acid",
            "HBr", "Hydrobromic acid",
            "HI", "Hydroiodic acid",
            "HF", "Hydrofluoric acid",
            "HNO3", "Nitric acid",
            "HNO2", "Nitrous acid",
            "H2SO4", "Sulfuric acid",
            "H2SO3", "Sulfurous acid",
            "H2CO3", "Carbonic acid",
            "H3PO4", "Phosphoric acid");

    private static final Map<String, String> POLYATOMIC_SUFFIXES = Map.of(
            "OH", "hydroxide",
            "NO3", "nitrate",
            "NO2", "nitrite",
            "SO4", "sulfate",
            "SO3", "sulfite",
            "CO3", "carbonate",
            "PO4", "phosphate",
            "C2H3O2", "acetate");

    private static final Set<String> FIXED_CHARGE_CATIONS = Set.of("Li", "Na", "K", "Rb", "Cs", "Fr", "Be", "Mg", "Ca", "Sr", "Ba", "Ra", "Al", "Zn", "Ag");
    private final PeriodicTableRepository periodicTableRepository;
    private final FormulaParser formulaParser = new DefaultFormulaParser();

    public RuleBasedNomenclatureService(PeriodicTableRepository periodicTableRepository) {
        this.periodicTableRepository = periodicTableRepository;
    }

    @Override
    public NomenclatureResult name(Species species) {
        String formula = species.canonicalFormula();
        if (formula.isBlank()) {
            return NomenclatureResult.unsupported("Species has no formula");
        }

        if (formula.contains("·")) {
            return nameHydrate(species);
        }
        if (EXACT_ACIDS.containsKey(formula) && species.charge() == 0) {
            return new NomenclatureResult(EXACT_ACIDS.get(formula), true, List.of(new ExplanationStep("Exact acid match", "Matched a curated acid rule.")), List.of());
        }
        if (species.composition().size() == 1) {
            String symbol = species.composition().keySet().iterator().next();
            return periodicTableRepository.findBySymbol(symbol)
                    .map(element -> new NomenclatureResult(element.name(), true, List.of(), List.of()))
                    .orElseGet(() -> NomenclatureResult.unsupported("Unknown monoatomic species: " + symbol));
        }
        if (species.composition().size() == 2) {
            return nameBinary(species);
        }
        NomenclatureResult polyatomic = nameCommonPolyatomic(species);
        if (polyatomic.supported()) {
            return polyatomic;
        }
        return NomenclatureResult.unsupported("The current naming engine only guarantees binary, hydrate, curated acid, and common polyatomic patterns.");
    }

    private NomenclatureResult nameHydrate(Species species) {
        String[] segments = species.canonicalFormula().split("·");
        if (segments.length != 2) {
            return NomenclatureResult.unsupported("Only simple hydrate patterns are supported right now.");
        }
        String hydrateSegment = segments[1];
        int waterCount = 1;
        String waterFormula = hydrateSegment;
        int index = 0;
        while (index < hydrateSegment.length() && Character.isDigit(hydrateSegment.charAt(index))) {
            index++;
        }
        if (index > 0) {
            waterCount = Integer.parseInt(hydrateSegment.substring(0, index));
            waterFormula = hydrateSegment.substring(index);
        }
        if (!"H2O".equals(waterFormula)) {
            return NomenclatureResult.unsupported("Hydrate naming is only implemented for water adducts.");
        }
        Species base = formulaParser.parseSpecies(segments[0]);
        NomenclatureResult baseName = name(base);
        if (!baseName.supported()) {
            return baseName;
        }
        return new NomenclatureResult(
                baseName.displayName() + " " + prefixFor(waterCount) + "hydrate",
                true,
                List.of(new ExplanationStep("Hydrate rule", "Detected a water adduct and applied the standard hydrate suffix.")),
                List.of());
    }

    private NomenclatureResult nameBinary(Species species) {
        List<String> symbols = new ArrayList<>(species.composition().keySet());
        String first = symbols.get(0);
        String second = symbols.get(1);
        ElementRecord firstElement = periodicTableRepository.findBySymbol(first).orElse(null);
        ElementRecord secondElement = periodicTableRepository.findBySymbol(second).orElse(null);
        if (firstElement == null || secondElement == null) {
            return NomenclatureResult.unsupported("Binary compound contains an unknown element.");
        }
        boolean firstMetal = isMetal(firstElement);
        boolean secondMetal = isMetal(secondElement);
        if (firstMetal && !secondMetal) {
            int cationCount = species.composition().get(first);
            int anionCount = species.composition().get(second);
            Integer anionCharge = typicalAnionCharge(second);
            if (anionCharge == null) {
                return NomenclatureResult.unsupported("No safe ionic charge rule exists for " + second + ".");
            }
            int cationChargeTimesCount = -anionCharge * anionCount;
            if (cationChargeTimesCount % cationCount != 0) {
                return NomenclatureResult.unsupported("Cannot derive an integer oxidation state for the cation.");
            }
            int cationCharge = cationChargeTimesCount / cationCount;
            String cationName = firstElement.name();
            if (!FIXED_CHARGE_CATIONS.contains(first)) {
                cationName += " (" + toRoman(cationCharge) + ")";
            }
            String anionName = ANION_NAMES.getOrDefault(second, secondElement.name().toLowerCase(Locale.ROOT) + "ide");
            return new NomenclatureResult(capitalize(cationName + " " + anionName), true, List.of(), List.of());
        }
        if (!firstMetal && !secondMetal) {
            int firstCount = species.composition().get(first);
            int secondCount = species.composition().get(second);
            String firstName = firstElement.name().toLowerCase(Locale.ROOT);
            String secondName = ANION_NAMES.getOrDefault(second, secondElement.name().toLowerCase(Locale.ROOT) + "ide");
            String firstPrefix = firstCount == 1 ? "" : prefixFor(firstCount);
            String name = (firstPrefix + firstName + " " + prefixFor(secondCount) + secondName).trim();
            name = name.replace("monooxide", "monoxide");
            name = name.replace("pentaoxide", "pentoxide");
            name = name.replace("tetraoxide", "tetroxide");
            return new NomenclatureResult(capitalize(name), true, List.of(), List.of());
        }
        return NomenclatureResult.unsupported("Binary metallic clusters are not named automatically.");
    }

    private NomenclatureResult nameCommonPolyatomic(Species species) {
        String formula = species.canonicalFormula();
        for (Map.Entry<String, String> entry : POLYATOMIC_SUFFIXES.entrySet()) {
            String ionFormula = entry.getKey();
            if (!formula.endsWith(ionFormula)) {
                continue;
            }
            Map<String, Integer> remaining = new LinkedHashMap<>(species.composition());
            for (Map.Entry<String, Integer> ionComponent : formulaParser.parseSpecies(ionFormula).composition().entrySet()) {
                remaining.computeIfPresent(ionComponent.getKey(), (ignored, count) -> count - ionComponent.getValue());
            }
            remaining.entrySet().removeIf(item -> item.getValue() <= 0);
            if (remaining.size() != 1) {
                continue;
            }
            String cationSymbol = remaining.keySet().iterator().next();
            ElementRecord cation = periodicTableRepository.findBySymbol(cationSymbol).orElse(null);
            if (cation == null) {
                continue;
            }
            return new NomenclatureResult(capitalize(cation.name() + " " + entry.getValue()), true, List.of(), List.of());
        }
        return NomenclatureResult.unsupported("No supported common polyatomic pattern matched.");
    }

    private boolean isMetal(ElementRecord element) {
        return element.category().contains("Metal") && !element.category().contains("Nonmetal");
    }

    private Integer typicalAnionCharge(String symbol) {
        return switch (symbol) {
            case "F", "Cl", "Br", "I", "H" -> -1;
            case "O", "S" -> -2;
            case "N", "P" -> -3;
            default -> null;
        };
    }

    private String prefixFor(int count) {
        return PREFIXES.getOrDefault(Integer.toString(count), count == 1 ? "mono" : Integer.toString(count));
    }

    private String capitalize(String value) {
        return value.substring(0, 1).toUpperCase(Locale.ROOT) + value.substring(1);
    }

    private String toRoman(int value) {
        return switch (value) {
            case 1 -> "I";
            case 2 -> "II";
            case 3 -> "III";
            case 4 -> "IV";
            case 5 -> "V";
            case 6 -> "VI";
            case 7 -> "VII";
            case 8 -> "VIII";
            default -> Integer.toString(value);
        };
    }
}
