package com.chemcalc.data;

import com.chemcalc.domain.CompoundStructureKind;
import com.chemcalc.domain.CompoundStructureRecord;
import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.SourceReference;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

final class CompoundStructureCatalog {
    private CompoundStructureCatalog() {
    }

    static List<CompoundStructureRecord> allRecords() {
        SourceRegistry sources = SourceRegistry.load();
        SourceReference internal = sources.lookup("INTERNAL");

        Map<String, CompoundStructureRecord> records = new LinkedHashMap<>();
        loadSeedRecords(sources).forEach(record -> records.put(record.id(), record));
        generateElementRecords(internal).forEach(record -> records.put(record.id(), record));
        generateElementMoleculeRecords(internal).forEach(record -> records.put(record.id(), record));
        generateCuratedDiscreteRecords(internal).forEach(record -> records.put(record.id(), record));
        generateRepresentativeFormulaUnits(internal).forEach(record -> records.put(record.id(), record));
        return List.copyOf(records.values());
    }

    private static List<CompoundStructureRecord> loadSeedRecords(SourceRegistry sources) {
        List<CompoundStructureRecord> loaded = new ArrayList<>();
        try (InputStream input = CompoundStructureCatalog.class.getResourceAsStream("/data/compound_structures.psv")) {
            if (input == null) {
                throw new IllegalStateException("Bundled compound structure snapshot not found");
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8))) {
                if (reader.readLine() == null) {
                    return List.of();
                }
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.isBlank() || line.startsWith("#")) {
                        continue;
                    }
                    String[] parts = line.split("\\|", -1);
                    if (parts.length < 7) {
                        continue;
                    }

                    String id;
                    String formula;
                    int charge;
                    String preferredName;
                    String canonicalSmiles;
                    String inchi;
                    String inchiKey;
                    CompoundStructureKind structureKind;
                    String sourceId;

                    if (parts.length >= 9) {
                        id = parts[0];
                        formula = parts[1];
                        charge = Integer.parseInt(parts[2]);
                        preferredName = parts[3];
                        canonicalSmiles = parts[4];
                        inchi = parts[5];
                        inchiKey = parts[6];
                        structureKind = parseStructureKind(parts[7]);
                        sourceId = parts[8];
                    } else {
                        formula = parts[0];
                        charge = Integer.parseInt(parts[1]);
                        preferredName = parts[2];
                        canonicalSmiles = parts[3];
                        inchi = parts[4];
                        inchiKey = parts[5];
                        structureKind = CompoundStructureKind.DISCRETE_MOLECULE;
                        sourceId = parts[6];
                        id = fallbackSeedId(formula, charge, preferredName);
                    }

                    loaded.add(new CompoundStructureRecord(
                            id,
                            formula,
                            charge,
                            preferredName,
                            canonicalSmiles,
                            inchi,
                            inchiKey,
                            structureKind,
                            sources.lookup(sourceId)));
                }
            }
        } catch (IOException exception) {
            throw new IllegalStateException("Failed to load bundled compound structure data", exception);
        }
        return List.copyOf(loaded);
    }

    private static List<CompoundStructureRecord> generateElementRecords(SourceReference source) {
        List<CompoundStructureRecord> generated = new ArrayList<>();
        for (ElementRecord element : new BundledPeriodicTableRepository().findAll()) {
            generated.add(new CompoundStructureRecord(
                    "ELM_" + element.symbol().toUpperCase(Locale.ROOT),
                    normalizeFormula(Map.of(element.symbol(), 1)),
                    0,
                    element.name(),
                    "[" + element.symbol() + "]",
                    "",
                    "",
                    CompoundStructureKind.ELEMENTAL_ATOM,
                    source));
        }
        return List.copyOf(generated);
    }

    private static List<CompoundStructureRecord> generateElementMoleculeRecords(SourceReference source) {
        return List.of(
                discreteRecord("MOL_H2", Map.of("H", 2), 0, "Hydrogen", "[H][H]", source),
                discreteRecord("MOL_N2", Map.of("N", 2), 0, "Nitrogen", "N#N", source),
                discreteRecord("MOL_O2", Map.of("O", 2), 0, "Oxygen", "O=O", source),
                discreteRecord("MOL_F2", Map.of("F", 2), 0, "Fluorine", "FF", source),
                discreteRecord("MOL_CL2", Map.of("Cl", 2), 0, "Chlorine", "ClCl", source),
                discreteRecord("MOL_BR2", Map.of("Br", 2), 0, "Bromine", "BrBr", source),
                discreteRecord("MOL_I2", Map.of("I", 2), 0, "Iodine", "II", source));
    }

    private static List<CompoundStructureRecord> generateCuratedDiscreteRecords(SourceReference source) {
        List<CompoundStructureRecord> generated = new ArrayList<>();
        generated.add(discreteRecord("CMP_HCL", Map.of("H", 1, "Cl", 1), 0, "Hydrogen chloride", "[H]Cl", source));
        generated.add(discreteRecord("CMP_HBR", Map.of("H", 1, "Br", 1), 0, "Hydrogen bromide", "[H]Br", source));
        generated.add(discreteRecord("CMP_HI", Map.of("H", 1, "I", 1), 0, "Hydrogen iodide", "[H]I", source));
        generated.add(discreteRecord("CMP_HF", Map.of("H", 1, "F", 1), 0, "Hydrogen fluoride", "[H]F", source));
        generated.add(discreteRecord("CMP_H2S", Map.of("H", 2, "S", 1), 0, "Hydrogen sulfide", "S", source));
        generated.add(discreteRecord("CMP_H2O2", Map.of("H", 2, "O", 2), 0, "Hydrogen peroxide", "OO", source));
        generated.add(discreteRecord("CMP_NO", Map.of("N", 1, "O", 1), 0, "Nitric oxide", "[N]=O", source));
        generated.add(discreteRecord("CMP_NO2", Map.of("N", 1, "O", 2), 0, "Nitrogen dioxide", "O=N=O", source));
        generated.add(discreteRecord("CMP_N2O", Map.of("N", 2, "O", 1), 0, "Nitrous oxide", "N#N=O", source));
        generated.add(discreteRecord("CMP_SO2", Map.of("S", 1, "O", 2), 0, "Sulfur dioxide", "O=S=O", source));
        generated.add(discreteRecord("CMP_SO3", Map.of("S", 1, "O", 3), 0, "Sulfur trioxide", "O=S(=O)=O", source));
        generated.add(discreteRecord("CMP_H2CO3", Map.of("H", 2, "C", 1, "O", 3), 0, "Carbonic acid", "OC(=O)O", source));
        generated.add(discreteRecord("CMP_HNO3", Map.of("H", 1, "N", 1, "O", 3), 0, "Nitric acid", "O[N+](=O)[O-]", source));
        generated.add(discreteRecord("CMP_HNO2", Map.of("H", 1, "N", 1, "O", 2), 0, "Nitrous acid", "ON=O", source));
        generated.add(discreteRecord("CMP_H2SO4", Map.of("H", 2, "S", 1, "O", 4), 0, "Sulfuric acid", "OS(=O)(=O)O", source));
        generated.add(discreteRecord("CMP_H2SO3", Map.of("H", 2, "S", 1, "O", 3), 0, "Sulfurous acid", "OS(=O)O", source));
        generated.add(discreteRecord("CMP_H3PO4", Map.of("H", 3, "P", 1, "O", 4), 0, "Phosphoric acid", "OP(=O)(O)O", source));
        generated.add(discreteRecord("CMP_H3PO3", Map.of("H", 3, "P", 1, "O", 3), 0, "Phosphorous acid", "OP(O)O", source));
        generated.add(discreteRecord("CMP_HCN", Map.of("H", 1, "C", 1, "N", 1), 0, "Hydrogen cyanide", "C#N", source));
        generated.add(discreteRecord("CMP_CO", Map.of("C", 1, "O", 1), 0, "Carbon monoxide", "[C-]#[O+]", source));
        generated.add(discreteRecord("CMP_CS2", Map.of("C", 1, "S", 2), 0, "Carbon disulfide", "S=C=S", source));
        generated.add(discreteRecord("CMP_NH3", Map.of("N", 1, "H", 3), 0, "Ammonia", "N", source));
        generated.add(discreteRecord("CMP_PH3", Map.of("P", 1, "H", 3), 0, "Phosphine", "P", source));
        generated.add(discreteRecord("CMP_SIH4", Map.of("Si", 1, "H", 4), 0, "Silane", "[SiH4]", source));
        generated.add(discreteRecord("CMP_BH3", Map.of("B", 1, "H", 3), 0, "Borane", "[BH3]", source));
        generated.add(discreteRecord("CMP_C2H4", Map.of("C", 2, "H", 4), 0, "Ethene", "C=C", source));
        generated.add(discreteRecord("CMP_C2H2", Map.of("C", 2, "H", 2), 0, "Ethyne", "C#C", source));
        generated.add(discreteRecord("CMP_C3H8", Map.of("C", 3, "H", 8), 0, "Propane", "CCC", source));
        generated.add(discreteRecord("CMP_C4H10", Map.of("C", 4, "H", 10), 0, "Butane", "CCCC", source));
        generated.add(discreteRecord("CMP_C2H6O2", Map.of("C", 2, "H", 6, "O", 2), 0, "Ethylene glycol", "OCCO", source));
        generated.add(discreteRecord("CMP_CH4O", Map.of("C", 1, "H", 4, "O", 1), 0, "Methanol", "CO", source));
        generated.add(discreteRecord("CMP_C3H8O", Map.of("C", 3, "H", 8, "O", 1), 0, "2-Propanol", "CC(O)C", source));
        generated.add(discreteRecord("CMP_C3H6O", Map.of("C", 3, "H", 6, "O", 1), 0, "Acetone", "CC(=O)C", source));
        generated.add(discreteRecord("CMP_CH2O", Map.of("C", 1, "H", 2, "O", 1), 0, "Formaldehyde", "C=O", source));
        generated.add(discreteRecord("CMP_CH2O2", Map.of("C", 1, "H", 2, "O", 2), 0, "Formic acid", "OC=O", source));
        generated.add(discreteRecord("CMP_C2H4O", Map.of("C", 2, "H", 4, "O", 1), 0, "Acetaldehyde", "CC=O", source));
        generated.add(discreteRecord("CMP_C6H12O6", Map.of("C", 6, "H", 12, "O", 6), 0, "Glucose", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O", source));
        generated.add(discreteRecord("CMP_C12H22O11", Map.of("C", 12, "H", 22, "O", 11), 0, "Sucrose", "OC[C@H]1O[C@@](O)(CO)[C@H](O)[C@@H](O)[C@H]1O[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O", source));
        generated.add(discreteRecord("CMP_C6H6", Map.of("C", 6, "H", 6), 0, "Benzene", "c1ccccc1", source));
        generated.add(discreteRecord("CMP_C7H8", Map.of("C", 7, "H", 8), 0, "Toluene", "Cc1ccccc1", source));
        generated.add(discreteRecord("CMP_C8H10", Map.of("C", 8, "H", 10), 0, "Ethylbenzene", "CCc1ccccc1", source));
        generated.add(discreteRecord("CMP_C6H12", Map.of("C", 6, "H", 12), 0, "Cyclohexane", "C1CCCCC1", source));
        generated.add(discreteRecord("CMP_C4H8O2", Map.of("C", 4, "H", 8, "O", 2), 0, "Ethyl acetate", "CCOC(C)=O", source));
        generated.add(discreteRecord("CMP_C3H6O2", Map.of("C", 3, "H", 6, "O", 2), 0, "Methyl acetate", "COC(C)=O", source));
        generated.add(discreteRecord("CMP_C2H4O2", Map.of("C", 2, "H", 4, "O", 2), 0, "Acetic acid", "CC(=O)O", source));
        generated.add(discreteRecord("CMP_CHCl3", Map.of("C", 1, "H", 1, "Cl", 3), 0, "Chloroform", "ClC(Cl)Cl", source));
        generated.add(discreteRecord("CMP_C2H2Cl4", Map.of("C", 2, "H", 2, "Cl", 4), 0, "1,1,2,2-Tetrachloroethane", "ClC(Cl)C(Cl)Cl", source));
        generated.add(discreteRecord("CMP_C6H14", Map.of("C", 6, "H", 14), 0, "Hexane", "CCCCCC", source));
        generated.add(discreteRecord("CMP_C7H16", Map.of("C", 7, "H", 16), 0, "Heptane", "CCCCCCC", source));
        generated.add(discreteRecord("CMP_C8H18", Map.of("C", 8, "H", 18), 0, "Octane", "CCCCCCCC", source));
        return List.copyOf(generated);
    }

    private static List<CompoundStructureRecord> generateRepresentativeFormulaUnits(SourceReference source) {
        List<IonDefinition> cations = List.of(
                ion("LI", "Lithium", 1, "[Li+]", Map.of("Li", 1)),
                ion("NA", "Sodium", 1, "[Na+]", Map.of("Na", 1)),
                ion("K", "Potassium", 1, "[K+]", Map.of("K", 1)),
                ion("RB", "Rubidium", 1, "[Rb+]", Map.of("Rb", 1)),
                ion("CS", "Cesium", 1, "[Cs+]", Map.of("Cs", 1)),
                ion("AG", "Silver", 1, "[Ag+]", Map.of("Ag", 1)),
                ion("NH4", "Ammonium", 1, "[NH4+]", Map.of("N", 1, "H", 4)),
                ion("BE", "Beryllium", 2, "[Be+2]", Map.of("Be", 1)),
                ion("MG", "Magnesium", 2, "[Mg+2]", Map.of("Mg", 1)),
                ion("CA", "Calcium", 2, "[Ca+2]", Map.of("Ca", 1)),
                ion("SR", "Strontium", 2, "[Sr+2]", Map.of("Sr", 1)),
                ion("BA", "Barium", 2, "[Ba+2]", Map.of("Ba", 1)),
                ion("ZN", "Zinc", 2, "[Zn+2]", Map.of("Zn", 1)),
                ion("CD", "Cadmium", 2, "[Cd+2]", Map.of("Cd", 1)),
                ion("CU2", "Copper(II)", 2, "[Cu+2]", Map.of("Cu", 1)),
                ion("FE2", "Iron(II)", 2, "[Fe+2]", Map.of("Fe", 1)),
                ion("CO2", "Cobalt(II)", 2, "[Co+2]", Map.of("Co", 1)),
                ion("NI2", "Nickel(II)", 2, "[Ni+2]", Map.of("Ni", 1)),
                ion("MN2", "Manganese(II)", 2, "[Mn+2]", Map.of("Mn", 1)),
                ion("PB2", "Lead(II)", 2, "[Pb+2]", Map.of("Pb", 1)),
                ion("SN2", "Tin(II)", 2, "[Sn+2]", Map.of("Sn", 1)),
                ion("AL3", "Aluminum", 3, "[Al+3]", Map.of("Al", 1)),
                ion("FE3", "Iron(III)", 3, "[Fe+3]", Map.of("Fe", 1)),
                ion("CR3", "Chromium(III)", 3, "[Cr+3]", Map.of("Cr", 1)));

        List<IonDefinition> anions = List.of(
                ion("F", "fluoride", -1, "[F-]", Map.of("F", 1)),
                ion("CL", "chloride", -1, "[Cl-]", Map.of("Cl", 1)),
                ion("BR", "bromide", -1, "[Br-]", Map.of("Br", 1)),
                ion("I", "iodide", -1, "[I-]", Map.of("I", 1)),
                ion("OH", "hydroxide", -1, "[OH-]", Map.of("O", 1, "H", 1)),
                ion("NO3", "nitrate", -1, "[O-][N+](=O)[O-]", Map.of("N", 1, "O", 3)),
                ion("NO2", "nitrite", -1, "[O-]N=O", Map.of("N", 1, "O", 2)),
                ion("C2H3O2", "acetate", -1, "CC(=O)[O-]", Map.of("C", 2, "H", 3, "O", 2)),
                ion("CLO", "hypochlorite", -1, "[O-]Cl", Map.of("Cl", 1, "O", 1)),
                ion("CLO3", "chlorate", -1, "[O-]Cl(=O)=O", Map.of("Cl", 1, "O", 3)),
                ion("CLO4", "perchlorate", -1, "[O-]Cl(=O)(=O)=O", Map.of("Cl", 1, "O", 4)),
                ion("CN", "cyanide", -1, "[C-]#N", Map.of("C", 1, "N", 1)),
                ion("SCN", "thiocyanate", -1, "[S-]C#N", Map.of("S", 1, "C", 1, "N", 1)),
                ion("HCO3", "hydrogen carbonate", -1, "OC(=O)[O-]", Map.of("H", 1, "C", 1, "O", 3)),
                ion("HSO4", "hydrogen sulfate", -1, "OS(=O)(=O)[O-]", Map.of("H", 1, "S", 1, "O", 4)),
                ion("HSO3", "hydrogen sulfite", -1, "OS(=O)[O-]", Map.of("H", 1, "S", 1, "O", 3)),
                ion("O", "oxide", -2, "[O-2]", Map.of("O", 1)),
                ion("S", "sulfide", -2, "[S-2]", Map.of("S", 1)),
                ion("CO3", "carbonate", -2, "[O-]C(=O)[O-]", Map.of("C", 1, "O", 3)),
                ion("SO3", "sulfite", -2, "[O-]S(=O)[O-]", Map.of("S", 1, "O", 3)),
                ion("SO4", "sulfate", -2, "[O-]S(=O)(=O)[O-]", Map.of("S", 1, "O", 4)),
                ion("C2O4", "oxalate", -2, "[O-]C(=O)C(=O)[O-]", Map.of("C", 2, "O", 4)),
                ion("CRO4", "chromate", -2, "[O-][Cr](=O)(=O)[O-]", Map.of("Cr", 1, "O", 4)),
                ion("CR2O7", "dichromate", -2, "[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-]", Map.of("Cr", 2, "O", 7)),
                ion("HPO4", "hydrogen phosphate", -2, "OP(=O)([O-])[O-]", Map.of("H", 1, "P", 1, "O", 4)),
                ion("PO4", "phosphate", -3, "[O-]P([O-])([O-])=O", Map.of("P", 1, "O", 4)),
                ion("BO3", "borate", -3, "[O-]B([O-])[O-]", Map.of("B", 1, "O", 3)));

        List<CompoundStructureRecord> generated = new ArrayList<>();
        for (IonDefinition cation : cations) {
            for (IonDefinition anion : anions) {
                generated.add(buildRepresentativeFormulaUnit(cation, anion, source));
            }
        }
        return List.copyOf(generated);
    }

    private static CompoundStructureRecord buildRepresentativeFormulaUnit(
            IonDefinition cation,
            IonDefinition anion,
            SourceReference source) {
        int cationCharge = cation.charge();
        int anionCharge = Math.abs(anion.charge());
        int ratio = leastCommonMultiple(cationCharge, anionCharge);
        int cationCount = ratio / cationCharge;
        int anionCount = ratio / anionCharge;

        Map<String, Integer> composition = new LinkedHashMap<>();
        mergeComposition(composition, cation.composition(), cationCount);
        mergeComposition(composition, anion.composition(), anionCount);

        StringBuilder smiles = new StringBuilder();
        appendRepeatedSmiles(smiles, cation.smiles(), cationCount);
        appendRepeatedSmiles(smiles, anion.smiles(), anionCount);

        return new CompoundStructureRecord(
                "IONIC_" + cation.id() + "_" + anion.id(),
                normalizeFormula(composition),
                0,
                cation.name() + " " + anion.name(),
                smiles.toString(),
                "",
                "",
                CompoundStructureKind.REPRESENTATIVE_FORMULA_UNIT,
                source);
    }

    private static void appendRepeatedSmiles(StringBuilder builder, String smiles, int count) {
        for (int index = 0; index < count; index++) {
            if (builder.length() > 0) {
                builder.append('.');
            }
            builder.append(smiles);
        }
    }

    private static CompoundStructureRecord discreteRecord(
            String id,
            Map<String, Integer> composition,
            int charge,
            String preferredName,
            String canonicalSmiles,
            SourceReference source) {
        CompoundStructureKind kind = isElementalMolecule(composition)
                ? CompoundStructureKind.ELEMENTAL_MOLECULE
                : CompoundStructureKind.DISCRETE_MOLECULE;
        return new CompoundStructureRecord(
                id,
                normalizeFormula(composition),
                charge,
                preferredName,
                canonicalSmiles,
                "",
                "",
                kind,
                source);
    }

    private static boolean isElementalMolecule(Map<String, Integer> composition) {
        return composition.size() == 1 && !composition.values().contains(1);
    }

    private static IonDefinition ion(String id, String name, int charge, String smiles, Map<String, Integer> composition) {
        return new IonDefinition(id, name, charge, smiles, Map.copyOf(composition));
    }

    private static void mergeComposition(Map<String, Integer> target, Map<String, Integer> source, int multiplier) {
        source.forEach((symbol, count) -> target.merge(symbol, count * multiplier, Integer::sum));
    }

    private static String normalizeFormula(Map<String, Integer> composition) {
        Map<String, Integer> ordered = new LinkedHashMap<>();
        if (composition.containsKey("C")) {
            ordered.put("C", composition.get("C"));
            if (composition.containsKey("H")) {
                ordered.put("H", composition.get("H"));
            }
        }
        composition.keySet().stream()
                .filter(symbol -> !ordered.containsKey(symbol))
                .sorted()
                .forEach(symbol -> ordered.put(symbol, composition.get(symbol)));

        StringBuilder builder = new StringBuilder();
        ordered.forEach((symbol, count) -> {
            builder.append(symbol);
            if (count != 1) {
                builder.append(count);
            }
        });
        return builder.toString();
    }

    private static CompoundStructureKind parseStructureKind(String value) {
        if (value == null || value.isBlank()) {
            return CompoundStructureKind.DISCRETE_MOLECULE;
        }
        return CompoundStructureKind.valueOf(value.trim().toUpperCase(Locale.ROOT));
    }

    private static String fallbackSeedId(String formula, int charge, String preferredName) {
        return ("SEED_" + formula + "_" + charge + "_" + preferredName)
                .replaceAll("[^A-Za-z0-9]+", "_")
                .replaceAll("_+", "_")
                .replaceAll("^_|_$", "")
                .toUpperCase(Locale.ROOT);
    }

    private static int leastCommonMultiple(int left, int right) {
        return left / greatestCommonDivisor(left, right) * right;
    }

    private static int greatestCommonDivisor(int left, int right) {
        int a = Math.abs(left);
        int b = Math.abs(right);
        while (b != 0) {
            int next = a % b;
            a = b;
            b = next;
        }
        return a;
    }

    private record IonDefinition(String id, String name, int charge, String smiles, Map<String, Integer> composition) {
    }
}
