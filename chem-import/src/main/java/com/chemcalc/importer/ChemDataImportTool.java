package com.chemcalc.importer;

import com.chemcalc.data.BundledCompoundStructureRepository;
import com.chemcalc.domain.CompoundStructureRecord;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;

public final class ChemDataImportTool {
    private ChemDataImportTool() {
    }

    public static void main(String[] args) throws Exception {
        String jdbcUrl = args.length > 0 ? args[0] : "jdbc:sqlite:build/generated/chem-data.db";
        try (Connection connection = DriverManager.getConnection(jdbcUrl)) {
            connection.setAutoCommit(false);
            applySchema(connection);
            importSources(connection);
            importElements(connection);
            importCompoundStructures(connection);
            connection.commit();
        }
        System.out.println("Imported bundled chemistry snapshot into " + jdbcUrl);
    }

    private static void applySchema(Connection connection) throws IOException, SQLException {
        String schema = readResource("/db/schema.sql");
        try (Statement statement = connection.createStatement()) {
            for (String command : schema.split(";")) {
                String trimmed = command.trim();
                if (!trimmed.isBlank()) {
                    statement.execute(trimmed);
                }
            }
        }
    }

    private static void importSources(Connection connection) throws IOException, SQLException {
        try (Statement statement = connection.createStatement()) {
            statement.executeUpdate("delete from sources");
        }
        try (PreparedStatement statement = connection.prepareStatement(
                "insert into sources(id, title, revision, url) values (?, ?, ?, ?)")) {
            for (String[] row : readRows("/data/sources.psv")) {
                statement.setString(1, row[0]);
                statement.setString(2, row[1]);
                statement.setString(3, row[2]);
                statement.setString(4, row[3]);
                statement.addBatch();
            }
            statement.executeBatch();
        }
    }

    private static void importElements(Connection connection) throws IOException, SQLException {
        try (Statement statement = connection.createStatement()) {
            statement.executeUpdate("delete from elements");
        }
        try (PreparedStatement statement = connection.prepareStatement("""
                insert into elements(
                    atomic_number, symbol, name, group_number, period_number, block, category,
                    standard_state, year_discovered, atomic_weight, atomic_weight_source
                ) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """)) {
            for (String[] row : readRows("/data/elements.psv")) {
                int atomicNumber = Integer.parseInt(row[0]);
                String symbol = row[1];
                String name = row[2];
                String atomicWeight = row[3];
                String atomicWeightSource = row[4];
                String standardState = row[5];
                String yearDiscovered = row[6];

                PeriodicTopology layout = PeriodicTopology.forAtomicNumber(atomicNumber, symbol);

                statement.setInt(1, atomicNumber);
                statement.setString(2, symbol);
                statement.setString(3, name);
                if (layout.groupNumber() == null) {
                    statement.setObject(4, null);
                } else {
                    statement.setInt(4, layout.groupNumber());
                }
                statement.setInt(5, layout.periodNumber());
                statement.setString(6, layout.block());
                statement.setString(7, layout.category());
                statement.setString(8, standardState);
                statement.setString(9, yearDiscovered);
                statement.setString(10, atomicWeight);
                statement.setString(11, atomicWeightSource);
                statement.addBatch();
            }
            statement.executeBatch();
        }
    }

    private static void importCompoundStructures(Connection connection) throws IOException, SQLException {
        try (Statement statement = connection.createStatement()) {
            statement.executeUpdate("delete from compound_structures");
        }
        try (PreparedStatement statement = connection.prepareStatement("""
                insert into compound_structures(
                    id, canonical_formula, charge, preferred_name, canonical_smiles, inchi, inchi_key, structure_kind, source_id
                ) values (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """)) {
            for (CompoundStructureRecord record : new BundledCompoundStructureRepository().findAll()) {
                statement.setString(1, record.id());
                statement.setString(2, record.canonicalFormula());
                statement.setInt(3, record.charge());
                statement.setString(4, record.preferredName());
                statement.setString(5, record.canonicalSmiles());
                statement.setString(6, record.inchi());
                statement.setString(7, record.inchiKey());
                statement.setString(8, record.structureKind().name());
                statement.setString(9, record.source().id());
                statement.addBatch();
            }
            statement.executeBatch();
        }
    }

    private static java.util.List<String[]> readRows(String resourcePath) throws IOException {
        try (InputStream input = ChemDataImportTool.class.getResourceAsStream(resourcePath)) {
            if (input == null) {
                throw new IOException("Missing resource " + resourcePath);
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8))) {
                java.util.List<String[]> rows = new java.util.ArrayList<>();
                String line = reader.readLine();
                while ((line = reader.readLine()) != null) {
                    if (line.isBlank() || line.startsWith("#")) {
                        continue;
                    }
                    rows.add(line.split("\\|", -1));
                }
                return rows;
            }
        }
    }

    private static String readResource(String resourcePath) throws IOException {
        try (InputStream input = ChemDataImportTool.class.getResourceAsStream(resourcePath)) {
            if (input == null) {
                throw new IOException("Missing resource " + resourcePath);
            }
            return new String(input.readAllBytes(), StandardCharsets.UTF_8);
        }
    }

    private record PeriodicTopology(Integer groupNumber, int periodNumber, String block, String category) {
        private static PeriodicTopology forAtomicNumber(int atomicNumber, String symbol) {
            Integer group;
            int period;

            if (atomicNumber == 1) {
                group = 1;
                period = 1;
            } else if (atomicNumber == 2) {
                group = 18;
                period = 1;
            } else if (atomicNumber <= 10) {
                period = 2;
                group = atomicNumber <= 4 ? atomicNumber - 2 : atomicNumber + 8;
            } else if (atomicNumber <= 18) {
                period = 3;
                group = atomicNumber <= 12 ? atomicNumber - 10 : atomicNumber;
            } else if (atomicNumber <= 36) {
                period = 4;
                group = atomicNumber - 18;
            } else if (atomicNumber <= 54) {
                period = 5;
                group = atomicNumber - 36;
            } else if (atomicNumber <= 56) {
                period = 6;
                group = atomicNumber - 54;
            } else if (atomicNumber <= 71) {
                period = 6;
                group = 3;
            } else if (atomicNumber <= 86) {
                period = 6;
                group = atomicNumber - 68;
            } else if (atomicNumber <= 88) {
                period = 7;
                group = atomicNumber - 86;
            } else if (atomicNumber <= 103) {
                period = 7;
                group = 3;
            } else {
                period = 7;
                group = atomicNumber - 100;
            }

            String block;
            if (atomicNumber >= 57 && atomicNumber <= 71 || atomicNumber >= 89 && atomicNumber <= 103) {
                block = "f";
            } else if ("He".equals(symbol) || group != null && group <= 2) {
                block = "s";
            } else if (group != null && group <= 12) {
                block = "d";
            } else {
                block = "p";
            }

            String category;
            if (atomicNumber >= 57 && atomicNumber <= 71) {
                category = "Lanthanide";
            } else if (atomicNumber >= 89 && atomicNumber <= 103) {
                category = "Actinide";
            } else if (java.util.Set.of("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og").contains(symbol)) {
                category = "Noble Gas";
            } else if (java.util.Set.of("F", "Cl", "Br", "I", "At", "Ts").contains(symbol)) {
                category = "Halogen";
            } else if (java.util.Set.of("Li", "Na", "K", "Rb", "Cs", "Fr").contains(symbol)) {
                category = "Alkali Metal";
            } else if (java.util.Set.of("Be", "Mg", "Ca", "Sr", "Ba", "Ra").contains(symbol)) {
                category = "Alkaline Earth Metal";
            } else if (java.util.Set.of("B", "Si", "Ge", "As", "Sb", "Te", "Po").contains(symbol)) {
                category = "Metalloid";
            } else if (java.util.Set.of("Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi", "Nh", "Fl", "Mc", "Lv").contains(symbol)) {
                category = "Post-Transition Metal";
            } else if (java.util.Set.of("H", "C", "N", "O", "P", "S", "Se").contains(symbol)) {
                category = "Nonmetal";
            } else if (group != null && group >= 3 && group <= 12) {
                category = "Transition Metal";
            } else {
                category = "Metal";
            }

            return new PeriodicTopology(group, period, block, category);
        }
    }
}
