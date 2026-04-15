package com.chemcalc.data;

import com.chemcalc.domain.CompoundStructureRecord;
import com.chemcalc.domain.CompoundStructureRepository;
import com.chemcalc.domain.CompoundStructureKind;
import com.chemcalc.domain.SourceReference;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public final class SqliteCompoundStructureRepository implements CompoundStructureRepository {
    private final String jdbcUrl;

    public SqliteCompoundStructureRepository(String jdbcUrl) {
        this.jdbcUrl = jdbcUrl;
    }

    @Override
    public List<CompoundStructureRecord> findAll() {
        String sql = """
                select id, canonical_formula, charge, preferred_name, canonical_smiles, inchi, inchi_key, structure_kind, source_id
                from compound_structures
                order by canonical_formula, charge, preferred_name
                """;
        try (Connection connection = DriverManager.getConnection(jdbcUrl);
             PreparedStatement statement = connection.prepareStatement(sql);
             ResultSet resultSet = statement.executeQuery()) {
            List<CompoundStructureRecord> records = new ArrayList<>();
            while (resultSet.next()) {
                String sourceId = resultSet.getString("source_id");
                records.add(new CompoundStructureRecord(
                        resultSet.getString("id"),
                        resultSet.getString("canonical_formula"),
                        resultSet.getInt("charge"),
                        resultSet.getString("preferred_name"),
                        resultSet.getString("canonical_smiles"),
                        resultSet.getString("inchi"),
                        resultSet.getString("inchi_key"),
                        CompoundStructureKind.valueOf(resultSet.getString("structure_kind")),
                        new SourceReference(sourceId, sourceId, "", "")));
            }
            return List.copyOf(records);
        } catch (SQLException exception) {
            throw new IllegalStateException("Failed to query SQLite compound structures", exception);
        }
    }

    @Override
    public List<CompoundStructureRecord> findByFormulaAndCharge(String canonicalFormula, int charge) {
        if (canonicalFormula == null || canonicalFormula.isBlank()) {
            return List.of();
        }
        String sql = """
                select id, canonical_formula, charge, preferred_name, canonical_smiles, inchi, inchi_key, structure_kind, source_id
                from compound_structures
                where upper(canonical_formula) = ? and charge = ?
                order by preferred_name
                """;
        try (Connection connection = DriverManager.getConnection(jdbcUrl);
             PreparedStatement statement = connection.prepareStatement(sql)) {
            statement.setString(1, canonicalFormula.trim().toUpperCase(Locale.ROOT));
            statement.setInt(2, charge);
            try (ResultSet resultSet = statement.executeQuery()) {
                List<CompoundStructureRecord> records = new ArrayList<>();
                while (resultSet.next()) {
                    String sourceId = resultSet.getString("source_id");
                    records.add(new CompoundStructureRecord(
                            resultSet.getString("id"),
                            resultSet.getString("canonical_formula"),
                            resultSet.getInt("charge"),
                            resultSet.getString("preferred_name"),
                            resultSet.getString("canonical_smiles"),
                            resultSet.getString("inchi"),
                            resultSet.getString("inchi_key"),
                            CompoundStructureKind.valueOf(resultSet.getString("structure_kind")),
                            new SourceReference(sourceId, sourceId, "", "")));
                }
                return List.copyOf(records);
            }
        } catch (SQLException exception) {
            throw new IllegalStateException("Failed to query SQLite compound structures", exception);
        }
    }
}
