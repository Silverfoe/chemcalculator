package com.chemcalc.data;

import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.PeriodicTableRepository;
import com.chemcalc.domain.PropertyKind;
import com.chemcalc.domain.PropertyValue;
import com.chemcalc.domain.SourceReference;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

public final class SqlitePeriodicTableRepository implements PeriodicTableRepository {
    private final String jdbcUrl;

    public SqlitePeriodicTableRepository(String jdbcUrl) {
        this.jdbcUrl = jdbcUrl;
    }

    @Override
    public List<ElementRecord> findAll() {
        String sql = """
                select atomic_number, symbol, name, group_number, period_number, block, category,
                       standard_state, year_discovered, atomic_weight, atomic_weight_source
                from elements
                order by atomic_number
                """;
        try (Connection connection = DriverManager.getConnection(jdbcUrl);
             PreparedStatement statement = connection.prepareStatement(sql);
             ResultSet resultSet = statement.executeQuery()) {
            List<ElementRecord> elements = new ArrayList<>();
            while (resultSet.next()) {
                Map<String, PropertyValue<String>> properties = new LinkedHashMap<>();
                String atomicWeight = resultSet.getString("atomic_weight");
                String sourceId = resultSet.getString("atomic_weight_source");
                if (atomicWeight != null && !atomicWeight.isBlank()) {
                    properties.put("standardAtomicWeight",
                            PropertyValue.of(atomicWeight, "g/mol", PropertyKind.TEXT, new SourceReference(sourceId, sourceId, "", "")));
                }
                elements.add(new ElementRecord(
                        resultSet.getInt("atomic_number"),
                        resultSet.getString("symbol"),
                        resultSet.getString("name"),
                        resultSet.getObject("group_number", Integer.class),
                        resultSet.getInt("period_number"),
                        resultSet.getString("block"),
                        resultSet.getString("category"),
                        resultSet.getString("standard_state"),
                        resultSet.getString("year_discovered"),
                        properties,
                        List.of()));
            }
            return List.copyOf(elements);
        } catch (SQLException exception) {
            throw new IllegalStateException("Failed to query SQLite periodic table", exception);
        }
    }

    @Override
    public Optional<ElementRecord> findByAtomicNumber(int atomicNumber) {
        return findAll().stream().filter(element -> element.atomicNumber() == atomicNumber).findFirst();
    }

    @Override
    public Optional<ElementRecord> findBySymbol(String symbol) {
        if (symbol == null || symbol.isBlank()) {
            return Optional.empty();
        }
        String normalized = symbol.trim().toUpperCase(Locale.ROOT);
        return findAll().stream()
                .filter(element -> element.symbol().toUpperCase(Locale.ROOT).equals(normalized))
                .findFirst();
    }

    @Override
    public List<ElementRecord> search(String query) {
        if (query == null || query.isBlank()) {
            return findAll();
        }
        String normalized = query.trim().toLowerCase(Locale.ROOT);
        return findAll().stream()
                .filter(element -> element.symbol().toLowerCase(Locale.ROOT).contains(normalized)
                        || element.name().toLowerCase(Locale.ROOT).contains(normalized))
                .toList();
    }
}
