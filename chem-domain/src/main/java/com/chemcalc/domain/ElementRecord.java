package com.chemcalc.domain;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;

public record ElementRecord(
        int atomicNumber,
        String symbol,
        String name,
        Integer group,
        int period,
        String block,
        String category,
        String standardState,
        String yearDiscovered,
        Map<String, PropertyValue<String>> properties,
        List<IsotopeRecord> isotopes) {

    public ElementRecord {
        symbol = Objects.requireNonNullElse(symbol, "").trim();
        name = Objects.requireNonNullElse(name, "").trim();
        block = Objects.requireNonNullElse(block, "").trim();
        category = Objects.requireNonNullElse(category, "").trim();
        standardState = Objects.requireNonNullElse(standardState, "").trim();
        yearDiscovered = Objects.requireNonNullElse(yearDiscovered, "").trim();
        properties = Map.copyOf(new LinkedHashMap<>(properties == null ? Map.of() : properties));
        isotopes = List.copyOf(isotopes == null ? List.of() : isotopes);
    }

    public Optional<PropertyValue<String>> property(String key) {
        return Optional.ofNullable(properties.get(key));
    }
}
