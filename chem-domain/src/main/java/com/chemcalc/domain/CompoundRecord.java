package com.chemcalc.domain;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;

public record CompoundRecord(
        String formula,
        String preferredName,
        Map<String, PropertyValue<String>> properties,
        SourceReference source) {

    public CompoundRecord {
        formula = Objects.requireNonNullElse(formula, "").trim();
        preferredName = Objects.requireNonNullElse(preferredName, "").trim();
        properties = Map.copyOf(new LinkedHashMap<>(properties == null ? Map.of() : properties));
        source = source == null ? SourceReference.unknown() : source;
    }
}
