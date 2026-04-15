package com.chemcalc.domain;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;

public record ThermoRecord(String subject, Map<String, PropertyValue<String>> properties, SourceReference source) {
    public ThermoRecord {
        subject = Objects.requireNonNullElse(subject, "").trim();
        properties = Map.copyOf(new LinkedHashMap<>(properties == null ? Map.of() : properties));
        source = source == null ? SourceReference.unknown() : source;
    }
}
