package com.chemcalc.domain;

import java.util.Objects;

public record IonRecord(String formula, String name, int charge, PropertyValue<String> notes) {
    public IonRecord {
        formula = Objects.requireNonNullElse(formula, "").trim();
        name = Objects.requireNonNullElse(name, "").trim();
        notes = notes == null ? PropertyValue.of("", "", PropertyKind.UNKNOWN, SourceReference.unknown()) : notes;
    }
}
