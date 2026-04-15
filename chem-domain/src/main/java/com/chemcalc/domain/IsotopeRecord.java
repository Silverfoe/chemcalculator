package com.chemcalc.domain;

import java.util.Objects;

public record IsotopeRecord(
        String symbol,
        String massNumber,
        PropertyValue<String> atomicMass,
        PropertyValue<String> abundance,
        String notes) {

    public IsotopeRecord {
        symbol = Objects.requireNonNullElse(symbol, "").trim();
        massNumber = Objects.requireNonNullElse(massNumber, "").trim();
        atomicMass = atomicMass == null ? PropertyValue.of("", "", PropertyKind.UNKNOWN, SourceReference.unknown()) : atomicMass;
        abundance = abundance == null ? PropertyValue.of("", "", PropertyKind.UNKNOWN, SourceReference.unknown()) : abundance;
        notes = Objects.requireNonNullElse(notes, "").trim();
    }
}
