package com.chemcalc.domain;

import java.util.Objects;

public record StructureCandidate(
        String id,
        String preferredName,
        String canonicalFormula,
        int charge,
        String canonicalSmiles,
        SourceReference source) {

    public StructureCandidate {
        id = Objects.requireNonNullElse(id, "").trim();
        preferredName = Objects.requireNonNullElse(preferredName, "").trim();
        canonicalFormula = Objects.requireNonNullElse(canonicalFormula, "").trim();
        canonicalSmiles = Objects.requireNonNullElse(canonicalSmiles, "").trim();
        source = source == null ? SourceReference.unknown() : source;
    }
}
