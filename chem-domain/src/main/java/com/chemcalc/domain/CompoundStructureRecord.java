package com.chemcalc.domain;

import java.util.Objects;

public record CompoundStructureRecord(
        String id,
        String canonicalFormula,
        int charge,
        String preferredName,
        String canonicalSmiles,
        String inchi,
        String inchiKey,
        CompoundStructureKind structureKind,
        SourceReference source) {

    public CompoundStructureRecord {
        id = Objects.requireNonNullElse(id, "").trim();
        canonicalFormula = Objects.requireNonNullElse(canonicalFormula, "").trim();
        preferredName = Objects.requireNonNullElse(preferredName, "").trim();
        canonicalSmiles = Objects.requireNonNullElse(canonicalSmiles, "").trim();
        inchi = Objects.requireNonNullElse(inchi, "").trim();
        inchiKey = Objects.requireNonNullElse(inchiKey, "").trim();
        structureKind = structureKind == null ? CompoundStructureKind.DISCRETE_MOLECULE : structureKind;
        source = source == null ? SourceReference.unknown() : source;
    }
}
