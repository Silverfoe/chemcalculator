package com.chemcalc.domain;

import java.util.Objects;

public record StructureRequest(
        StructureInputFormat format,
        String input,
        String selectedCandidateKey) {

    public StructureRequest {
        format = format == null ? StructureInputFormat.FORMULA : format;
        input = Objects.requireNonNullElse(input, "").trim();
        selectedCandidateKey = Objects.requireNonNullElse(selectedCandidateKey, "").trim();
    }

    public StructureRequest(StructureInputFormat format, String input) {
        this(format, input, "");
    }
}
