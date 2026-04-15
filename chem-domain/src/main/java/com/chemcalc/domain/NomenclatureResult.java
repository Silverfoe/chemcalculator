package com.chemcalc.domain;

import java.util.List;
import java.util.Objects;

public record NomenclatureResult(
        String displayName,
        boolean supported,
        List<ExplanationStep> explanation,
        List<ValidationMessage> messages) {

    public NomenclatureResult {
        displayName = Objects.requireNonNullElse(displayName, "").trim();
        explanation = List.copyOf(explanation == null ? List.of() : explanation);
        messages = List.copyOf(messages == null ? List.of() : messages);
    }

    public static NomenclatureResult unsupported(String message) {
        return new NomenclatureResult(
                "",
                false,
                List.of(),
                List.of(ValidationMessage.warning("UNSUPPORTED_NOMENCLATURE", message)));
    }
}
