package com.chemcalc.domain;

import java.util.Objects;

public record ExplanationStep(String title, String detail) {
    public ExplanationStep {
        title = Objects.requireNonNullElse(title, "").trim();
        detail = Objects.requireNonNullElse(detail, "").trim();
    }
}
