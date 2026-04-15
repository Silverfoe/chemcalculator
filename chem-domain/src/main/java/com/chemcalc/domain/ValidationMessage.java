package com.chemcalc.domain;

import java.util.Objects;

public record ValidationMessage(ValidationSeverity severity, String code, String message) {
    public ValidationMessage {
        severity = severity == null ? ValidationSeverity.INFO : severity;
        code = Objects.requireNonNullElse(code, "").trim();
        message = Objects.requireNonNullElse(message, "").trim();
    }

    public static ValidationMessage info(String code, String message) {
        return new ValidationMessage(ValidationSeverity.INFO, code, message);
    }

    public static ValidationMessage warning(String code, String message) {
        return new ValidationMessage(ValidationSeverity.WARNING, code, message);
    }

    public static ValidationMessage error(String code, String message) {
        return new ValidationMessage(ValidationSeverity.ERROR, code, message);
    }
}
