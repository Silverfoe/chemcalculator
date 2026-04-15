package com.chemcalc.domain;

import java.util.List;
import java.util.Objects;

public record BalanceNormalizationTrace(
        String leastCommonMultiple,
        List<String> scaledVector,
        boolean signFlipped,
        List<String> signNormalizedVector,
        String greatestCommonDivisor,
        List<Integer> coefficients) {

    public static BalanceNormalizationTrace empty() {
        return new BalanceNormalizationTrace("", List.of(), false, List.of(), "", List.of());
    }

    public BalanceNormalizationTrace {
        leastCommonMultiple = Objects.requireNonNullElse(leastCommonMultiple, "").trim();
        scaledVector = List.copyOf(scaledVector == null ? List.of() : scaledVector);
        signNormalizedVector = List.copyOf(signNormalizedVector == null ? List.of() : signNormalizedVector);
        greatestCommonDivisor = Objects.requireNonNullElse(greatestCommonDivisor, "").trim();
        coefficients = List.copyOf(coefficients == null ? List.of() : coefficients);
    }
}
