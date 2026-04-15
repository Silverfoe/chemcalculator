package com.chemcalc.domain;

import java.util.List;
import java.util.Objects;

public record BalanceSolutionTrace(
        List<String> pivotVariables,
        List<String> freeVariables,
        List<String> nullspaceVector,
        String detail) {

    public static BalanceSolutionTrace empty() {
        return new BalanceSolutionTrace(List.of(), List.of(), List.of(), "");
    }

    public BalanceSolutionTrace {
        pivotVariables = List.copyOf(pivotVariables == null ? List.of() : pivotVariables);
        freeVariables = List.copyOf(freeVariables == null ? List.of() : freeVariables);
        nullspaceVector = List.copyOf(nullspaceVector == null ? List.of() : nullspaceVector);
        detail = Objects.requireNonNullElse(detail, "").trim();
    }
}
