package com.chemcalc.domain;

import java.util.List;

public record FormulaSegment(int multiplier, List<FormulaPart> parts) {
    public FormulaSegment {
        if (multiplier < 1) {
            throw new IllegalArgumentException("Segment multiplier must be >= 1");
        }
        parts = List.copyOf(parts == null ? List.of() : parts);
    }

    public String toFormulaString() {
        StringBuilder builder = new StringBuilder();
        if (multiplier != 1) {
            builder.append(multiplier);
        }
        for (FormulaPart part : parts) {
            builder.append(part.toFormulaString());
        }
        return builder.toString();
    }
}
