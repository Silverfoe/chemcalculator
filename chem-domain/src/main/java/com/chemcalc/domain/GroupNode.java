package com.chemcalc.domain;

import java.util.List;

public record GroupNode(char openingDelimiter, char closingDelimiter, List<FormulaPart> parts, int multiplier)
        implements FormulaPart {

    public GroupNode {
        parts = List.copyOf(parts == null ? List.of() : parts);
        if (multiplier < 1) {
            throw new IllegalArgumentException("Group multiplier must be >= 1");
        }
    }

    @Override
    public String toFormulaString() {
        StringBuilder builder = new StringBuilder();
        builder.append(openingDelimiter);
        for (FormulaPart part : parts) {
            builder.append(part.toFormulaString());
        }
        builder.append(closingDelimiter);
        if (multiplier != 1) {
            builder.append(multiplier);
        }
        return builder.toString();
    }
}
