package com.chemcalc.domain;

import java.util.Objects;

public record ElementNode(String symbol, int count) implements FormulaPart {
    public ElementNode {
        symbol = Objects.requireNonNullElse(symbol, "").trim();
        if (count < 1) {
            throw new IllegalArgumentException("Element count must be >= 1");
        }
    }

    @Override
    public String toFormulaString() {
        return symbol + (count == 1 ? "" : Integer.toString(count));
    }
}
