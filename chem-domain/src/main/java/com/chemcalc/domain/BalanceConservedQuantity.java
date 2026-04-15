package com.chemcalc.domain;

import java.util.Objects;

public record BalanceConservedQuantity(String label, String kind) {
    public BalanceConservedQuantity {
        label = Objects.requireNonNullElse(label, "").trim();
        kind = Objects.requireNonNullElse(kind, "").trim();
    }
}
