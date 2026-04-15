package com.chemcalc.domain;

import java.util.Objects;

public record BalanceConservationCheck(String label, int reactantTotal, int productTotal) {
    public BalanceConservationCheck {
        label = Objects.requireNonNullElse(label, "").trim();
    }
}
