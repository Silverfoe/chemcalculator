package com.chemcalc.domain;

import java.util.Objects;

public record BalanceRowOperation(String label, String detail, BalanceMatrixSnapshot snapshot) {
    public BalanceRowOperation {
        label = Objects.requireNonNullElse(label, "").trim();
        detail = Objects.requireNonNullElse(detail, "").trim();
        snapshot = snapshot == null ? BalanceMatrixSnapshot.empty() : snapshot;
    }
}
