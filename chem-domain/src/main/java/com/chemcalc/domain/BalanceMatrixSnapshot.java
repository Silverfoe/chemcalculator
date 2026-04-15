package com.chemcalc.domain;

import java.util.List;

public record BalanceMatrixSnapshot(
        List<String> rowLabels,
        List<String> columnLabels,
        List<List<String>> values) {

    public static BalanceMatrixSnapshot empty() {
        return new BalanceMatrixSnapshot(List.of(), List.of(), List.of());
    }

    public BalanceMatrixSnapshot {
        rowLabels = List.copyOf(rowLabels == null ? List.of() : rowLabels);
        columnLabels = List.copyOf(columnLabels == null ? List.of() : columnLabels);
        values = List.copyOf((values == null ? List.<List<String>>of() : values).stream()
                .map(row -> List.copyOf(row == null ? List.<String>of() : row))
                .toList());
    }
}
