package com.chemcalc.domain;

import java.util.List;

public record BalanceTrace(
        List<BalanceConservedQuantity> quantities,
        BalanceMatrixSnapshot initialMatrix,
        List<BalanceRowOperation> rowOperations,
        BalanceSolutionTrace solution,
        BalanceNormalizationTrace normalization,
        List<BalanceConservationCheck> conservationChecks) {

    public static BalanceTrace empty() {
        return new BalanceTrace(
                List.of(),
                BalanceMatrixSnapshot.empty(),
                List.of(),
                BalanceSolutionTrace.empty(),
                BalanceNormalizationTrace.empty(),
                List.of());
    }

    public BalanceTrace {
        quantities = List.copyOf(quantities == null ? List.of() : quantities);
        initialMatrix = initialMatrix == null ? BalanceMatrixSnapshot.empty() : initialMatrix;
        rowOperations = List.copyOf(rowOperations == null ? List.of() : rowOperations);
        solution = solution == null ? BalanceSolutionTrace.empty() : solution;
        normalization = normalization == null ? BalanceNormalizationTrace.empty() : normalization;
        conservationChecks = List.copyOf(conservationChecks == null ? List.of() : conservationChecks);
    }
}
