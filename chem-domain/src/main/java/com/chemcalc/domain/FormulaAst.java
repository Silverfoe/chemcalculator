package com.chemcalc.domain;

import java.util.List;

public record FormulaAst(List<FormulaSegment> segments) {
    public FormulaAst {
        segments = List.copyOf(segments == null ? List.of() : segments);
    }

    public String toCanonicalFormula() {
        return segments.stream()
                .map(FormulaSegment::toFormulaString)
                .reduce((left, right) -> left + "·" + right)
                .orElse("");
    }
}
