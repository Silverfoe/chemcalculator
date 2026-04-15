package com.chemcalc.domain;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Objects;
import java.util.Collections;

public record Species(
        int coefficient,
        String originalInput,
        FormulaAst formulaAst,
        Map<String, Integer> composition,
        int charge,
        Phase phase) {

    public Species {
        if (coefficient < 1) {
            throw new IllegalArgumentException("Species coefficient must be >= 1");
        }
        originalInput = Objects.requireNonNullElse(originalInput, "").trim();
        formulaAst = formulaAst == null ? new FormulaAst(java.util.List.of()) : formulaAst;
        composition = Collections.unmodifiableMap(new LinkedHashMap<>(composition == null ? Map.of() : composition));
        phase = phase == null ? Phase.UNKNOWN : phase;
    }

    public String canonicalFormula() {
        return formulaAst.toCanonicalFormula();
    }

    public String displayFormula() {
        StringBuilder builder = new StringBuilder(canonicalFormula());
        if (charge != 0) {
            int magnitude = Math.abs(charge);
            if (magnitude != 1 && !canonicalFormula().isEmpty() && Character.isDigit(canonicalFormula().charAt(canonicalFormula().length() - 1))) {
                builder.append("^").append(magnitude);
            } else if (magnitude != 1) {
                builder.append(magnitude);
            }
            builder.append(charge > 0 ? "+" : "-");
        }
        if (phase != Phase.UNKNOWN) {
            builder.append("(").append(phase.token()).append(")");
        }
        return builder.toString();
    }

    public Species withCoefficient(int newCoefficient) {
        return new Species(newCoefficient, originalInput, formulaAst, composition, charge, phase);
    }
}
