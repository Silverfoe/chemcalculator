package com.chemcalc.domain;

import java.util.List;

public record BalancedReaction(
        Reaction reaction,
        List<Integer> reactantCoefficients,
        List<Integer> productCoefficients,
        List<ExplanationStep> explanation,
        List<ValidationMessage> messages,
        BalanceTrace trace) {

    public BalancedReaction(
            Reaction reaction,
            List<Integer> reactantCoefficients,
            List<Integer> productCoefficients,
            List<ExplanationStep> explanation,
            List<ValidationMessage> messages) {
        this(reaction, reactantCoefficients, productCoefficients, explanation, messages, BalanceTrace.empty());
    }

    public BalancedReaction {
        reactantCoefficients = List.copyOf(reactantCoefficients == null ? List.of() : reactantCoefficients);
        productCoefficients = List.copyOf(productCoefficients == null ? List.of() : productCoefficients);
        explanation = List.copyOf(explanation == null ? List.of() : explanation);
        messages = List.copyOf(messages == null ? List.of() : messages);
        trace = trace == null ? BalanceTrace.empty() : trace;
    }

    public String format() {
        String left = formatSide(reaction.reactants(), reactantCoefficients);
        String right = formatSide(reaction.products(), productCoefficients);
        return left + " " + reaction.arrow().token() + " " + right;
    }

    private String formatSide(List<Species> speciesList, List<Integer> coefficients) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < speciesList.size(); i++) {
            if (i > 0) {
                builder.append(" + ");
            }
            int coefficient = coefficients.isEmpty() ? speciesList.get(i).coefficient() : coefficients.get(i);
            if (coefficient != 1) {
                builder.append(coefficient).append(" ");
            }
            builder.append(speciesList.get(i).displayFormula());
        }
        return builder.toString();
    }
}
