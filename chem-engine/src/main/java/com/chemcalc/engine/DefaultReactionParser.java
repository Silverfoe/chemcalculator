package com.chemcalc.engine;

import com.chemcalc.domain.Reaction;
import com.chemcalc.domain.ReactionArrow;
import com.chemcalc.domain.Species;
import com.chemcalc.domain.ValidationMessage;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class DefaultReactionParser implements ReactionParser {
    private static final Pattern ARROW_PATTERN = Pattern.compile("(<->|⇌|->|→|=)");
    private final FormulaParser formulaParser;

    public DefaultReactionParser(FormulaParser formulaParser) {
        this.formulaParser = formulaParser;
    }

    @Override
    public Reaction parse(String input) {
        if (input == null || input.isBlank()) {
            throw new IllegalArgumentException("Reaction must not be empty");
        }
        Matcher matcher = ARROW_PATTERN.matcher(input);
        if (!matcher.find()) {
            throw new IllegalArgumentException("Reaction arrow not found");
        }
        String arrowToken = matcher.group(1);
        String left = input.substring(0, matcher.start()).trim();
        String right = input.substring(matcher.end()).trim();
        if (left.isBlank() || right.isBlank()) {
            throw new IllegalArgumentException("Reaction must have reactants and products");
        }

        List<Species> reactants = splitSide(left).stream().map(formulaParser::parseSpecies).toList();
        List<Species> products = splitSide(right).stream().map(formulaParser::parseSpecies).toList();
        ReactionArrow arrow = "<->".equals(arrowToken) || "⇌".equals(arrowToken)
                ? ReactionArrow.REVERSIBLE
                : ReactionArrow.FORWARD;
        return new Reaction(reactants, products, arrow, List.of(ValidationMessage.info("PARSED", "Reaction parsed successfully")));
    }

    private List<String> splitSide(String side) {
        List<String> tokens = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        int depth = 0;
        for (int index = 0; index < side.length(); index++) {
            char token = side.charAt(index);
            if (token == '(' || token == '[' || token == '{') {
                depth++;
                current.append(token);
                continue;
            }
            if (token == ')' || token == ']' || token == '}') {
                depth--;
                current.append(token);
                continue;
            }
            if (token == '+' && depth == 0) {
                char nextSignificant = nextSignificant(side, index + 1);
                if (nextSignificant == '\0' || nextSignificant == '+' || nextSignificant == '-' || nextSignificant == '>') {
                    current.append(token);
                } else {
                    String segment = current.toString().trim();
                    if (!segment.isBlank()) {
                        tokens.add(segment);
                    }
                    current.setLength(0);
                }
                continue;
            }
            current.append(token);
        }
        String tail = current.toString().trim();
        if (!tail.isBlank()) {
            tokens.add(tail);
        }
        return tokens;
    }

    private char nextSignificant(String text, int start) {
        for (int index = start; index < text.length(); index++) {
            if (!Character.isWhitespace(text.charAt(index))) {
                return text.charAt(index);
            }
        }
        return '\0';
    }
}
