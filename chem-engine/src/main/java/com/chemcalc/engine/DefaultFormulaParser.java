package com.chemcalc.engine;

import com.chemcalc.domain.ElementNode;
import com.chemcalc.domain.FormulaAst;
import com.chemcalc.domain.FormulaPart;
import com.chemcalc.domain.FormulaSegment;
import com.chemcalc.domain.GroupNode;
import com.chemcalc.domain.Phase;
import com.chemcalc.domain.Species;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class DefaultFormulaParser implements FormulaParser {
    private static final Pattern STATE_PATTERN = Pattern.compile("^(.*)\\((aq|s|l|g)\\)\\s*$", Pattern.CASE_INSENSITIVE);

    @Override
    public FormulaAst parseFormula(String input) {
        return parseFormulaInternal(normalize(input));
    }

    @Override
    public Species parseSpecies(String input) {
        String normalized = normalize(input);
        int coefficient = 1;
        int index = 0;
        while (index < normalized.length() && Character.isDigit(normalized.charAt(index))) {
            index++;
        }
        if (index > 0 && index < normalized.length()) {
            char next = normalized.charAt(index);
            if (Character.isUpperCase(next) || next == '(' || next == '[' || next == '{') {
                coefficient = Integer.parseInt(normalized.substring(0, index));
                normalized = normalized.substring(index).trim();
            }
        }

        Phase phase = Phase.UNKNOWN;
        Matcher stateMatcher = STATE_PATTERN.matcher(normalized);
        if (stateMatcher.matches()) {
            phase = Phase.fromToken(stateMatcher.group(2));
            normalized = stateMatcher.group(1).trim();
        }

        ChargeParse chargeParse = stripCharge(normalized);
        FormulaAst formulaAst = parseFormulaInternal(chargeParse.formula());
        Map<String, Integer> composition = toComposition(formulaAst);
        return new Species(coefficient, input, formulaAst, composition, chargeParse.charge(), phase);
    }

    private FormulaAst parseFormulaInternal(String input) {
        if (input.isBlank()) {
            throw new IllegalArgumentException("Formula must not be empty");
        }
        Cursor cursor = new Cursor(input);
        List<FormulaSegment> segments = new ArrayList<>();
        while (!cursor.atEnd()) {
            int segmentMultiplier = cursor.parseSegmentMultiplier();
            List<FormulaPart> parts = parseParts(cursor, null);
            if (parts.isEmpty()) {
                throw new IllegalArgumentException("Invalid formula segment in '" + input + "'");
            }
            segments.add(new FormulaSegment(segmentMultiplier, parts));
            if (cursor.atEnd()) {
                break;
            }
            char separator = cursor.peek();
            if (separator == '.' || separator == '·') {
                cursor.advance();
                continue;
            }
            throw new IllegalArgumentException("Unexpected token '" + separator + "' in formula '" + input + "'");
        }
        return new FormulaAst(segments);
    }

    private List<FormulaPart> parseParts(Cursor cursor, Character expectedClosing) {
        List<FormulaPart> parts = new ArrayList<>();
        while (!cursor.atEnd()) {
            char current = cursor.peek();
            if (expectedClosing != null && current == expectedClosing) {
                cursor.advance();
                return parts;
            }
            if (current == '.' || current == '·') {
                if (expectedClosing != null) {
                    throw new IllegalArgumentException("Hydrate separator is not allowed inside a grouped section");
                }
                return parts;
            }
            if (isOpeningDelimiter(current)) {
                char opening = current;
                char closing = matchingDelimiter(opening);
                cursor.advance();
                List<FormulaPart> nested = parseParts(cursor, closing);
                int multiplier = cursor.parseNumberOrDefault(1);
                parts.add(new GroupNode(opening, closing, nested, multiplier));
                continue;
            }
            if (Character.isUpperCase(current)) {
                String symbol = cursor.parseElementSymbol();
                int count = cursor.parseNumberOrDefault(1);
                parts.add(new ElementNode(symbol, count));
                continue;
            }
            throw new IllegalArgumentException("Unexpected token '" + current + "' in formula");
        }
        if (expectedClosing != null) {
            throw new IllegalArgumentException("Unmatched delimiter '" + expectedClosing + "'");
        }
        return parts;
    }

    private Map<String, Integer> toComposition(FormulaAst formulaAst) {
        Map<String, Integer> composition = new LinkedHashMap<>();
        for (FormulaSegment segment : formulaAst.segments()) {
            for (FormulaPart part : segment.parts()) {
                accumulate(part, segment.multiplier(), composition);
            }
        }
        return composition;
    }

    private void accumulate(FormulaPart part, int multiplier, Map<String, Integer> composition) {
        if (part instanceof ElementNode elementNode) {
            composition.merge(elementNode.symbol(), elementNode.count() * multiplier, Integer::sum);
            return;
        }
        GroupNode groupNode = (GroupNode) part;
        int nestedMultiplier = multiplier * groupNode.multiplier();
        for (FormulaPart child : groupNode.parts()) {
            accumulate(child, nestedMultiplier, composition);
        }
    }

    private ChargeParse stripCharge(String input) {
        String normalized = input.trim();
        if (normalized.isEmpty()) {
            return new ChargeParse("", 0);
        }
        char last = normalized.charAt(normalized.length() - 1);
        if (last != '+' && last != '-') {
            return new ChargeParse(normalized, 0);
        }
        int sign = last == '+' ? 1 : -1;

        int caretIndex = normalized.lastIndexOf('^');
        if (caretIndex >= 0 && caretIndex < normalized.length() - 1) {
            String suffix = normalized.substring(caretIndex + 1);
            if (suffix.matches("\\d*[+-]")) {
                String magnitudeText = suffix.substring(0, suffix.length() - 1);
                int magnitude = magnitudeText.isBlank() ? 1 : Integer.parseInt(magnitudeText);
                return new ChargeParse(normalized.substring(0, caretIndex), sign * magnitude);
            }
        }

        int digitStart = normalized.length() - 1;
        while (digitStart > 0 && Character.isDigit(normalized.charAt(digitStart - 1))) {
            digitStart--;
        }
        String digits = normalized.substring(digitStart, normalized.length() - 1);
        String prefix = normalized.substring(0, digitStart);
        if (!digits.isBlank() && looksLikeChargeMagnitudePrefix(prefix)) {
            return new ChargeParse(prefix, sign * Integer.parseInt(digits));
        }
        return new ChargeParse(normalized.substring(0, normalized.length() - 1), sign);
    }

    private boolean looksLikeChargeMagnitudePrefix(String prefix) {
        if (prefix.isBlank()) {
            return false;
        }
        char last = prefix.charAt(prefix.length() - 1);
        if (last == ')' || last == ']' || last == '}') {
            return true;
        }
        return prefix.matches("[A-Z][a-z]?");
    }

    private boolean isOpeningDelimiter(char token) {
        return token == '(' || token == '[' || token == '{';
    }

    private char matchingDelimiter(char opening) {
        return switch (opening) {
            case '(' -> ')';
            case '[' -> ']';
            case '{' -> '}';
            default -> throw new IllegalArgumentException("Unsupported delimiter: " + opening);
        };
    }

    private String normalize(String input) {
        if (input == null) {
            throw new IllegalArgumentException("Input must not be null");
        }
        return input.trim().replace(" ", "");
    }

    private record ChargeParse(String formula, int charge) {
    }

    private static final class Cursor {
        private final String input;
        private int index;

        private Cursor(String input) {
            this.input = input;
        }

        boolean atEnd() {
            return index >= input.length();
        }

        char peek() {
            return input.charAt(index);
        }

        void advance() {
            index++;
        }

        int parseNumberOrDefault(int defaultValue) {
            int start = index;
            while (!atEnd() && Character.isDigit(peek())) {
                advance();
            }
            return start == index ? defaultValue : Integer.parseInt(input.substring(start, index));
        }

        int parseSegmentMultiplier() {
            int start = index;
            while (!atEnd() && Character.isDigit(peek())) {
                advance();
            }
            if (start == index) {
                return 1;
            }
            if (atEnd()) {
                throw new IllegalArgumentException("Dangling segment multiplier");
            }
            char next = peek();
            if (Character.isUpperCase(next) || next == '(' || next == '[' || next == '{') {
                return Integer.parseInt(input.substring(start, index));
            }
            index = start;
            return 1;
        }

        String parseElementSymbol() {
            StringBuilder symbol = new StringBuilder();
            symbol.append(peek());
            advance();
            while (!atEnd() && Character.isLowerCase(peek())) {
                symbol.append(peek());
                advance();
            }
            return symbol.toString();
        }
    }
}
