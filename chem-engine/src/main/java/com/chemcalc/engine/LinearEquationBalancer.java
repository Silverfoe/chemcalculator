package com.chemcalc.engine;

import com.chemcalc.domain.BalanceConservationCheck;
import com.chemcalc.domain.BalanceConservedQuantity;
import com.chemcalc.domain.BalanceMatrixSnapshot;
import com.chemcalc.domain.BalanceNormalizationTrace;
import com.chemcalc.domain.BalanceRowOperation;
import com.chemcalc.domain.BalanceSolutionTrace;
import com.chemcalc.domain.BalanceTrace;
import com.chemcalc.domain.BalancedReaction;
import com.chemcalc.domain.ExplanationStep;
import com.chemcalc.domain.Reaction;
import com.chemcalc.domain.Species;
import com.chemcalc.domain.ValidationMessage;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public final class LinearEquationBalancer implements BalancerService {
    @Override
    public BalancedReaction balance(Reaction reaction) {
        List<Species> species = new ArrayList<>();
        species.addAll(reaction.reactants());
        species.addAll(reaction.products());

        Set<String> elementSet = new LinkedHashSet<>();
        for (Species item : species) {
            elementSet.addAll(item.composition().keySet());
        }
        List<String> elements = List.copyOf(elementSet);
        List<BalanceConservedQuantity> quantities = buildQuantities(elements);
        List<String> rowLabels = quantities.stream().map(BalanceConservedQuantity::label).toList();
        List<String> columnLabels = buildVariableLabels(reaction);

        Rational[][] matrix = buildConservationMatrix(reaction, species, elements);
        BalanceMatrixSnapshot initialMatrix = snapshot(rowLabels, columnLabels, matrix);

        SolveResult solveResult = solveNullspace(matrix, rowLabels, columnLabels);
        NormalizationResult normalizationResult = normalizeToSmallestIntegers(solveResult.solution());
        int[] coefficients = normalizationResult.coefficients();

        List<Integer> reactantCoefficients = new ArrayList<>();
        List<Integer> productCoefficients = new ArrayList<>();
        for (int index = 0; index < reaction.reactants().size(); index++) {
            reactantCoefficients.add(coefficients[index]);
        }
        for (int index = reaction.reactants().size(); index < coefficients.length; index++) {
            productCoefficients.add(coefficients[index]);
        }

        List<BalanceConservationCheck> conservationChecks = buildConservationChecks(
                elements, reaction, reactantCoefficients, productCoefficients);

        BalanceTrace trace = new BalanceTrace(
                quantities,
                initialMatrix,
                solveResult.rowOperations(),
                buildSolutionTrace(columnLabels, solveResult),
                buildNormalizationTrace(normalizationResult),
                conservationChecks);

        List<ExplanationStep> explanation = buildExplanation(trace, reaction, coefficients);
        List<ValidationMessage> messages = new ArrayList<>();
        if (!conserves(conservationChecks)) {
            messages.add(ValidationMessage.error("BALANCE_FAILED", "Balanced coefficients do not conserve atoms and charge."));
        } else {
            messages.add(ValidationMessage.info("BALANCED", "Equation balances atoms and charge."));
        }

        return new BalancedReaction(
                reaction,
                reactantCoefficients,
                productCoefficients,
                explanation,
                messages,
                trace);
    }

    private Rational[][] buildConservationMatrix(Reaction reaction, List<Species> species, List<String> elements) {
        int rows = elements.size() + 1;
        int columns = species.size();
        Rational[][] matrix = new Rational[rows][columns];
        for (int row = 0; row < rows; row++) {
            Arrays.fill(matrix[row], Rational.ZERO);
        }

        for (int column = 0; column < columns; column++) {
            Species current = species.get(column);
            int sign = column < reaction.reactants().size() ? 1 : -1;
            for (int row = 0; row < elements.size(); row++) {
                int count = current.composition().getOrDefault(elements.get(row), 0);
                matrix[row][column] = Rational.of((long) sign * count);
            }
            matrix[rows - 1][column] = Rational.of((long) sign * current.charge());
        }
        return matrix;
    }

    private List<BalanceConservedQuantity> buildQuantities(List<String> elements) {
        List<BalanceConservedQuantity> quantities = new ArrayList<>();
        for (String element : elements) {
            quantities.add(new BalanceConservedQuantity(element, "element"));
        }
        quantities.add(new BalanceConservedQuantity("Charge", "charge"));
        return List.copyOf(quantities);
    }

    private List<String> buildVariableLabels(Reaction reaction) {
        List<String> labels = new ArrayList<>();
        int variableIndex = 1;
        for (Species species : reaction.reactants()) {
            labels.add("x" + variableIndex++ + " " + species.displayFormula());
        }
        for (Species species : reaction.products()) {
            labels.add("x" + variableIndex++ + " " + species.displayFormula());
        }
        return List.copyOf(labels);
    }

    private SolveResult solveNullspace(Rational[][] inputMatrix, List<String> rowLabels, List<String> columnLabels) {
        int rows = inputMatrix.length;
        int columns = inputMatrix[0].length;
        Rational[][] matrix = copyMatrix(inputMatrix);
        List<BalanceRowOperation> operations = new ArrayList<>();
        List<Integer> pivotColumns = new ArrayList<>();
        int pivotRow = 0;

        for (int column = 0; column < columns && pivotRow < rows; column++) {
            int selectedRow = -1;
            for (int row = pivotRow; row < rows; row++) {
                if (!matrix[row][column].isZero()) {
                    selectedRow = row;
                    break;
                }
            }
            if (selectedRow == -1) {
                continue;
            }

            if (selectedRow != pivotRow) {
                swapRows(matrix, pivotRow, selectedRow);
                operations.add(new BalanceRowOperation(
                        "Swap R" + (pivotRow + 1) + " and R" + (selectedRow + 1),
                        "Move a non-zero pivot into the " + rowLabels.get(pivotRow) + " row.",
                        snapshot(rowLabels, columnLabels, matrix)));
            }

            Rational pivot = matrix[pivotRow][column];
            if (!pivot.equals(Rational.ONE)) {
                for (int currentColumn = column; currentColumn < columns; currentColumn++) {
                    matrix[pivotRow][currentColumn] = matrix[pivotRow][currentColumn].divide(pivot);
                }
                operations.add(new BalanceRowOperation(
                        "Normalize R" + (pivotRow + 1),
                        "Divide row " + (pivotRow + 1) + " by the pivot " + pivot + ".",
                        snapshot(rowLabels, columnLabels, matrix)));
            }

            for (int row = 0; row < rows; row++) {
                if (row == pivotRow || matrix[row][column].isZero()) {
                    continue;
                }
                Rational factor = matrix[row][column];
                for (int currentColumn = column; currentColumn < columns; currentColumn++) {
                    matrix[row][currentColumn] = matrix[row][currentColumn]
                            .subtract(factor.multiply(matrix[pivotRow][currentColumn]));
                }
                operations.add(new BalanceRowOperation(
                        "Eliminate R" + (row + 1),
                        "R" + (row + 1) + " = R" + (row + 1) + " - (" + factor + " * R" + (pivotRow + 1) + ").",
                        snapshot(rowLabels, columnLabels, matrix)));
            }

            pivotColumns.add(column);
            pivotRow++;
        }

        List<Integer> freeColumns = new ArrayList<>();
        for (int column = 0; column < columns; column++) {
            if (!pivotColumns.contains(column)) {
                freeColumns.add(column);
            }
        }
        if (freeColumns.isEmpty()) {
            throw new IllegalStateException("No free variable found while balancing reaction");
        }

        Rational[] solution = new Rational[columns];
        Arrays.fill(solution, Rational.ZERO);
        for (Integer freeColumn : freeColumns) {
            solution[freeColumn] = Rational.ONE;
        }
        for (int row = pivotColumns.size() - 1; row >= 0; row--) {
            int column = pivotColumns.get(row);
            Rational value = Rational.ZERO;
            for (int currentColumn = column + 1; currentColumn < columns; currentColumn++) {
                if (!matrix[row][currentColumn].isZero()) {
                    value = value.add(matrix[row][currentColumn].multiply(solution[currentColumn]));
                }
            }
            solution[column] = value.negate();
        }

        return new SolveResult(solution, List.copyOf(pivotColumns), List.copyOf(freeColumns), List.copyOf(operations));
    }

    private NormalizationResult normalizeToSmallestIntegers(Rational[] vector) {
        BigInteger leastCommonMultiple = BigInteger.ONE;
        for (Rational value : vector) {
            leastCommonMultiple = leastCommonMultiple.multiply(value.denominator())
                    .divide(leastCommonMultiple.gcd(value.denominator()));
        }

        BigInteger[] scaled = new BigInteger[vector.length];
        for (int index = 0; index < vector.length; index++) {
            scaled[index] = vector[index].numerator().multiply(leastCommonMultiple).divide(vector[index].denominator());
        }

        BigInteger[] signNormalized = Arrays.copyOf(scaled, scaled.length);
        boolean signFlipped = Arrays.stream(signNormalized).allMatch(value -> value.signum() <= 0);
        if (signFlipped) {
            for (int index = 0; index < signNormalized.length; index++) {
                signNormalized[index] = signNormalized[index].negate();
            }
        }

        for (BigInteger value : signNormalized) {
            if (value.signum() <= 0) {
                throw new IllegalStateException("Unable to derive a positive coefficient vector");
            }
        }

        BigInteger gcd = Arrays.stream(signNormalized)
                .reduce(BigInteger.ZERO, (left, right) -> left.equals(BigInteger.ZERO) ? right.abs() : left.gcd(right.abs()));
        int[] coefficients = new int[signNormalized.length];
        for (int index = 0; index < signNormalized.length; index++) {
            coefficients[index] = signNormalized[index].divide(gcd).intValueExact();
        }

        return new NormalizationResult(
                leastCommonMultiple,
                Arrays.copyOf(scaled, scaled.length),
                signFlipped,
                signNormalized,
                gcd,
                coefficients);
    }

    private List<BalanceConservationCheck> buildConservationChecks(
            List<String> elements,
            Reaction reaction,
            List<Integer> reactantCoefficients,
            List<Integer> productCoefficients) {
        Map<String, Integer> left = accumulate(reaction.reactants(), reactantCoefficients);
        Map<String, Integer> right = accumulate(reaction.products(), productCoefficients);
        List<BalanceConservationCheck> checks = new ArrayList<>();
        for (String element : elements) {
            checks.add(new BalanceConservationCheck(
                    element,
                    left.getOrDefault(element, 0),
                    right.getOrDefault(element, 0)));
        }
        checks.add(new BalanceConservationCheck(
                "Charge",
                charge(reaction.reactants(), reactantCoefficients),
                charge(reaction.products(), productCoefficients)));
        return List.copyOf(checks);
    }

    private List<ExplanationStep> buildExplanation(BalanceTrace trace, Reaction reaction, int[] coefficients) {
        List<ExplanationStep> steps = new ArrayList<>();
        steps.add(new ExplanationStep(
                "Build conservation matrix",
                "Tracked " + trace.quantities().stream().map(BalanceConservedQuantity::label).collect(Collectors.joining(", "))
                        + " across " + (reaction.reactants().size() + reaction.products().size()) + " species."));
        steps.add(new ExplanationStep(
                "Row-reduce the system",
                "Performed " + trace.rowOperations().size() + " exact row operations to isolate "
                        + trace.solution().pivotVariables().size() + " pivot variables and "
                        + trace.solution().freeVariables().size() + " free variables."));
        steps.add(new ExplanationStep(
                "Solve coefficient ratios",
                trace.solution().detail() + " Nullspace vector: " + formatVector(trace.solution().nullspaceVector()) + "."));
        steps.add(new ExplanationStep(
                "Normalize integers",
                "LCM " + trace.normalization().leastCommonMultiple()
                        + ", GCD " + trace.normalization().greatestCommonDivisor()
                        + ", final coefficients " + formatIntVector(coefficients) + "."));
        steps.add(new ExplanationStep(
                "Verify conservation",
                trace.conservationChecks().stream()
                        .map(check -> check.label() + " " + check.reactantTotal() + "=" + check.productTotal())
                        .collect(Collectors.joining("; "))));
        return List.copyOf(steps);
    }

    private BalanceSolutionTrace buildSolutionTrace(List<String> columnLabels, SolveResult solveResult) {
        List<String> pivotVariables = solveResult.pivotColumns().stream().map(columnLabels::get).toList();
        List<String> freeVariables = solveResult.freeColumns().stream().map(columnLabels::get).toList();
        List<String> nullspaceVector = Arrays.stream(solveResult.solution()).map(Rational::toString).toList();
        String detail = "Set " + String.join(", ", freeVariables) + " to 1 and back-substituted the pivot variables.";
        return new BalanceSolutionTrace(pivotVariables, freeVariables, nullspaceVector, detail);
    }

    private BalanceNormalizationTrace buildNormalizationTrace(NormalizationResult normalizationResult) {
        return new BalanceNormalizationTrace(
                normalizationResult.leastCommonMultiple().toString(),
                Arrays.stream(normalizationResult.scaledVector()).map(BigInteger::toString).toList(),
                normalizationResult.signFlipped(),
                Arrays.stream(normalizationResult.signNormalizedVector()).map(BigInteger::toString).toList(),
                normalizationResult.greatestCommonDivisor().toString(),
                Arrays.stream(normalizationResult.coefficients()).boxed().toList());
    }

    private BalanceMatrixSnapshot snapshot(List<String> rowLabels, List<String> columnLabels, Rational[][] matrix) {
        List<List<String>> values = new ArrayList<>();
        for (Rational[] row : matrix) {
            values.add(Arrays.stream(row).map(Rational::toString).toList());
        }
        return new BalanceMatrixSnapshot(rowLabels, columnLabels, values);
    }

    private Rational[][] copyMatrix(Rational[][] matrix) {
        Rational[][] copy = new Rational[matrix.length][matrix[0].length];
        for (int row = 0; row < matrix.length; row++) {
            System.arraycopy(matrix[row], 0, copy[row], 0, matrix[row].length);
        }
        return copy;
    }

    private void swapRows(Rational[][] matrix, int first, int second) {
        if (first == second) {
            return;
        }
        Rational[] temporary = matrix[first];
        matrix[first] = matrix[second];
        matrix[second] = temporary;
    }

    private boolean conserves(List<BalanceConservationCheck> checks) {
        return checks.stream().allMatch(check -> check.reactantTotal() == check.productTotal());
    }

    private Map<String, Integer> accumulate(List<Species> side, List<Integer> coefficients) {
        Map<String, Integer> total = new LinkedHashMap<>();
        for (int index = 0; index < side.size(); index++) {
            Species species = side.get(index);
            int coefficient = coefficients.get(index);
            for (Map.Entry<String, Integer> entry : species.composition().entrySet()) {
                total.merge(entry.getKey(), entry.getValue() * coefficient, Integer::sum);
            }
        }
        return total;
    }

    private int charge(List<Species> side, List<Integer> coefficients) {
        int total = 0;
        for (int index = 0; index < side.size(); index++) {
            total += side.get(index).charge() * coefficients.get(index);
        }
        return total;
    }

    private String formatVector(List<String> values) {
        return values.isEmpty() ? "[]" : "[" + String.join(", ", values) + "]";
    }

    private String formatIntVector(int[] coefficients) {
        return Arrays.stream(coefficients).mapToObj(Integer::toString).collect(Collectors.joining(", ", "[", "]"));
    }

    private record SolveResult(
            Rational[] solution,
            List<Integer> pivotColumns,
            List<Integer> freeColumns,
            List<BalanceRowOperation> rowOperations) {
    }

    private record NormalizationResult(
            BigInteger leastCommonMultiple,
            BigInteger[] scaledVector,
            boolean signFlipped,
            BigInteger[] signNormalizedVector,
            BigInteger greatestCommonDivisor,
            int[] coefficients) {
    }
}
