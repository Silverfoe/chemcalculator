import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.math.BigInteger;
import java.util.*;
import java.util.List;
import java.util.regex.*;

public class ChemicalGFMCalculatorTestBuild {
    // Atomic weights of elements (g/mol) - preserved from GFM calculator
    private static final Map<String, Double> atomicWeights = loadAtomicWeights();
    // Element names for naming (partial list)
    private static final Map<String, String> elementNames = new HashMap<>();
    // Monatomic anion names
    private static final Map<String, String> anionNames = new HashMap<>();
    // Common polyatomic ions (formula -> name and charge)
    private static class PolyIon {
        String name;
        int charge;
        PolyIon(String name, int charge) {
            this.name = name;
            this.charge = charge;
        }
    }
    private static final Map<String, PolyIon> polyIonMap = new HashMap<>();

    static {
        // Populate elementNames (common elements)
        elementNames.put("H", "Hydrogen");
        elementNames.put("He", "Helium");
        elementNames.put("Li", "Lithium");
        elementNames.put("Be", "Beryllium");
        elementNames.put("B", "Boron");
        elementNames.put("C", "Carbon");
        elementNames.put("N", "Nitrogen");
        elementNames.put("O", "Oxygen");
        elementNames.put("F", "Fluorine");
        elementNames.put("Ne", "Neon");
        elementNames.put("Na", "Sodium");
        elementNames.put("Mg", "Magnesium");
        elementNames.put("Al", "Aluminum");
        elementNames.put("Si", "Silicon");
        elementNames.put("P", "Phosphorus");
        elementNames.put("S", "Sulfur");
        elementNames.put("Cl", "Chlorine");
        elementNames.put("Ar", "Argon");
        elementNames.put("K", "Potassium");
        elementNames.put("Ca", "Calcium");
        elementNames.put("Fe", "Iron");
        elementNames.put("Cu", "Copper");
        elementNames.put("Zn", "Zinc");
        elementNames.put("Ag", "Silver");
        elementNames.put("Au", "Gold");
        elementNames.put("Hg", "Mercury");
        elementNames.put("Pb", "Lead");
        elementNames.put("I", "Iodine");
        elementNames.put("Br", "Bromine");
        // Monatomic anion names
        anionNames.put("H", "Hydride");
        anionNames.put("F", "Fluoride");
        anionNames.put("Cl", "Chloride");
        anionNames.put("Br", "Bromide");
        anionNames.put("I", "Iodide");
        anionNames.put("O", "Oxide");
        anionNames.put("S", "Sulfide");
        anionNames.put("N", "Nitride");
        anionNames.put("P", "Phosphide");
        anionNames.put("C", "Carbide");
        // Common polyatomic ions
        polyIonMap.put("NH4", new PolyIon("Ammonium", +1));
        polyIonMap.put("OH", new PolyIon("Hydroxide", -1));
        polyIonMap.put("NO3", new PolyIon("Nitrate", -1));
        polyIonMap.put("NO2", new PolyIon("Nitrite", -1));
        polyIonMap.put("SO4", new PolyIon("Sulfate", -2));
        polyIonMap.put("SO3", new PolyIon("Sulfite", -2));
        polyIonMap.put("CO3", new PolyIon("Carbonate", -2));
        polyIonMap.put("PO4", new PolyIon("Phosphate", -3));
        polyIonMap.put("PO3", new PolyIon("Phosphite", -3));
        polyIonMap.put("ClO4", new PolyIon("Perchlorate", -1));
        polyIonMap.put("ClO3", new PolyIon("Chlorate", -1));
        polyIonMap.put("ClO2", new PolyIon("Chlorite", -1));
        polyIonMap.put("ClO", new PolyIon("Hypochlorite", -1));
        polyIonMap.put("MnO4", new PolyIon("Permanganate", -1));
        polyIonMap.put("Cr2O7", new PolyIon("Dichromate", -2));
        polyIonMap.put("CrO4", new PolyIon("Chromate", -2));
        polyIonMap.put("CN", new PolyIon("Cyanide", -1));
        polyIonMap.put("O2", new PolyIon("Peroxide", -2));
        polyIonMap.put("C2H3O2", new PolyIon("Acetate", -1));
        polyIonMap.put("SCN", new PolyIon("Thiocyanate", -1));
    }

    // Utility class for Fraction arithmetic with BigInteger for exact ratios
    private static class Fraction {
        BigInteger num;
        BigInteger den;
        public Fraction(long numerator, long denominator) {
            this(BigInteger.valueOf(numerator), BigInteger.valueOf(denominator));
        }
        public Fraction(BigInteger numerator, BigInteger denominator) {
            if (denominator.signum() == 0) {
                throw new ArithmeticException("Denominator zero in fraction");
            }
            // normalize sign
            if (denominator.signum() < 0) {
                numerator = numerator.negate();
                denominator = denominator.negate();
            }
            BigInteger gcd = numerator.gcd(denominator);
            if (gcd.signum() < 0) gcd = gcd.negate();
            // reduce
            if (!gcd.equals(BigInteger.ZERO)) {
                numerator = numerator.divide(gcd);
                denominator = denominator.divide(gcd);
            }
            this.num = numerator;
            this.den = denominator;
        }
        public Fraction(int numerator) {
            this(BigInteger.valueOf(numerator), BigInteger.ONE);
        }
        public Fraction add(Fraction other) {
            BigInteger newNum = this.num.multiply(other.den).add(other.num.multiply(this.den));
            BigInteger newDen = this.den.multiply(other.den);
            return new Fraction(newNum, newDen);
        }
        public Fraction subtract(Fraction other) {
            BigInteger newNum = this.num.multiply(other.den).subtract(other.num.multiply(this.den));
            BigInteger newDen = this.den.multiply(other.den);
            return new Fraction(newNum, newDen);
        }
        public Fraction multiply(Fraction other) {
            BigInteger newNum = this.num.multiply(other.num);
            BigInteger newDen = this.den.multiply(other.den);
            return new Fraction(newNum, newDen);
        }
        public Fraction divide(Fraction other) {
            BigInteger newNum = this.num.multiply(other.den);
            BigInteger newDen = this.den.multiply(other.num);
            return new Fraction(newNum, newDen);
        }
        public BigInteger getNum() { return num; }
        public BigInteger getDen() { return den; }
        @Override
        public String toString() {
            return num + "/" + den;
        }
    }

    // Compound class to store formula info
    private static class Compound {
        String formula; // original formula string (without coefficient)
        Map<String, Integer> composition;
        int charge;
        Compound(String formula, Map<String,Integer> comp, int charge) {
            this.formula = formula;
            this.composition = comp;
            this.charge = charge;
        }
    }

    public static void main(String[] args) {
        // Build GUI
        JFrame frame = new JFrame("Chemistry Tool - By: Jacob Wyrozebski");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800, 600);
        // Use tabbed pane for features
        JTabbedPane tabs = new JTabbedPane();

        // Panel for Equation Balancer
        JPanel balancePanel = new JPanel();
        balancePanel.setLayout(new BorderLayout());
        JPanel inputPanel1 = new JPanel(new FlowLayout());
        JLabel eqLabel = new JLabel("Enter unbalanced equation:");
        JTextField eqField = new JTextField(40);
        JButton balanceButton = new JButton("Balance");
        inputPanel1.add(eqLabel);
        inputPanel1.add(eqField);
        inputPanel1.add(balanceButton);
        balancePanel.add(inputPanel1, BorderLayout.NORTH);
        JTextArea balanceOutput = new JTextArea();
        balanceOutput.setEditable(false);
        balanceOutput.setLineWrap(true);
        balanceOutput.setWrapStyleWord(true);
        JScrollPane scroll1 = new JScrollPane(balanceOutput);
        balancePanel.add(scroll1, BorderLayout.CENTER);
        tabs.addTab("Equation Balancer", balancePanel);

        // Panel for GFM Calculator
        JPanel gfmPanel = new JPanel();
        gfmPanel.setLayout(new BorderLayout());
        JPanel inputPanel2 = new JPanel(new FlowLayout());
        JLabel gfmLabel = new JLabel("Enter formula(s) (comma-separated):");
        JTextField gfmField = new JTextField(30);
        JButton gfmButton = new JButton("Calculate GFM");
        inputPanel2.add(gfmLabel);
        inputPanel2.add(gfmField);
        inputPanel2.add(gfmButton);
        gfmPanel.add(inputPanel2, BorderLayout.NORTH);
        JTextArea gfmOutput = new JTextArea();
        gfmOutput.setEditable(false);
        gfmOutput.setLineWrap(true);
        gfmOutput.setWrapStyleWord(true);
        JScrollPane scroll2 = new JScrollPane(gfmOutput);
        gfmPanel.add(scroll2, BorderLayout.CENTER);
        tabs.addTab("GFM Calculator", gfmPanel);

        // Panel for Periodic Table
        JPanel tablePanel = new JPanel();
        tablePanel.setLayout(new BorderLayout());
        JTextArea tableArea = new JTextArea();
        tableArea.setEditable(false);
        // Build periodic table info
        StringBuilder tableText = new StringBuilder();
        int atomicNumber = 1;
        for (Map.Entry<String, Double> entry : atomicWeights.entrySet()) {
            String sym = entry.getKey();
            double wt = entry.getValue();
            tableText.append(atomicNumber + ". " + sym + " - " + wt + " g/mol\n");
            atomicNumber++;
        }
        tableArea.setText(tableText.toString());
        JScrollPane scroll3 = new JScrollPane(tableArea);
        tablePanel.add(scroll3, BorderLayout.CENTER);
        tabs.addTab("Periodic Table", tablePanel);

        // Panel for Compound Naming
        JPanel namePanel = new JPanel();
        namePanel.setLayout(new BorderLayout());
        JPanel inputPanel3 = new JPanel(new FlowLayout());
        JLabel nameLabel = new JLabel("Enter chemical formula:");
        JTextField nameField = new JTextField(30);
        JButton nameButton = new JButton("Name Compound");
        inputPanel3.add(nameLabel);
        inputPanel3.add(nameField);
        inputPanel3.add(nameButton);
        namePanel.add(inputPanel3, BorderLayout.NORTH);
        JTextArea nameOutput = new JTextArea();
        nameOutput.setEditable(false);
        nameOutput.setLineWrap(true);
        nameOutput.setWrapStyleWord(true);
        JScrollPane scroll4 = new JScrollPane(nameOutput);
        namePanel.add(scroll4, BorderLayout.CENTER);
        tabs.addTab("Compound Naming (Alpha)", namePanel);

        frame.add(tabs);
        frame.setVisible(true);

        // Action Listeners for buttons
        balanceButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String eq = eqField.getText().trim();
                if (eq.isEmpty()) {
                    return;
                }
                List<String> steps = balanceChemicalEquation(eq);
                // Display steps in output area
                balanceOutput.setText("");
                for (String step : steps) {
                    balanceOutput.append(step + "\n");
                }
            }
        });
        gfmButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String input = gfmField.getText().trim();
                if (input.isEmpty()) {
                    return;
                }
                StringBuilder result = new StringBuilder();
                String[] formulas = input.split("\\s*,\\s*");
                for (String formula : formulas) {
                    if (formula.isEmpty()) continue;
                    result.append("Formula: " + formula + "\n");
                    try {
                        Map<String,Integer> comp = parseFormulaComposition(formula);
                        double totalMass = 0.0;
                        for (Map.Entry<String,Integer> entry : comp.entrySet()) {
                            String element = entry.getKey();
                            int count = entry.getValue();
                            double atomicMass = atomicWeights.getOrDefault(element, 0.0);
                            double mass = atomicMass * count;
                            totalMass += mass;
                            result.append("  " + element + ": " + String.format("%.3f", mass) + " g/mol (x" + count + ")\n");
                        }
                        result.append("  Total GFM: " + String.format("%.3f", totalMass) + " g/mol\n\n");
                    } catch (Exception ex) {
                        result.append("  Error parsing formula.\n\n");
                    }
                }
                gfmOutput.setText(result.toString());
            }
        });
        nameButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String formula = nameField.getText().trim();
                if (formula.isEmpty()) {
                    return;
                }
                try {
                    Compound comp = parseCompound(formula);
                    String name = nameCompound(comp);
                    nameOutput.setText(name);
                } catch (Exception ex) {
                    nameOutput.setText("Unable to name the compound. Please check the formula.");
                }
            }
        });
    }

    // Balancing chemical equation main logic
    private static List<String> balanceChemicalEquation(String equation) {
        List<String> steps = new ArrayList<>();
        // Parse equation into reactant and product compounds
        String[] sides = equation.split("->|=");
        if (sides.length != 2) {
            steps.add("Error: Equation must have a single '->' separating reactants and products.");
            return steps;
        }
        String reactantsStr = sides[0].trim();
        String productsStr = sides[1].trim();
        if (reactantsStr.isEmpty() || productsStr.isEmpty()) {
            steps.add("Error: Reactant or product side is empty.");
            return steps;
        }
        // Split by + (assuming plus signs separate species)
        String[] reactantTokens = reactantsStr.split("\\s*\\+\\s*");
        String[] productTokens = productsStr.split("\\s*\\+\\s*");
        List<Compound> reactants = new ArrayList<>();
        List<Compound> products = new ArrayList<>();
        try {
            for (String token : reactantTokens) {
                if (token.isEmpty()) continue;
                Compound c = parseCompound(token.trim());
                reactants.add(c);
            }
            for (String token : productTokens) {
                if (token.isEmpty()) continue;
                Compound c = parseCompound(token.trim());
                products.add(c);
            }
        } catch (Exception ex) {
            steps.add("Error: Failed to parse equation components.");
            return steps;
        }
        // Identify if redox
        boolean isRedox = false;
        // Map each element to oxidation numbers on each side
        Map<String, Integer> oxLeft = new HashMap<>();
        Map<String, Integer> oxRight = new HashMap<>();
        // Compute oxidation states for each compound and record any changes
        for (Compound comp : reactants) {
            Map<String, Integer> oxMap = assignOxidationNumbers(comp);
            for (Map.Entry<String,Integer> entry : oxMap.entrySet()) {
                String elem = entry.getKey();
                int ox = entry.getValue();
                if (!oxLeft.containsKey(elem)) {
                    oxLeft.put(elem, ox);
                } else {
                    // If element appears in multiple reactants with different states, mark redox
                    if (oxLeft.get(elem) != ox) {
                        isRedox = true;
                    }
                }
            }
        }
        for (Compound comp : products) {
            Map<String, Integer> oxMap = assignOxidationNumbers(comp);
            for (Map.Entry<String,Integer> entry : oxMap.entrySet()) {
                String elem = entry.getKey();
                int ox = entry.getValue();
                if (!oxRight.containsKey(elem)) {
                    oxRight.put(elem, ox);
                } else {
                    if (oxRight.get(elem) != ox) {
                        isRedox = true;
                    }
                }
            }
        }
        // Compare oxidation states between left and right
        for (String elem : oxLeft.keySet()) {
            if (!oxRight.containsKey(elem)) continue;
            if (oxLeft.get(elem).intValue() != oxRight.get(elem).intValue()) {
                isRedox = true;
            }
        }
        if (isRedox) {
            steps.add("Redox reaction detected (oxidation states change). Using half-reaction method:");
            steps.addAll(balanceRedoxReaction(reactants, products));
        } else {
            steps.add("No change in oxidation numbers detected (non-redox reaction). Using algebraic balancing:");
            steps.addAll(balanceNonRedoxReaction(reactants, products));
        }
        return steps;
    }

    // Balance redox reactions via half-reaction method
    private static List<String> balanceRedoxReaction(List<Compound> reactants, List<Compound> products) {
        List<String> steps = new ArrayList<>();
        // Identify elements that change oxidation number
        Map<String,Integer> startOx = new HashMap<>();
        Map<String,Integer> endOx = new HashMap<>();
        for (Compound comp : reactants) {
            Map<String,Integer> ox = assignOxidationNumbers(comp);
            startOx.putAll(ox);
        }
        for (Compound comp : products) {
            Map<String,Integer> ox = assignOxidationNumbers(comp);
            endOx.putAll(ox);
        }
        List<String> oxidizedElements = new ArrayList<>();
        List<String> reducedElements = new ArrayList<>();
        for (String elem : startOx.keySet()) {
            if (!endOx.containsKey(elem)) continue;
            int oxStart = startOx.get(elem);
            int oxEnd = endOx.get(elem);
            if (oxEnd > oxStart) {
                oxidizedElements.add(elem);
            } else if (oxEnd < oxStart) {
                reducedElements.add(elem);
            }
        }
        // Remove duplicates
        Set<String> oxSet = new HashSet<>(oxidizedElements);
        Set<String> redSet = new HashSet<>(reducedElements);
        oxidizedElements = new ArrayList<>(oxSet);
        reducedElements = new ArrayList<>(redSet);
        // If the same element appears as both (disproportionation), ensure it's in both lists
        for (String elem : oxidizedElements) {
            if (reducedElements.contains(elem)) {
                // already in both, do nothing
            }
        }
        // Prepare half-reactions
        class HalfReaction {
            Compound reactant;
            Compound product;
            List<Compound> leftExtras = new ArrayList<>();
            List<Compound> rightExtras = new ArrayList<>();
            int electrons = 0;
            boolean electronsOnLeft = false;
        }
        List<HalfReaction> halfReactions = new ArrayList<>();
        if (oxidizedElements.size() == 1 && reducedElements.size() == 1 && oxidizedElements.get(0).equals(reducedElements.get(0))) {
            // Disproportionation: one element is both oxidized and reduced
            String elem = oxidizedElements.get(0);
            Compound reactantComp = null;
            for (Compound comp : reactants) {
                if (comp.composition.containsKey(elem)) {
                    reactantComp = comp;
                    break;
                }
            }
            List<Compound> prodComps = new ArrayList<>();
            for (Compound comp : products) {
                if (comp.composition.containsKey(elem)) {
                    prodComps.add(comp);
                }
            }
            if (reactantComp != null && prodComps.size() >= 2) {
                Compound prod1 = prodComps.get(0);
                Compound prod2 = prodComps.get(1);
                int oxReact = assignOxidationNumbers(reactantComp).get(elem);
                int oxP1 = assignOxidationNumbers(prod1).get(elem);
                int oxP2 = assignOxidationNumbers(prod2).get(elem);
                HalfReaction half1 = new HalfReaction();
                HalfReaction half2 = new HalfReaction();
                half1.reactant = reactantComp;
                half2.reactant = reactantComp;
                // Assign which product is oxidation vs reduction form
                if (Math.abs(oxP1 - oxReact) > Math.abs(oxP2 - oxReact)) {
                    // assume one with higher ox state difference is oxidation
                }
                // Determine by value:
                if (oxP1 > oxReact && oxP2 < oxReact) {
                    half1.product = prod1; // oxidation product
                    half2.product = prod2; // reduction product
                } else if (oxP2 > oxReact && oxP1 < oxReact) {
                    half1.product = prod2;
                    half2.product = prod1;
                } else {
                    // if both are either higher or lower, just assign arbitrarily
                    half1.product = prod1;
                    half2.product = prod2;
                }
                halfReactions.add(half1);
                halfReactions.add(half2);
            }
        } else {
            // Typical case: separate elements
            for (String elem : oxidizedElements) {
                Compound src = null, dst = null;
                for (Compound comp : reactants) {
                    if (comp.composition.containsKey(elem)) { src = comp; break; }
                }
                for (Compound comp : products) {
                    if (comp.composition.containsKey(elem)) { dst = comp; break; }
                }
                if (src != null && dst != null) {
                    HalfReaction half = new HalfReaction();
                    half.reactant = src;
                    half.product = dst;
                    halfReactions.add(half);
                }
            }
            for (String elem : reducedElements) {
                Compound src = null, dst = null;
                for (Compound comp : reactants) {
                    if (comp.composition.containsKey(elem)) { src = comp; break; }
                }
                for (Compound comp : products) {
                    if (comp.composition.containsKey(elem)) { dst = comp; break; }
                }
                if (src != null && dst != null) {
                    HalfReaction half = new HalfReaction();
                    half.reactant = src;
                    half.product = dst;
                    halfReactions.add(half);
                }
            }
        }
        // Balance each half-reaction
        // Determine medium (acidic or basic)
        String medium = "acidic";
        for (Compound comp : reactants) {
            if (comp.composition.containsKey("OH")) { medium = "basic"; break; }
        }
        if (medium.equals("acidic")) {
            for (Compound comp : products) {
                if (comp.composition.containsKey("OH")) { medium = "basic"; break; }
            }
        }
        for (HalfReaction half : halfReactions) {
            // Identify the key element being transferred
            String keyElem = null;
            for (String e : half.reactant.composition.keySet()) {
                if (half.product.composition.containsKey(e)) {
                    keyElem = e;
                    break;
                }
            }
            if (keyElem == null) {
                // if no common element, pick any element from reactant as key
                keyElem = half.reactant.composition.keySet().iterator().next();
            }
            steps.add("Half-reaction: " + half.reactant.formula + " -> " + half.product.formula);
            // Balance the key element
            int rCount = half.reactant.composition.getOrDefault(keyElem, 0);
            int pCount = half.product.composition.getOrDefault(keyElem, 0);
            if (rCount != pCount && rCount > 0 && pCount > 0) {
                if (rCount < pCount) {
                    // multiply reactant molecule
                    // (We conceptually adjust the coefficient, but will reflect in final combination)
                } else {
                    // multiply product molecule
                }
                // Note: We show the step but actual coefficient adjustment is applied later when combining half-reactions
                if (rCount != pCount) {
                    steps.add("Balance " + keyElem + ": " + half.reactant.formula + " -> " + half.product.formula);
                }
            }
            // Balance oxygen by adding H2O
            int oxyLeft = half.reactant.composition.getOrDefault("O", 0) + extraElementCount(half.leftExtras, "O");
            int oxyRight = half.product.composition.getOrDefault("O", 0) + extraElementCount(half.rightExtras, "O");
            if (oxyLeft != oxyRight) {
                if (oxyLeft < oxyRight) {
                    int diff = oxyRight - oxyLeft;
                    Compound water = parseCompound(diff + "H2O");
                    half.leftExtras.add(water);
                } else {
                    int diff = oxyLeft - oxyRight;
                    Compound water = parseCompound(diff + "H2O");
                    half.rightExtras.add(water);
                }
                steps.add("Balance O with H2O: " + formatHalfReaction(half));
            }
            // Balance hydrogen by adding H+ (acidic) or H2O/OH- (basic)
            int hydLeft = half.reactant.composition.getOrDefault("H", 0) + extraElementCount(half.leftExtras, "H");
            int hydRight = half.product.composition.getOrDefault("H", 0) + extraElementCount(half.rightExtras, "H");
            if (hydLeft != hydRight) {
                if (medium.equals("acidic")) {
                    if (hydLeft < hydRight) {
                        int diff = hydRight - hydLeft;
                        Compound proton = parseCompound(diff + "H+");
                        half.leftExtras.add(proton);
                    } else {
                        int diff = hydLeft - hydRight;
                        Compound proton = parseCompound(diff + "H+");
                        half.rightExtras.add(proton);
                    }
                    steps.add("Balance H with H+: " + formatHalfReaction(half));
                } else { // basic
                    if (hydLeft < hydRight) {
                        int diff = hydRight - hydLeft;
                        // In basic solution, add equal OH- to both sides for each H+ that would be needed
                        // Here, we effectively add H2O to left (to supply H) and OH- to right
                        Compound water = parseCompound(diff + "H2O");
                        Compound hydroxide = parseCompound(diff + "OH-");
                        half.leftExtras.add(water);
                        half.rightExtras.add(hydroxide);
                    } else {
                        int diff = hydLeft - hydRight;
                        Compound water = parseCompound(diff + "H2O");
                        Compound hydroxide = parseCompound(diff + "OH-");
                        half.rightExtras.add(water);
                        half.leftExtras.add(hydroxide);
                    }
                    steps.add("Balance H in basic solution (H2O/OH-): " + formatHalfReaction(half));
                }
            }
            // Balance charge by adding electrons
            int chargeLeft = totalCharge(half.reactant) + extraCharge(half.leftExtras);
            int chargeRight = totalCharge(half.product) + extraCharge(half.rightExtras);
            if (chargeLeft != chargeRight) {
                if (chargeLeft > chargeRight) {
                    half.electrons = chargeLeft - chargeRight;
                    half.electronsOnLeft = true;
                } else {
                    half.electrons = chargeRight - chargeLeft;
                    half.electronsOnLeft = false;
                }
                steps.add("Balance charge with e-: " + formatHalfReaction(half));
            }
        }
        // Now equalize electrons between the two half-reactions
        if (halfReactions.size() == 2) {
            HalfReaction h1 = halfReactions.get(0);
            HalfReaction h2 = halfReactions.get(1);
            int e1 = h1.electrons;
            int e2 = h2.electrons;
            if (e1 != e2 && e1 != 0 && e2 != 0) {
                int lcm = lcm(e1, e2);
                int factor1 = lcm / e1;
                int factor2 = lcm / e2;
                steps.add("Multiply half-reactions to equalize electrons: oxidation x" + factor1 + ", reduction x" + factor2);
            }
            // Combine half-reactions
            Map<String, Integer> leftMap = new LinkedHashMap<>();
            Map<String, Integer> rightMap = new LinkedHashMap<>();
            for (HalfReaction half : halfReactions) {
                addSpecies(leftMap, half.reactant.formula, 1);
                for (Compound extra : half.leftExtras) {
                    addSpecies(leftMap, extra.formula, 1);
                }
                addSpecies(rightMap, half.product.formula, 1);
                for (Compound extra : half.rightExtras) {
                    addSpecies(rightMap, extra.formula, 1);
                }
                if (half.electrons > 0) {
                    String eSymbol = "e-";
                    if (half.electronsOnLeft) {
                        addSpecies(leftMap, eSymbol, half.electrons);
                    } else {
                        addSpecies(rightMap, eSymbol, half.electrons);
                    }
                }
            }
            // Cancel electrons and other common species
            cancelSpecies(leftMap, rightMap, "e-");
            cancelSpecies(leftMap, rightMap, "H2O");
            cancelSpecies(leftMap, rightMap, "H+");
            cancelSpecies(leftMap, rightMap, "OH-");
            // Format final balanced equation
            String finalEq = formatEquation(leftMap, rightMap);
            steps.add("Balanced Equation: " + finalEq);
        } else {
            // Fallback to algebraic if half-reaction method did not produce exactly two halves
            steps.addAll(balanceNonRedoxReaction(reactants, products));
        }
        return steps;
    }

    // Balance non-redox reaction using linear algebra (matrix method)
    private static List<String> balanceNonRedoxReaction(List<Compound> reactants, List<Compound> products) {
        List<String> steps = new ArrayList<>();
        // Combine all compounds and list all elements
        List<Compound> allCompounds = new ArrayList<>();
        allCompounds.addAll(reactants);
        allCompounds.addAll(products);
        int n = allCompounds.size();
        int reactantCount = reactants.size();
        Set<String> elementSet = new HashSet<>();
        for (Compound comp : allCompounds) {
            elementSet.addAll(comp.composition.keySet());
        }
        List<String> elements = new ArrayList<>(elementSet);
        int m = elements.size();
        // Matrix (m x (n-1)) for unknown coefficients (excluding first compound)
        Fraction[][] A = new Fraction[m][n-1];
        Fraction[] B = new Fraction[m];
        // We fix the first compound’s coefficient = 1 and move its contribution to constants
        for (int i = 0; i < m; i++) {
            String elem = elements.get(i);
            // Compute total from first compound
            int contribFirst = 0;
            if (reactants.contains(allCompounds.get(0))) {
                contribFirst = allCompounds.get(0).composition.getOrDefault(elem, 0);
            } else {
                contribFirst = - allCompounds.get(0).composition.getOrDefault(elem, 0);
            }
            B[i] = new Fraction(-contribFirst, 1);
            // Fill coefficients for j = 1..n-1 compounds
            for (int j = 1; j < n; j++) {
                Compound comp = allCompounds.get(j);
                int count = comp.composition.getOrDefault(elem, 0);
                if (reactants.contains(comp)) {
                    A[i][j-1] = new Fraction(count, 1);
                } else {
                    A[i][j-1] = new Fraction(-count, 1);
                }
            }
        }
        // Gaussian elimination to solve A * x = B
        int unknowns = n - 1;
        int row = 0, col = 0;
        int[] pivotCol = new int[m];
        Arrays.fill(pivotCol, -1);
        while (row < m && col < unknowns) {
            int pivot = row;
            while (pivot < m && A[pivot][col].num.equals(BigInteger.ZERO)) {
                pivot++;
            }
            if (pivot == m) {
                col++;
                continue;
            }
            if (pivot != row) {
                // swap rows
                Fraction[] tempRow = A[row];
                A[row] = A[pivot];
                A[pivot] = tempRow;
                Fraction tempB = B[row];
                B[row] = B[pivot];
                B[pivot] = tempB;
            }
            pivotCol[row] = col;
            // Normalize pivot to 1
            Fraction pivotVal = A[row][col];
            for (int j = col; j < unknowns; j++) {
                A[row][j] = A[row][j].divide(pivotVal);
            }
            B[row] = B[row].divide(pivotVal);
            // Eliminate other rows
            for (int i = 0; i < m; i++) {
                if (i != row && !A[i][col].num.equals(BigInteger.ZERO)) {
                    Fraction factor = A[i][col];
                    for (int j = col; j < unknowns; j++) {
                        A[i][j] = A[i][j].subtract(factor.multiply(A[row][j]));
                    }
                    B[i] = B[i].subtract(factor.multiply(B[row]));
                }
            }
            row++;
            col++;
        }
        // Back-substitution
        Fraction[] solution = new Fraction[unknowns];
        Arrays.fill(solution, new Fraction(0));
        for (int i = m - 1; i >= 0; i--) {
            if (pivotCol[i] != -1) {
                int pc = pivotCol[i];
                Fraction sum = new Fraction(0);
                for (int j = pc + 1; j < unknowns; j++) {
                    sum = sum.add(A[i][j].multiply(solution[j]));
                }
                solution[pc] = B[i].subtract(sum);
            }
        }
        // Include first coefficient as 1 (in fraction form)
        Fraction firstCoef = new Fraction(1);
        // Scale all coefficients by LCM of denominators
        BigInteger lcmDen = BigInteger.ONE;
        for (Fraction frac : solution) {
            lcmDen = lcmDen.multiply(frac.den).divide(lcmDen.gcd(frac.den));
        }
        BigInteger scale = lcmDen;
        BigInteger[] coeffs = new BigInteger[n];
        coeffs[0] = scale.multiply(firstCoef.num).divide(firstCoef.den);
        for (int j = 1; j < n; j++) {
            coeffs[j] = scale.multiply(solution[j-1].num).divide(solution[j-1].den);
        }
        // Simplify by GCD of all coefficients
        BigInteger gcdAll = coeffs[0].abs();
        for (int j = 1; j < n; j++) {
            gcdAll = gcdAll.gcd(coeffs[j].abs());
        }
        if (!gcdAll.equals(BigInteger.ZERO)) {
            for (int j = 0; j < n; j++) {
                coeffs[j] = coeffs[j].divide(gcdAll);
            }
        }
        // Construct balanced equation string
        StringBuilder eqStr = new StringBuilder();
        for (int j = 0; j < reactantCount; j++) {
            Compound comp = reactants.get(j);
            BigInteger coef = coeffs[j];
            if (j > 0) eqStr.append(" + ");
            if (!coef.equals(BigInteger.ONE)) eqStr.append(coef + " ");
            eqStr.append(comp.formula);
        }
        eqStr.append(" -> ");
        for (int k = 0; k < products.size(); k++) {
            Compound comp = products.get(k);
            BigInteger coef = coeffs[reactantCount + k];
            if (k > 0) eqStr.append(" + ");
            if (!coef.equals(BigInteger.ONE)) eqStr.append(coef + " ");
            eqStr.append(comp.formula);
        }
        steps.add("Balanced Equation: " + eqStr.toString());
        return steps;
    }

    // Parse a formula string into a Compound (with composition and charge)
    private static Compound parseCompound(String formulaStr) {
        String formula = formulaStr.trim();
        // Remove leading coefficient if present (e.g., "2 H2O" -> "H2O")
        if (!formula.isEmpty() && Character.isDigit(formula.charAt(0))) {
            int idx = 0;
            while (idx < formula.length() && Character.isDigit(formula.charAt(idx))) idx++;
            if (idx < formula.length() && formula.charAt(idx) == ' ') idx++;
            formula = formula.substring(idx);
        }
        // Extract charge at end if present
        int charge = 0;
        Pattern chargePattern = Pattern.compile("(.*?)(?:\\^(\\d+)?([+-])|\\((\\d+)([+-])\\)|([+-]))$");
        Matcher matcher = chargePattern.matcher(formula);
        String coreFormula = formula;
        if (matcher.matches()) {
            coreFormula = matcher.group(1);
            String number = null;
            String sign = null;
            if (matcher.group(3) != null) { // case with caret
                number = matcher.group(2);
                sign = matcher.group(3);
            } else if (matcher.group(5) != null) { // case with parentheses
                number = matcher.group(4);
                sign = matcher.group(5);
            } else if (matcher.group(6) != null) { // single + or -
                number = "";
                sign = matcher.group(6);
            }
            int magnitude = (number == null || number.isEmpty()) ? 1 : Integer.parseInt(number);
            if (sign != null) {
                charge = sign.equals("+") ? magnitude : -magnitude;
            }
        }
        Map<String,Integer> composition = parseFormulaComposition(coreFormula);
        return new Compound(coreFormula, composition, charge);
    }

    // Parse formula (with parentheses) into element composition
    private static Map<String,Integer> parseFormulaComposition(String formula) {
        Map<String,Integer> comp = new HashMap<>();
        parseFormulaRecursive(formula, 1, comp);
        return comp;
    }
    private static void parseFormulaRecursive(String formula, int multiplier, Map<String,Integer> comp) {
        int i = 0;
        while (i < formula.length()) {
            char ch = formula.charAt(i);
            if (ch == '(') {
                // find matching parenthesis
                int depth = 1;
                int j = i + 1;
                while (j < formula.length() && depth > 0) {
                    if (formula.charAt(j) == '(') depth++;
                    if (formula.charAt(j) == ')') depth--;
                    j++;
                }
                if (depth != 0) throw new IllegalArgumentException("Unmatched parentheses in formula");
                String subformula = formula.substring(i+1, j-1);
                // find if there's a numeric multiplier after ')'
                int k = j;
                StringBuilder numBuilder = new StringBuilder();
                while (k < formula.length() && Character.isDigit(formula.charAt(k))) {
                    numBuilder.append(formula.charAt(k));
                    k++;
                }
                int count = numBuilder.length() > 0 ? Integer.parseInt(numBuilder.toString()) : 1;
                parseFormulaRecursive(subformula, multiplier * count, comp);
                i = k;
            } else if (Character.isUpperCase(ch)) {
                StringBuilder elem = new StringBuilder();
                elem.append(ch);
                i++;
                while (i < formula.length() && Character.isLowerCase(formula.charAt(i))) {
                    elem.append(formula.charAt(i));
                    i++;
                }
                // parse number after element
                StringBuilder numBuilder = new StringBuilder();
                while (i < formula.length() && Character.isDigit(formula.charAt(i))) {
                    numBuilder.append(formula.charAt(i));
                    i++;
                }
                int count = numBuilder.length() > 0 ? Integer.parseInt(numBuilder.toString()) : 1;
                String element = elem.toString();
                comp.put(element, comp.getOrDefault(element, 0) + count * multiplier);
            } else {
                // Skip any unexpected characters (e.g., spaces or + which should not appear here)
                i++;
            }
        }
    }

    // Compute oxidation numbers for each element in a compound
    private static Map<String,Integer> assignOxidationNumbers(Compound comp) {
        Map<String,Integer> oxMap = new HashMap<>();
        Map<String,Integer> compMap = comp.composition;
        int netCharge = comp.charge;
        if (compMap.size() == 0) return oxMap;
        if (compMap.size() == 1) {
            // Single-element species
            String elem = compMap.keySet().iterator().next();
            if (netCharge == 0) {
                // Elemental form
                oxMap.put(elem, 0);
            } else {
                int count = compMap.get(elem);
                oxMap.put(elem, netCharge / count);
            }
            return oxMap;
        }
        int sumKnown = 0;
        List<String> unknowns = new ArrayList<>();
        for (String elem : compMap.keySet()) {
            int count = compMap.get(elem);
            int ox = 0;
            if (elem.equals("H")) {
                // Hydrogen: +1 (with nonmetals) or -1 (with metals)
                boolean withMetal = false;
                for (String e2 : compMap.keySet()) {
                    if (!e2.equals("H") && isMetal(e2)) { withMetal = true; break; }
                }
                ox = withMetal ? -1 : +1;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("O")) {
                // Oxygen: assume -2 initially
                ox = -2;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("F")) {
                ox = -1;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (isGroup1(elem)) {
                ox = +1;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (isGroup2(elem)) {
                ox = +2;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("Al")) {
                ox = +3;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("Zn")) {
                ox = +2;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("Ag")) {
                ox = +1;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (anionNames.containsKey(elem)) {
                // Halogens or other nonmetals
                if ((elem.equals("Cl")||elem.equals("Br")||elem.equals("I")) && compMap.containsKey("O")) {
                    // Halogen in presence of O (likely positive oxidation state, leave unknown)
                    unknowns.add(elem);
                } else if ((elem.equals("S")||elem.equals("Se")||elem.equals("Te")) && compMap.containsKey("O")) {
                    // Chalcogen in oxyanion
                    unknowns.add(elem);
                } else if ((elem.equals("N")||elem.equals("P")) && compMap.containsKey("O")) {
                    unknowns.add(elem);
                } else {
                    // Otherwise assign typical negative charge
                    if (elem.equals("Cl")||elem.equals("Br")||elem.equals("I")) ox = -1;
                    else if (elem.equals("S")||elem.equals("Se")||elem.equals("Te")) ox = -2;
                    else if (elem.equals("N")||elem.equals("P")) ox = -3;
                    else ox = -1;
                    oxMap.put(elem, ox);
                    sumKnown += ox * count;
                }
            } else {
                // Element with variable oxidation (transition metal or others)
                unknowns.add(elem);
            }
        }
        if (!unknowns.isEmpty()) {
            if (unknowns.size() == 1) {
                String elem = unknowns.get(0);
                int count = compMap.get(elem);
                int ox = (netCharge - sumKnown) / count;
                oxMap.put(elem, ox);
            } else {
                // If multiple unknown oxidation states (e.g., organic molecules), we may not determine all.
                for (String elem : unknowns) {
                    // As a fallback, assign 0 for unknown (they might change in redox)
                    oxMap.put(elem, 0);
                }
            }
        }
        // Adjust for special cases like peroxides:
        if (compMap.containsKey("O")) {
            int Ocount = compMap.get("O");
            if (Ocount == 2 && oxMap.getOrDefault("O", -2) == -2) {
                // Check if total charge doesn't match, possibly O is -1 each
                int totalCalc = 0;
                for (String elem : compMap.keySet()) {
                    totalCalc += oxMap.getOrDefault(elem, 0) * compMap.get(elem);
                }
                if (totalCalc != netCharge) {
                    oxMap.put("O", -1);
                    // Recompute one unknown if exists
                    for (String elem : unknowns) {
                        if (!elem.equals("O")) {
                            int count = compMap.get(elem);
                            int ox = (netCharge - (totalCalc + Ocount)) / count;
                            oxMap.put(elem, ox);
                        }
                    }
                }
            }
        }
        return oxMap;
    }

    private static boolean isMetal(String element) {
        // Define metals by exclusion of known nonmetals:
        Set<String> nonmetals = new HashSet<>(Arrays.asList(
            "H","He","C","N","O","F","Ne","P","S","Cl","Ar","Se","Br","Kr","I","Xe","At","Rn","Og"
        ));
        return !nonmetals.contains(element);
    }
    private static boolean isGroup1(String element) {
        return Arrays.asList("Li","Na","K","Rb","Cs","Fr").contains(element);
    }
    private static boolean isGroup2(String element) {
        return Arrays.asList("Be","Mg","Ca","Sr","Ba","Ra").contains(element);
    }

    // Helper to compute total charge of a compound or extras list
    private static int totalCharge(Compound comp) {
        return comp.charge;
    }
    private static int extraCharge(List<Compound> extras) {
        int sum = 0;
        for (Compound c : extras) {
            sum += c.charge;
        }
        return sum;
    }
    private static int extraElementCount(List<Compound> extras, String element) {
        int sum = 0;
        for (Compound c : extras) {
            sum += c.composition.getOrDefault(element, 0);
        }
        return sum;
    }

    // Format half-reaction for output
    private static String formatHalfReaction(Object halfObj) {
        if (!(halfObj instanceof HalfReaction)) {
            return "Invalid object passed to formatHalfReaction.";
        }
        HalfReaction half = (HalfReaction) halfObj;
    
        StringBuilder left = new StringBuilder(half.reactant.formula);
        for (Compound extra : half.leftExtras) {
            left.append(" + ").append(extra.formula);
        }
    
        StringBuilder right = new StringBuilder(half.product.formula);
        for (Compound extra : half.rightExtras) {
            right.append(" + ").append(extra.formula);
        }
    
        if (half.electrons > 0) {
            if (half.electronsOnLeft) {
                left.append(" + ").append(half.electrons).append("e-");
            } else {
                right.append(" + ").append(half.electrons).append("e-");
            }
        }
    
        return left + " -> " + right;
    }

    // Add a species count to a side map
    private static void addSpecies(Map<String,Integer> map, String species, int count) {
        if (species == null || species.isEmpty()) return;
        map.put(species, map.getOrDefault(species, 0) + count);
    }
    // Cancel common species between sides
    private static void cancelSpecies(Map<String,Integer> left, Map<String,Integer> right, String species) {
        int leftCount = left.getOrDefault(species, 0);
        int rightCount = right.getOrDefault(species, 0);
        if (leftCount > 0 && rightCount > 0) {
            int cancel = Math.min(leftCount, rightCount);
            left.put(species, leftCount - cancel);
            right.put(species, rightCount - cancel);
        }
    }
    // Format full equation from left and right species maps
    private static String formatEquation(Map<String,Integer> left, Map<String,Integer> right) {
        StringBuilder sb = new StringBuilder();
        boolean firstTerm = true;
        for (Map.Entry<String,Integer> entry : left.entrySet()) {
            int count = entry.getValue();
            String species = entry.getKey();
            if (count == 0) continue;
            if (!firstTerm) sb.append(" + ");
            if (count > 1) sb.append(count + " ");
            sb.append(species);
            firstTerm = false;
        }
        sb.append(" -> ");
        firstTerm = true;
        for (Map.Entry<String,Integer> entry : right.entrySet()) {
            int count = entry.getValue();
            String species = entry.getKey();
            if (count == 0) continue;
            if (!firstTerm) sb.append(" + ");
            if (count > 1) sb.append(count + " ");
            sb.append(species);
            firstTerm = false;
        }
        return sb.toString();
    }

    // Compute least common multiple of two integers
    private static int lcm(int a, int b) {
        if (a == 0 || b == 0) return 0;
        int gcd = BigInteger.valueOf(a).gcd(BigInteger.valueOf(b)).intValue();
        return Math.abs(a / gcd * b);
    }

    // Load atomic weights into a LinkedHashMap (to preserve order by atomic number)
    private static Map<String, Double> loadAtomicWeights() {
        Map<String, Double> weights = new LinkedHashMap<>();
        weights.put("H", 1.008);
        weights.put("He", 4.0026);
        weights.put("Li", 6.94);
        weights.put("Be", 9.0122);
        weights.put("B", 10.81);
        weights.put("C", 12.011);
        weights.put("N", 14.007);
        weights.put("O", 15.999);
        weights.put("F", 18.998);
        weights.put("Ne", 20.180);
        weights.put("Na", 22.990);
        weights.put("Mg", 24.305);
        weights.put("Al", 26.982);
        weights.put("Si", 28.085);
        weights.put("P", 30.974);
        weights.put("S", 32.06);
        weights.put("Cl", 35.45);
        weights.put("K", 39.098);
        weights.put("Ar", 39.948);
        weights.put("Ca", 40.078);
        weights.put("Sc", 44.956);
        weights.put("Ti", 47.867);
        weights.put("V", 50.942);
        weights.put("Cr", 51.996);
        weights.put("Mn", 54.938);
        weights.put("Fe", 55.845);
        weights.put("Co", 58.933);
        weights.put("Ni", 58.693);
        weights.put("Cu", 63.546);
        weights.put("Zn", 65.38);
        weights.put("Ga", 69.723);
        weights.put("Ge", 72.630);
        weights.put("As", 74.922);
        weights.put("Se", 78.971);
        weights.put("Br", 79.904);
        weights.put("Kr", 83.798);
        weights.put("Rb", 85.468);
        weights.put("Sr", 87.62);
        weights.put("Y", 88.906);
        weights.put("Zr", 91.224);
        weights.put("Nb", 92.906);
        weights.put("Mo", 95.95);
        weights.put("Tc", 98.0);
        weights.put("Ru", 101.07);
        weights.put("Rh", 102.91);
        weights.put("Pd", 106.42);
        weights.put("Ag", 107.87);
        weights.put("Cd", 112.41);
        weights.put("In", 114.82);
        weights.put("Sn", 118.71);
        weights.put("Sb", 121.76);
        weights.put("Te", 127.60);
        weights.put("I", 126.90);
        weights.put("Xe", 131.29);
        weights.put("Cs", 132.91);
        weights.put("Ba", 137.33);
        weights.put("La", 138.91);
        weights.put("Ce", 140.12);
        weights.put("Pr", 140.91);
        weights.put("Nd", 144.24);
        weights.put("Pm", 145.0);
        weights.put("Sm", 150.36);
        weights.put("Eu", 151.96);
        weights.put("Gd", 157.25);
        weights.put("Tb", 158.93);
        weights.put("Dy", 162.50);
        weights.put("Ho", 164.93);
        weights.put("Er", 167.26);
        weights.put("Tm", 168.93);
        weights.put("Yb", 173.04);
        weights.put("Lu", 174.97);
        weights.put("Hf", 178.49);
        weights.put("Ta", 180.95);
        weights.put("W", 183.84);
        weights.put("Re", 186.21);
        weights.put("Os", 190.23);
        weights.put("Ir", 192.22);
        weights.put("Pt", 195.08);
        weights.put("Au", 196.97);
        weights.put("Hg", 200.59);
        weights.put("Tl", 204.38);
        weights.put("Pb", 207.2);
        weights.put("Bi", 208.98);
        weights.put("Po", 209.0);
        weights.put("At", 210.0);
        weights.put("Rn", 222.0);
        weights.put("Fr", 223.0);
        weights.put("Ra", 226.0);
        weights.put("Ac", 227.0);
        weights.put("Th", 232.04);
        weights.put("Pa", 231.04);
        weights.put("U", 238.03);
        weights.put("Np", 237.0);
        weights.put("Pu", 244.0);
        weights.put("Am", 243.0);
        weights.put("Cm", 247.0);
        weights.put("Bk", 247.0);
        weights.put("Cf", 251.0);
        weights.put("Es", 252.0);
        weights.put("Fm", 257.0);
        weights.put("Md", 258.0);
        weights.put("No", 259.0);
        weights.put("Lr", 262.0);
        weights.put("Rf", 267.0);
        weights.put("Db", 270.0);
        weights.put("Sg", 271.0);
        weights.put("Bh", 270.0);
        weights.put("Hs", 277.0);
        weights.put("Mt", 278.0);
        weights.put("Ds", 281.0);
        weights.put("Rg", 282.0);
        weights.put("Cn", 285.0);
        weights.put("Nh", 286.0);
        weights.put("Fl", 289.0);
        weights.put("Mc", 290.0);
        weights.put("Lv", 293.0);
        weights.put("Ts", 294.0);
        weights.put("Og", 294.0);
        return weights;
    }

    // Generate a name for a given Compound
    private static String nameCompound(Compound comp) {
        // Check for acids
        if (comp.composition.containsKey("H")) {
            // Oxyacid (contains oxygen and hydrogen)
            if (comp.composition.size() > 1 && comp.composition.containsKey("O")) {
                // Derive polyatomic anion by removing all H
                Map<String,Integer> anionComp = new HashMap<>(comp.composition);
                if (anionComp.containsKey("H")) anionComp.remove("H");
                // Construct anion formula string
                StringBuilder anionFormula = new StringBuilder();
                for (Map.Entry<String,Integer> entry : anionComp.entrySet()) {
                    String elem = entry.getKey();
                    int count = entry.getValue();
                    anionFormula.append(elem);
                    if (count > 1) anionFormula.append(count);
                }
                String anionStr = anionFormula.toString();
                if (polyIonMap.containsKey(anionStr)) {
                    String anionName = polyIonMap.get(anionStr).name;
                    if (anionName.endsWith("ate")) {
                        String base = anionName.substring(0, anionName.length()-3);
                        if (base.endsWith("sulf")) base = "sulfur";
                        if (base.endsWith("phosph")) base = "phosphor";
                        return capitalize(base + "ic acid");
                    } else if (anionName.endsWith("ite")) {
                        String base = anionName.substring(0, anionName.length()-3);
                        if (base.endsWith("sulf")) base = "sulfur";
                        if (base.endsWith("phosph")) base = "phosphor";
                        return capitalize(base + "ous acid");
                    }
                }
            }
            // Binary acid (no O)
            if (comp.composition.size() == 2) {
                String otherElem = null;
                for (String e : comp.composition.keySet()) {
                    if (!e.equals("H")) { otherElem = e; break; }
                }
                if (otherElem != null && anionNames.containsKey(otherElem)) {
                    String anionBase = anionNames.get(otherElem);
                    if (anionBase.endsWith("ide")) {
                        String base = anionBase.substring(0, anionBase.length()-3);
                        if (base.equals("sulf")) base = "sulfur";
                        if (base.equals("phosph")) base = "phosphor";
                        return "Hydro" + base + "ic acid";
                    }
                }
            }
        }
        // If ionic (metal or polyatomic):
        boolean containsMetal = false;
        for (String elem : comp.composition.keySet()) {
            if (isMetal(elem) && !elem.equals("H")) {
                containsMetal = true;
                break;
            }
        }
        if (containsMetal || comp.charge != 0) {
            // If it's a polyatomic ion by itself
            if (comp.charge != 0 && comp.composition.size() > 1) {
                String formulaStr = comp.formula;
                if (polyIonMap.containsKey(formulaStr)) {
                    return polyIonMap.get(formulaStr).name + " ion";
                }
            }
            if (comp.charge != 0 && comp.composition.size() == 1) {
                // Monatomic ion
                String elem = comp.composition.keySet().iterator().next();
                if (comp.charge > 0) {
                    String elemName = elementNames.getOrDefault(elem, elem);
                    // If multiple possible charges (transition metals), add Roman numeral
                    if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg").contains(elem)) {
                        return elemName + " (" + comp.charge + (comp.charge>0?"+":"-") + ") ion";
                    } else {
                        return elemName + " ion";
                    }
                } else {
                    if (anionNames.containsKey(elem)) {
                        return anionNames.get(elem) + " ion";
                    } else {
                        String elemName = elementNames.getOrDefault(elem, elem);
                        return elemName + " ion";
                    }
                }
            }
            // For neutral ionic compounds, split into cation and anion
            String cationName = "", anionName = "";
            // Try to find a polyatomic anion in the formula
            for (String polyForm : polyIonMap.keySet()) {
                PolyIon poly = polyIonMap.get(polyForm);
                if (poly.charge < 0) {
                    Map<String,Integer> polyComp = parseFormulaComposition(polyForm);
                    int possibleCount = Integer.MAX_VALUE;
                    for (String e : polyComp.keySet()) {
                        if (!comp.composition.containsKey(e)) {
                            possibleCount = 0;
                            break;
                        }
                        possibleCount = Math.min(possibleCount, comp.composition.get(e) / polyComp.get(e));
                    }
                    if (possibleCount > 0 && possibleCount != Integer.MAX_VALUE) {
                        // Remove that many poly groups from composition
                        Map<String,Integer> remaining = new HashMap<>(comp.composition);
                        for (String e : polyComp.keySet()) {
                            remaining.put(e, remaining.get(e) - polyComp.get(e)*possibleCount);
                            if (remaining.get(e) == 0) remaining.remove(e);
                        }
                        if (remaining.isEmpty()) {
                            // The whole compound is just that poly ion multiplied
                            return poly.name;
                        }
                        // Remaining is cation part
                        if (remaining.size() == 1) {
                            String catElem = remaining.keySet().iterator().next();
                            int catCount = remaining.get(catElem);
                            String baseName = elementNames.getOrDefault(catElem, catElem);
                            // Determine cation charge from anion charge balance
                            int totalAnionCharge = poly.charge * possibleCount;
                            int cationCharge = - totalAnionCharge / catCount;
                            if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg").contains(catElem) && cationCharge != 0) {
                                baseName += " (" + cationCharge + "+)";
                            }
                            cationName = baseName;
                            anionName = poly.name;
                            break;
                        }
                    }
                }
            }
            if (!cationName.isEmpty() && !anionName.isEmpty()) {
                return cationName + " " + anionName;
            }
            // If no polyatomic found, assume binary ionic (metal + nonmetal)
            if (comp.composition.size() == 2) {
                String metalElem = null;
                String nonmetalElem = null;
                for (String e : comp.composition.keySet()) {
                    if (isMetal(e)) metalElem = e;
                    else nonmetalElem = e;
                }
                if (metalElem != null && nonmetalElem != null) {
                    String metalName = elementNames.getOrDefault(metalElem, metalElem);
                    int metalCount = comp.composition.get(metalElem);
                    int nonCount = comp.composition.get(nonmetalElem);
                    int nonCharge = -1;
                    if (nonmetalElem.equals("O")) nonCharge = -2;
                    else if (nonmetalElem.equals("N")) nonCharge = -3;
                    else if (nonmetalElem.equals("S")) nonCharge = -2;
                    else if (nonmetalElem.equals("P")) nonCharge = -3;
                    else nonCharge = -1;
                    int metalCharge = 0;
                    if (nonCharge != -1) {
                        metalCharge = (- nonCharge * nonCount) / metalCount;
                    }
                    if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg").contains(metalElem) && metalCharge != 0) {
                        metalName += " (" + metalCharge + "+)";
                    }
                    String anName = anionNames.getOrDefault(nonmetalElem, nonmetalElem);
                    return metalName + " " + anName;
                }
            }
        }
        // If covalent (nonmetals):
        if (comp.charge == 0) {
            if (comp.composition.size() == 2) {
                Iterator<String> it = comp.composition.keySet().iterator();
                String e1 = it.next();
                String e2 = it.next();
                // Order elements so that the more electropositive comes first
                List<String> order = Arrays.asList("C","P","N","H","Si","B","S","I","Br","Cl","O","F");
                if (order.indexOf(e2) < order.indexOf(e1)) {
                    String temp = e1; e1 = e2; e2 = temp;
                }
                int c1 = comp.composition.get(e1);
                int c2 = comp.composition.get(e2);
                String name1 = elementNames.getOrDefault(e1, e1);
                String base2;
                if (anionNames.containsKey(e2)) {
                    String anName = anionNames.get(e2);
                    base2 = anName.endsWith("ide") ? anName.substring(0, anName.length()-3) : anName;
                } else {
                    base2 = elementNames.getOrDefault(e2, e2);
                }
                String prefix1 = prefixForNumber(c1);
                String prefix2 = prefixForNumber(c2);
                if (prefix1.equals("mono")) prefix1 = "";  // no "mono" for first element
                String name = "";
                if (!prefix1.isEmpty()) name += prefix1;
                name += name1.toLowerCase();
                name += " " + prefix2 + base2.toLowerCase() + "ide";
                return capitalize(name);
            }
        }
        // Fallback: just return formula
        return comp.formula;
    }

    private static String capitalize(String s) {
        if (s == null || s.isEmpty()) return s;
        return Character.toUpperCase(s.charAt(0)) + s.substring(1);
    }

    private static String prefixForNumber(int n) {
        switch(n) {
            case 1: return "mono";
            case 2: return "di";
            case 3: return "tri";
            case 4: return "tetra";
            case 5: return "penta";
            case 6: return "hexa";
            case 7: return "hepta";
            case 8: return "octa";
            case 9: return "nona";
            case 10: return "deca";
            default: return "";
        }
    }
}
