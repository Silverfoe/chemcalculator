import java.awt.*;
import java.awt.event.*;
import java.math.BigInteger;
import java.util.*;
import java.util.List;
import java.util.regex.*;
import javax.swing.*;

public class ChemistryTool {
    // --- Data Structures for Chemical Data ---
    // Atomic weights (g/mol) for elements, stored in atomic number order (1=H, 2=He, ..., 118=Og)
    private static final Map<String, Double> atomicWeights = loadAtomicWeights();
    // Element full names for a selection of elements (for naming and tooltips)
    private static final Map<String, String> elementNames = new HashMap<>();
    // Monatomic anion names (for naming binary ionic compounds)
    private static final Map<String, String> anionNames = new HashMap<>();
    // Common polyatomic ions (formula -> name and charge)
    private static class PolyIon { String name; int charge; PolyIon(String name,int charge){this.name=name; this.charge=charge;} }
    private static final Map<String, PolyIon> polyIonMap = new HashMap<>();
    static {
        // Populate elementNames for common elements
        elementNames.put("H", "Hydrogen");    elementNames.put("He", "Helium");
        elementNames.put("Li", "Lithium");    elementNames.put("Be", "Beryllium");
        elementNames.put("B", "Boron");       elementNames.put("C", "Carbon");
        elementNames.put("N", "Nitrogen");    elementNames.put("O", "Oxygen");
        elementNames.put("F", "Fluorine");    elementNames.put("Ne", "Neon");
        elementNames.put("Na", "Sodium");     elementNames.put("Mg", "Magnesium");
        elementNames.put("Al", "Aluminum");   elementNames.put("Si", "Silicon");
        elementNames.put("P", "Phosphorus");  elementNames.put("S", "Sulfur");
        elementNames.put("Cl", "Chlorine");   elementNames.put("Ar", "Argon");
        elementNames.put("K", "Potassium");   elementNames.put("Ca", "Calcium");
        elementNames.put("Sc", "Scandium");   elementNames.put("Ti", "Titanium");
        elementNames.put("V", "Vanadium");    elementNames.put("Cr", "Chromium");
        elementNames.put("Mn", "Manganese");  elementNames.put("Fe", "Iron");
        elementNames.put("Co", "Cobalt");     elementNames.put("Ni", "Nickel");
        elementNames.put("Cu", "Copper");     elementNames.put("Zn", "Zinc");
        elementNames.put("Ag", "Silver");     elementNames.put("Au", "Gold");
        elementNames.put("Hg", "Mercury");    elementNames.put("Pb", "Lead");
        elementNames.put("I", "Iodine");      elementNames.put("U", "Uranium");
        elementNames.put("Pu", "Plutonium");  // ... (etc., can add more as needed)
        // Monatomic anion names (negative ions of single elements)
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
        // Common polyatomic ions with their names and charges
        polyIonMap.put("NH4", new PolyIon("Ammonium", +1));
        polyIonMap.put("OH",  new PolyIon("Hydroxide", -1));
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
        polyIonMap.put("ClO",  new PolyIon("Hypochlorite", -1));
        polyIonMap.put("MnO4", new PolyIon("Permanganate", -1));
        polyIonMap.put("Cr2O7",new PolyIon("Dichromate", -2));
        polyIonMap.put("CrO4", new PolyIon("Chromate", -2));
        polyIonMap.put("CN",   new PolyIon("Cyanide", -1));
        polyIonMap.put("O2",   new PolyIon("Peroxide", -2));
        polyIonMap.put("C2H3O2", new PolyIon("Acetate", -1));
        polyIonMap.put("SCN",  new PolyIon("Thiocyanate", -1));
        // (Additional polyatomic ions can be added as needed)
    }

    // Compound representation (parsed formula)
    private static class Compound {
        String formula;               // Formula string (without leading coefficient)
        Map<String,Integer> composition; // Element counts
        int charge;                   // Net charge
        Compound(String formula, Map<String,Integer> comp, int charge) {
            this.formula = formula;
            this.composition = comp;
            this.charge = charge;
        }
    }

    // Utility class for exact fraction arithmetic (for algebraic balancing)
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
            // Normalize sign
            if (denominator.signum() < 0) {
                numerator = numerator.negate();
                denominator = denominator.negate();
            }
            BigInteger gcd = numerator.gcd(denominator);
            if (gcd.signum() < 0) gcd = gcd.negate();
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
        @Override
        public String toString() {
            return num + "/" + den;
        }
    }

    // HalfReaction class for balancing redox equations
    static class HalfReaction {
        Compound reactant;
        Compound product;
        List<Compound> leftExtras = new ArrayList<>();   // species added to reactant side
        List<Compound> rightExtras = new ArrayList<>();  // species added to product side
        int electrons = 0;
        boolean electronsOnLeft = false;
    }

    public static void main(String[] args) {
        // Set Nimbus Look-and-Feel for modern UI
        try {
            for (UIManager.LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (Exception e) {
            // If Nimbus not available, fall back to default
        }

        // Create main application window
        JFrame frame = new JFrame("Chemistry Tool");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1000, 650);  // a bit wider to accommodate periodic table text

        // Use a tabbed pane to organize features
        JTabbedPane tabs = new JTabbedPane();

        // --- Equation Balancer Tab ---
        JPanel balancePanel = new JPanel(new BorderLayout());
        // Input area for equation
        JPanel inputPanel1 = new JPanel(new FlowLayout());
        JLabel eqLabel = new JLabel("Enter unbalanced equation:");
        JTextField eqField = new JTextField(40);
        JButton balanceButton = new JButton("Balance");
        inputPanel1.add(eqLabel);
        inputPanel1.add(eqField);
        inputPanel1.add(balanceButton);
        balancePanel.add(inputPanel1, BorderLayout.NORTH);
        // Output area for steps and balanced equation
        JTextArea balanceOutput = new JTextArea();
        balanceOutput.setEditable(false);
        balanceOutput.setLineWrap(true);
        balanceOutput.setWrapStyleWord(true);
        JScrollPane scrollBalance = new JScrollPane(balanceOutput);
        balancePanel.add(scrollBalance, BorderLayout.CENTER);
        tabs.addTab("Equation Balancer", balancePanel);

        // --- Gram Formula Mass (GFM) Calculator Tab ---
        JPanel gfmPanel = new JPanel(new BorderLayout());
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
        JScrollPane scrollGFM = new JScrollPane(gfmOutput);
        gfmPanel.add(scrollGFM, BorderLayout.CENTER);
        tabs.addTab("GFM Calculator", gfmPanel);

        // --- Periodic Table Viewer Tab ---
        JPanel tablePanel = new JPanel(new BorderLayout());
        // Create panel for the main periodic table grid
        JPanel mainTablePanel = new JPanel(new GridLayout(7, 18, 2, 2));  // 7 periods x 18 groups
        // Define colors for element categories
        Color alkaliColor       = new Color(255, 182, 193);  // light pink (Group 1 except H)
        Color alkalineColor     = new Color(255, 228, 181);  // light orange (Group 2)
        Color transitionColor   = new Color(255, 255, 153);  // light yellow (d-block metals)
        Color postTransitionColor = new Color(176, 224, 230); // powder blue (post-transition metals)
        Color metalloidColor    = Color.LIGHT_GRAY;         // metalloids (stair-step line elements)
        Color nonmetalColor     = new Color(255, 255, 224);  // light yellow (other nonmetals + H)
        Color halogenColor      = new Color(144, 238, 144);  // light green (halogens)
        Color nobleGasColor     = new Color(224, 255, 255);  // light cyan (noble gases)
        Color lanthColor        = new Color(221, 160, 221);  // plum (lanthanides)
        Color actColor          = new Color(216, 191, 216);  // thistle (actinides)

        // Tooltip and detail helper: on hover or focus, show element details
        final ToolTipManager ttm = ToolTipManager.sharedInstance();
        ttm.setInitialDelay(100);  // show quickly
        ttm.setReshowDelay(50);

        // Build the main periodic table grid
        for (int period = 1; period <= 7; period++) {
            for (int group = 1; group <= 18; group++) {
                JButton cell = new JButton();
                cell.setMargin(new Insets(1,1,1,1));
                cell.setFocusPainted(false);
                cell.setOpaque(true);
                cell.setForeground(Color.BLACK);
                cell.setFont(new Font("SansSerif", Font.BOLD, 11));
                String symbol = "";
                String tooltip = null;
                Color bgColor = null;
                // Determine which element (if any) goes in this period/group
                int atomicNum = 0;
                if (period == 1) {
                    if (group == 1) atomicNum = 1;       // H
                    else if (group == 18) atomicNum = 2; // He
                } else if (period == 2) {
                    if (group == 1) atomicNum = 3;   // Li
                    else if (group == 2) atomicNum = 4;   // Be
                    else if (group >= 13) atomicNum = group - 8;  // B(5) at 13 ... Ne(10) at 18
                } else if (period == 3) {
                    if (group == 1) atomicNum = 11;  // Na
                    else if (group == 2) atomicNum = 12;  // Mg
                    else if (group >= 13) atomicNum = group;      // Al(13) ... Ar(18)
                } else if (period == 4) {
                    atomicNum = 18 + group;  // K(19) at 1 ... Kr(36) at 18
                } else if (period == 5) {
                    atomicNum = 36 + group;  // Rb(37) ... Xe(54)
                } else if (period == 6) {
                    if (group == 1) atomicNum = 55;       // Cs
                    else if (group == 2) atomicNum = 56;  // Ba
                    else if (group == 3) {
                        // Lanthanide placeholder
                        symbol = "La-Lu";
                        bgColor = lanthColor;
                    } else if (group >= 4) {
                        atomicNum = 68 + group;  // Hf(72) at 4 ... Rn(86) at 18  (offset 68)
                    }
                } else if (period == 7) {
                    if (group == 1) atomicNum = 87;       // Fr
                    else if (group == 2) atomicNum = 88;  // Ra
                    else if (group == 3) {
                        // Actinide placeholder
                        symbol = "Ac-Lr";
                        bgColor = actColor;
                    } else if (group >= 4) {
                        atomicNum = 100 + group; // Rf(104) at 4 ... Og(118) at 18 (offset 100)
                    }
                }
                if (atomicNum > 0 && atomicNum <= 118) {
                    // There is an element in this cell
                    // Get element symbol by atomic number from atomicWeights map (which is in order)
                    // We need to fetch by index in insertion order. Let's pre-build an array of symbols for quick access.
                    // (Alternatively, maintain a static list of symbols in order of atomic number.)
                    // Build static list of element symbols if not already done:
                    if (ELEMENT_SYMBOLS == null) {
                        ELEMENT_SYMBOLS = atomicWeights.keySet().toArray(new String[0]);
                    }
                    symbol = ELEMENT_SYMBOLS[atomicNum - 1];  // atomicNum is 1-indexed in list
                    // Determine category and set color
                    if (symbol.equals("H")) {
                        // Hydrogen (treated as nonmetal)
                        bgColor = nonmetalColor;
                    } else if (isElementIn(symbol, "Li","Na","K","Rb","Cs","Fr")) {
                        bgColor = alkaliColor;
                    } else if (isElementIn(symbol, "Be","Mg","Ca","Sr","Ba","Ra")) {
                        bgColor = alkalineColor;
                    } else if (atomicNum >= 21 && atomicNum <= 30 || atomicNum >= 39 && atomicNum <= 48 ||
                               atomicNum >= 72 && atomicNum <= 80 || atomicNum >= 104 && atomicNum <= 112) {
                        // Transition metals: Sc-Zn, Y-Cd, Hf-Hg, Rf-Cn
                        bgColor = transitionColor;
                    } else if ((isElementIn(symbol,"Al","Ga","In","Tl","Nh")) ||
                               (isElementIn(symbol,"Sn","Pb","Fl")) ||
                               (isElementIn(symbol,"Bi","Mc")) ||
                               (isElementIn(symbol,"Po","Lv"))) {
                        // Post-transition metals (a.k.a. poor metals)
                        bgColor = postTransitionColor;
                    } else if (isElementIn(symbol, "B","Si","Ge","As","Sb","Te")) {
                        // Metalloids (staircase elements)
                        bgColor = metalloidColor;
                    } else if (isElementIn(symbol, "F","Cl","Br","I","At","Ts")) {
                        // Halogens
                        bgColor = halogenColor;
                    } else if (isElementIn(symbol, "He","Ne","Ar","Kr","Xe","Rn","Og")) {
                        // Noble gases
                        bgColor = nobleGasColor;
                    } else if (atomicNum >= 57 && atomicNum <= 71) {
                        // Lanthanides (should not appear individually in main grid in this design)
                        bgColor = lanthColor;
                    } else if (atomicNum >= 89 && atomicNum <= 103) {
                        // Actinides (should not appear individually in main grid)
                        bgColor = actColor;
                    } else {
                        // Other nonmetals (like C, N, O, P, S, Se)
                        bgColor = nonmetalColor;
                    }
                    // Prepare tooltip with element details
                    String name = elementNames.getOrDefault(symbol, "Element " + atomicNum);
                    Double mass = atomicWeights.get(symbol);
                    tooltip = "<html><b>" + symbol + " - " + name + "</b><br>" +
                              "Atomic Number: " + atomicNum + "<br>" +
                              "Atomic Mass: " + (mass != null ? mass : "?") + " g/mol" + "</html>";
                }
                if (!symbol.isEmpty()) {
                    cell.setText(symbol);
                }
                if (bgColor != null) {
                    cell.setBackground(bgColor);
                } else {
                    cell.setBackground(tablePanel.getBackground()); // default background for empties
                    cell.setEnabled(false);
                }
                if (tooltip != null) {
                    cell.setToolTipText(tooltip);
                }
                mainTablePanel.add(cell);
            }
        }
        // Create panel for Lanthanide and Actinide series below main table
        JPanel lanthActPanel = new JPanel();
        lanthActPanel.setLayout(new BoxLayout(lanthActPanel, BoxLayout.Y_AXIS));
        // Lanthanides row
        JPanel lanthRowPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        lanthRowPanel.add(new JLabel("Lanthanides: "));
        for (int Z = 57; Z <= 71; Z++) {
            String sym = ELEMENT_SYMBOLS[Z - 1];
            JButton elemBtn = new JButton(sym);
            elemBtn.setMargin(new Insets(1,2,1,2));
            elemBtn.setOpaque(true);
            elemBtn.setBackground(lanthColor);
            elemBtn.setForeground(Color.BLACK);
            elemBtn.setFocusPainted(false);
            elemBtn.setFont(new Font("SansSerif", Font.BOLD, 11));
            String name = elementNames.getOrDefault(sym, "");
            Double mass = atomicWeights.get(sym);
            elemBtn.setToolTipText("<html><b>" + sym + " - " + name + "</b><br>" +
                                    "Atomic Number: " + Z + "<br>" +
                                    "Atomic Mass: " + (mass != null ? mass : "?") + " g/mol" + "</html>");
            lanthRowPanel.add(elemBtn);
        }
        // Actinides row
        JPanel actRowPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        actRowPanel.add(new JLabel("Actinides: "));
        for (int Z = 89; Z <= 103; Z++) {
            String sym = ELEMENT_SYMBOLS[Z - 1];
            JButton elemBtn = new JButton(sym);
            elemBtn.setMargin(new Insets(1,2,1,2));
            elemBtn.setOpaque(true);
            elemBtn.setBackground(actColor);
            elemBtn.setForeground(Color.BLACK);
            elemBtn.setFocusPainted(false);
            elemBtn.setFont(new Font("SansSerif", Font.BOLD, 11));
            String name = elementNames.getOrDefault(sym, "");
            Double mass = atomicWeights.get(sym);
            elemBtn.setToolTipText("<html><b>" + sym + " - " + name + "</b><br>" +
                                    "Atomic Number: " + Z + "<br>" +
                                    "Atomic Mass: " + (mass != null ? mass : "?") + " g/mol" + "</html>");
            actRowPanel.add(elemBtn);
        }
        // Add the sub-panels to lanthActPanel
        lanthActPanel.add(lanthRowPanel);
        lanthActPanel.add(actRowPanel);
        // Combine main table and lanth/act panels
        tablePanel.add(mainTablePanel, BorderLayout.CENTER);
        tablePanel.add(lanthActPanel, BorderLayout.SOUTH);
        tabs.addTab("Periodic Table", tablePanel);

        // --- Compound Naming Tab ---
        JPanel namePanel = new JPanel(new BorderLayout());
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
        JScrollPane scrollName = new JScrollPane(nameOutput);
        namePanel.add(scrollName, BorderLayout.CENTER);
        tabs.addTab("Compound Naming", namePanel);

        // Add tabs to frame and show
        frame.add(tabs);
        frame.setLocationRelativeTo(null);  // center on screen
        frame.setVisible(true);

        // --- Action Listeners for Buttons ---
        // Balance equation button:
        balanceButton.addActionListener((ActionEvent e) -> {
            String eq = eqField.getText().trim();
            if (eq.isEmpty()) return;
            List<String> steps = balanceChemicalEquation(eq);
            // Display each step on a new line
            balanceOutput.setText("");
            for (String step : steps) {
                balanceOutput.append(step + "\n");
            }
        });
        // GFM calculate button:
        gfmButton.addActionListener((ActionEvent e) -> {
            String input = gfmField.getText().trim();
            if (input.isEmpty()) return;
            StringBuilder result = new StringBuilder();
            String[] formulas = input.split("\\s*,\\s*");
            for (String formula : formulas) {
                if (formula.isEmpty()) continue;
                result.append("Formula: ").append(formula).append("\n");
                try {
                    Map<String,Integer> comp = parseFormulaComposition(formula);
                    double totalMass = 0.0;
                    for (Map.Entry<String,Integer> entry : comp.entrySet()) {
                        String element = entry.getKey();
                        int count = entry.getValue();
                        double atomicMass = atomicWeights.getOrDefault(element, 0.0);
                        double mass = atomicMass * count;
                        totalMass += mass;
                        result.append(String.format("  %s: %.3f g/mol (x%d)\n", element, mass, count));
                    }
                    result.append(String.format("  Total GFM: %.3f g/mol\n\n", totalMass));
                } catch (Exception ex) {
                    result.append("  Error parsing formula.\n\n");
                }
            }
            gfmOutput.setText(result.toString());
        });
        // Name compound button:
        nameButton.addActionListener((ActionEvent e) -> {
            String formulaInput = nameField.getText().trim();
            if (formulaInput.isEmpty()) return;
            try {
                Compound comp = parseCompound(formulaInput);
                String name = nameCompound(comp);
                nameOutput.setText(name);
            } catch (Exception ex) {
                nameOutput.setText("Unable to name the compound. Please check the formula.");
            }
        });
    }

    // --- Chemical Equation Balancing Logic ---
    private static List<String> balanceChemicalEquation(String equation) {
        List<String> steps = new ArrayList<>();
        // Split into reactant and product part
        String[] sides = equation.replaceAll("<->", "->").split("->|=");
        if (sides.length != 2) {
            steps.add("Error: Equation must have a single '->' (or '=') separating reactants and products.");
            return steps;
        }
        String reactantsStr = sides[0].trim();
        String productsStr = sides[1].trim();
        if (reactantsStr.isEmpty() || productsStr.isEmpty()) {
            steps.add("Error: Reactant or product side is empty.");
            return steps;
        }
        // Split compounds by '+' and parse each compound
        String[] reactTokens = reactantsStr.split("\\s*\\+\\s*");
        String[] prodTokens = productsStr.split("\\s*\\+\\s*");
        List<Compound> reactants = new ArrayList<>();
        List<Compound> products = new ArrayList<>();
        try {
            for (String token : reactTokens) {
                if (token.isBlank()) continue;
                reactants.add(parseCompound(token.trim()));
            }
            for (String token : prodTokens) {
                if (token.isBlank()) continue;
                products.add(parseCompound(token.trim()));
            }
        } catch (Exception ex) {
            steps.add("Error: Failed to parse the equation components.");
            return steps;
        }
        // Identify if redox by checking changes in oxidation numbers
        boolean isRedox = false;
        Map<String, Integer> startOx = new HashMap<>();
        Map<String, Integer> endOx = new HashMap<>();
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
            int ox1 = startOx.get(elem);
            int ox2 = endOx.getOrDefault(elem, ox1);
            if (ox2 > ox1) {
                oxidizedElements.add(elem);
                isRedox = true;
            } else if (ox2 < ox1) {
                reducedElements.add(elem);
                isRedox = true;
            }
        }

        if (isRedox && !oxidizedElements.isEmpty() && !reducedElements.isEmpty()) {
            // Attempt half-reaction balancing
            steps.add("Redox reaction detected. Using half-reaction method:");
            List<HalfReaction> halfReactions = new ArrayList<>();
            // If one element is both oxidized and reduced (disproportionation), handle specially
            if (oxidizedElements.size() == 1 && reducedElements.size() == 1 && oxidizedElements.get(0).equals(reducedElements.get(0))) {
                String elem = oxidizedElements.get(0);
                // Find one reactant compound containing elem and two product compounds containing elem
                Compound reactantComp = null;
                for (Compound comp : reactants) {
                    if (comp.composition.containsKey(elem)) { reactantComp = comp; break; }
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
                    // Determine which product is oxidation vs reduction
                    if (oxP1 > oxReact && oxP2 < oxReact) {
                        half1.product = prod1; // oxidation
                        half2.product = prod2; // reduction
                    } else if (oxP2 > oxReact && oxP1 < oxReact) {
                        half1.product = prod2;
                        half2.product = prod1;
                    } else {
                        // If uncertain, assign arbitrarily
                        half1.product = prod1;
                        half2.product = prod2;
                    }
                    halfReactions.add(half1);
                    halfReactions.add(half2);
                }
            } else {
                // Typical case: at least one oxidized and one reduced element
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
            // Balance each half reaction
            // Determine if medium is acidic or basic (presence of OH- suggests basic)
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
                // Identify key element that changes (present in both reactant and product)
                String keyElem = null;
                for (String e : half.reactant.composition.keySet()) {
                    if (half.product.composition.containsKey(e)) {
                        keyElem = e;
                        break;
                    }
                }
                if (keyElem == null) {
                    // if no common element, pick any from reactant as key
                    keyElem = half.reactant.composition.keySet().iterator().next();
                }
                steps.add("Half-reaction: " + half.reactant.formula + " -> " + half.product.formula);
                // Balance the key element by adjusting molecule count (conceptually)
                int rCount = half.reactant.composition.getOrDefault(keyElem, 0);
                int pCount = half.product.composition.getOrDefault(keyElem, 0);
                if (rCount != pCount && rCount > 0 && pCount > 0) {
                    // Show step, actual scalar multiplication will be handled when combining halves
                    steps.add("Balance " + keyElem + ": " + half.reactant.formula + " -> " + half.product.formula);
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
                // Balance hydrogen by adding H+ (acidic) or H2O + OH- (basic)
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
                    } else { // basic solution
                        if (hydLeft < hydRight) {
                            int diff = hydRight - hydLeft;
                            // In basic medium, add H2O to left (to supply H) and OH- to right for each needed H
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
                // Balance charge by adding electrons (e-)
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
            // If two half-reactions (common case), combine them
            if (halfReactions.size() == 2) {
                HalfReaction h1 = halfReactions.get(0);
                HalfReaction h2 = halfReactions.get(1);
                int e1 = h1.electrons;
                int e2 = h2.electrons;
                if (e1 != 0 && e2 != 0 && e1 != e2) {
                    int lcm = lcm(e1, e2);
                    int factor1 = lcm / e1;
                    int factor2 = lcm / e2;
                    steps.add("Multiply half-reactions to equalize electrons: oxidation x" + factor1 + ", reduction x" + factor2);
                    // Multiply species counts accordingly (conceptually combine later)
                    // (We won't explicitly multiply comp counts here, just note it in steps)
                }
                // Combine half-reactions into full reaction
                Map<String,Integer> leftMap = new LinkedHashMap<>();
                Map<String,Integer> rightMap = new LinkedHashMap<>();
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
                        if (half.electronsOnLeft) addSpecies(leftMap, eSymbol, half.electrons);
                        else                      addSpecies(rightMap, eSymbol, half.electrons);
                    }
                }
                // Cancel out electrons and other species that appear on both sides
                cancelSpecies(leftMap, rightMap, "e-");
                cancelSpecies(leftMap, rightMap, "H2O");
                cancelSpecies(leftMap, rightMap, "H+");
                cancelSpecies(leftMap, rightMap, "OH-");
                // Format the final balanced equation
                String finalEq = formatEquation(leftMap, rightMap);
                steps.add("Balanced Equation: " + finalEq);
            } else {
                // If not exactly two half-reactions, fall back to algebraic method
                steps.addAll(balanceNonRedoxReaction(reactants, products));
            }
        } else {
            // Not a redox reaction, or could not identify redox changes; use algebraic balancing
            steps.addAll(balanceNonRedoxReaction(reactants, products));
        }
        return steps;
    }

    // Balance a reaction via linear algebra (if not using redox half-reaction method)
    private static List<String> balanceNonRedoxReaction(List<Compound> reactants, List<Compound> products) {
        List<String> steps = new ArrayList<>();
        steps.add("Using algebraic method for balancing:");
        // Combine reactants and products into one list for indexing
        List<Compound> allCompounds = new ArrayList<>();
        allCompounds.addAll(reactants);
        allCompounds.addAll(products);
        int n = allCompounds.size();
        int reactantCount = reactants.size();
        // Collect unique elements involved
        Set<String> elementSet = new HashSet<>();
        for (Compound comp : allCompounds) {
            elementSet.addAll(comp.composition.keySet());
        }
        List<String> elements = new ArrayList<>(elementSet);
        int m = elements.size();
        // Set up matrix equation Ax = 0 (with one coefficient fixed to 1)
        Fraction[][] A = new Fraction[m][n - 1];  // coefficients matrix (exclude first compound)
        Fraction[] B = new Fraction[m];          // constants (from first compound's contribution)
        // We fix the first compound’s coefficient = 1 and move its contribution to constants
        Compound first = allCompounds.get(0);
        for (int i = 0; i < m; i++) {
            String elem = elements.get(i);
            int contribFirst = first.composition.getOrDefault(elem, 0);
            if (!reactants.contains(first)) {
                // If first compound is a product, treat its contributions as negative
                contribFirst = -contribFirst;
            }
            B[i] = new Fraction(-contribFirst, 1);
            // Fill matrix coefficients for compounds 2..n
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
        // Solve A * x = B using Gaussian elimination (on fractions)
        int unknowns = n - 1;
        int row = 0, col = 0;
        int[] pivotCol = new int[m];
        Arrays.fill(pivotCol, -1);
        while (row < m && col < unknowns) {
            // Find pivot row for this column
            int pivot = row;
            while (pivot < m && A[pivot][col].num.equals(BigInteger.ZERO)) {
                pivot++;
            }
            if (pivot == m) {
                col++;
                continue;
            }
            if (pivot != row) {
                // Swap rows
                Fraction[] tempRow = A[row];
                A[row] = A[pivot];
                A[pivot] = tempRow;
                Fraction tempB = B[row];
                B[row] = B[pivot];
                B[pivot] = tempB;
            }
            pivotCol[row] = col;
            // Normalize pivot row
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
        // Back substitution
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
        // Include the first compound's coefficient as 1 (as a Fraction)
        Fraction firstCoef = new Fraction(1);
        // Scale all coefficients by LCM of denominators to get integer coefficients
        BigInteger lcmDen = BigInteger.ONE;
        for (Fraction frac : solution) {
            lcmDen = lcmDen.multiply(frac.den).divide(lcmDen.gcd(frac.den));
        }
        BigInteger scale = lcmDen;
        int[] intCoeffs = new int[n];
        intCoeffs[0] = firstCoef.num.multiply(scale).divide(firstCoef.den).intValue();
        for (int j = 1; j < n; j++) {
            Fraction frac = solution[j-1];
            BigInteger num = frac.num.multiply(scale).divide(frac.den);
            intCoeffs[j] = num.intValue();
        }
        // Normalize coefficients (find gcd and divide to get smallest integers)
        int gcdCoeffs = 0;
        for (int coef : intCoeffs) {
            gcdCoeffs = gcdCoeffs == 0 ? coef : gcd(gcdCoeffs, coef);
        }
        if (gcdCoeffs > 0) {
            for (int j = 0; j < n; j++) {
                intCoeffs[j] /= gcdCoeffs;
            }
        }
        // Build formatted balanced equation string
        StringBuilder sb = new StringBuilder();
        // Reactants
        for (int j = 0; j < reactantCount; j++) {
            int coef = intCoeffs[j];
            Compound comp = reactants.get(j);
            if (j > 0) sb.append(" + ");
            if (coef != 1) sb.append(coef).append(" ");
            sb.append(comp.formula);
        }
        sb.append(" -> ");
        // Products
        for (int j = reactantCount; j < n; j++) {
            int coef = intCoeffs[j];
            Compound comp = products.get(j - reactantCount);
            if (j > reactantCount) sb.append(" + ");
            if (coef != 1) sb.append(coef).append(" ");
            sb.append(comp.formula);
        }
        steps.add("Balanced Equation: " + sb.toString());
        return steps;
    }

    // --- Compound Naming Logic ---
    private static String nameCompound(Compound comp) {
        // Handle acids first
        if (comp.composition.containsKey("H")) {
            // Oxyacid (contains H and O)
            if (comp.composition.size() > 1 && comp.composition.containsKey("O")) {
                // Form anion by removing H
                Map<String,Integer> anionComp = new HashMap<>(comp.composition);
                anionComp.remove("H");
                StringBuilder anionFormula = new StringBuilder();
                for (Map.Entry<String,Integer> entry : anionComp.entrySet()) {
                    anionFormula.append(entry.getKey());
                    int count = entry.getValue();
                    if (count > 1) anionFormula.append(count);
                }
                String anionStr = anionFormula.toString();
                if (polyIonMap.containsKey(anionStr)) {
                    String anionName = polyIonMap.get(anionStr).name;
                    if (anionName.endsWith("ate")) {
                        // e.g. sulfate -> sulfuric acid
                        String base = anionName.substring(0, anionName.length()-3);
                        if (base.endsWith("sulf")) base = "sulfur";
                        if (base.endsWith("phosph")) base = "phosphor";
                        return capitalize(base + "ic acid");
                    } else if (anionName.endsWith("ite")) {
                        // e.g. sulfite -> sulfurous acid
                        String base = anionName.substring(0, anionName.length()-3);
                        if (base.endsWith("sulf")) base = "sulfur";
                        if (base.endsWith("phosph")) base = "phosphor";
                        return capitalize(base + "ous acid");
                    }
                }
            }
            // Binary acid (no oxygen, e.g. HCl)
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
        // If ionic (contains metal or polyatomic or overall charge)
        boolean containsMetal = false;
        for (String elem : comp.composition.keySet()) {
            if (isMetal(elem) && !elem.equals("H")) {
                containsMetal = true;
                break;
            }
        }
        if (containsMetal || comp.charge != 0) {
            // If compound itself is a polyatomic ion (more than one element and has charge)
            if (comp.charge != 0 && comp.composition.size() > 1) {
                if (polyIonMap.containsKey(comp.formula)) {
                    return polyIonMap.get(comp.formula).name + " ion";
                }
            }
            // Monatomic ion
            if (comp.charge != 0 && comp.composition.size() == 1) {
                String elem = comp.composition.keySet().iterator().next();
                if (comp.charge > 0) {
                    // Cation (positive ion)
                    String elemName = elementNames.getOrDefault(elem, elem);
                    // If it's a transition metal with multiple possible charges, include Roman numeral
                    if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg","Cr","Mn").contains(elem) && comp.charge != 0) {
                        return elemName + " (" + comp.charge + (comp.charge > 0 ? "+" : "-") + ") ion";
                    } else {
                        return elemName + " ion";
                    }
                } else {
                    // Anion (negative ion)
                    if (anionNames.containsKey(elem)) {
                        return anionNames.get(elem) + " ion";
                    } else {
                        String elemName = elementNames.getOrDefault(elem, elem);
                        return elemName + " ion";
                    }
                }
            }
            // Neutral ionic compound (metal + nonmetal or polyatomic)
            String cationName = "", anionName = "";
            // Check if the formula contains a polyatomic anion
            for (String polyForm : polyIonMap.keySet()) {
                PolyIon poly = polyIonMap.get(polyForm);
                if (poly.charge < 0) {
                    // See if polyForm appears in comp (in stoichiometric proportion)
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
                        // Found that comp includes polyForm * possibleCount
                        // Remove that many poly groups from a copy of comp to isolate cation
                        Map<String,Integer> remaining = new HashMap<>(comp.composition);
                        for (String e : polyComp.keySet()) {
                            remaining.put(e, remaining.get(e) - polyComp.get(e) * possibleCount);
                            if (remaining.get(e) == 0) remaining.remove(e);
                        }
                        if (remaining.isEmpty()) {
                            // The entire compound is just the polyatomic ion repeated
                            return poly.name; // e.g. "Sulfate" for SO4
                        }
                        // Remaining part is the cation
                        if (remaining.size() == 1) {
                            String catElem = remaining.keySet().iterator().next();
                            int catCount = remaining.get(catElem);
                            String baseName = elementNames.getOrDefault(catElem, catElem);
                            // Determine cation charge by charge balance
                            int totalAnionCharge = poly.charge * possibleCount;
                            // Total positive charge should balance negative: charge * count * polyCount + cationCharge*catCount = 0
                            int cationCharge = - totalAnionCharge / catCount;
                            if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg","Cr","Mn").contains(catElem) && cationCharge != 0) {
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
            // If no polyatomic anion found, assume binary ionic (metal + nonmetal)
            if (comp.composition.size() == 2) {
                String metalElem = null, nonmetalElem = null;
                for (String e : comp.composition.keySet()) {
                    if (isMetal(e)) metalElem = e;
                    else nonmetalElem = e;
                }
                if (metalElem != null && nonmetalElem != null) {
                    String metalName = elementNames.getOrDefault(metalElem, metalElem);
                    int metalCount = comp.composition.get(metalElem);
                    int nonCount = comp.composition.get(nonmetalElem);
                    // Determine nonmetal's typical charge
                    int nonCharge;
                    if      (nonmetalElem.equals("O")) nonCharge = -2;
                    else if (nonmetalElem.equals("N")) nonCharge = -3;
                    else if (nonmetalElem.equals("S")) nonCharge = -2;
                    else if (nonmetalElem.equals("P")) nonCharge = -3;
                    else nonCharge = -1;
                    // Calculate metal charge needed to balance
                    int metalCharge = 0;
                    if (nonCharge != 0) {
                        metalCharge = (- nonCharge * nonCount) / metalCount;
                    }
                    if (Arrays.asList("Fe","Cu","Co","Sn","Pb","Hg","Cr","Mn").contains(metalElem) && metalCharge != 0) {
                        metalName += " (" + metalCharge + "+)";
                    }
                    String baseAnion = anionNames.getOrDefault(nonmetalElem, nonmetalElem);
                    return metalName + " " + baseAnion;
                }
            }
            // If reached here, fallback to formula as name (unlikely)
        }
        // If covalent (all nonmetals, no overall charge)
        if (comp.charge == 0) {
            // Use Greek prefixes for number of atoms (binary covalent compounds)
            if (comp.composition.size() == 2) {
                Iterator<String> it = comp.composition.keySet().iterator();
                String e1 = it.next();
                String e2 = it.next();
                // Order elements: more electropositive (lower group number, except that halogens and others: use a common order list)
                List<String> covOrder = Arrays.asList("C","P","N","H","Si","B","S","I","Br","Cl","O","F");
                if (covOrder.indexOf(e2) < covOrder.indexOf(e1)) {
                    String temp = e1; e1 = e2; e2 = temp;
                }
                int c1 = comp.composition.get(e1);
                int c2 = comp.composition.get(e2);
                String name1 = elementNames.getOrDefault(e1, e1);
                String base2;
                if (anionNames.containsKey(e2)) {
                    String anName = anionNames.get(e2);  // e.g. "Oxide"
                    base2 = anName.endsWith("ide") ? anName.substring(0, anName.length()-3) : anName;
                } else {
                    base2 = elementNames.getOrDefault(e2, e2);
                }
                String prefix1 = prefixForNumber(c1);
                String prefix2 = prefixForNumber(c2);
                if (prefix1.equals("mono")) prefix1 = "";  // no "mono" for first element
                String compoundName = "";
                if (!prefix1.isEmpty()) compoundName += prefix1;
                compoundName += name1.toLowerCase();
                compoundName += " " + prefix2 + base2.toLowerCase() + "ide";
                return capitalize(compoundName);
            }
        }
        // Fallback: return formula as name if no rules matched
        return comp.formula;
    }

    private static String capitalize(String s) {
        if (s == null || s.isEmpty()) return s;
        return Character.toUpperCase(s.charAt(0)) + s.substring(1);
    }
    private static String prefixForNumber(int n) {
        return switch (n) {
            case 1 -> "mono";
            case 2 -> "di";
            case 3 -> "tri";
            case 4 -> "tetra";
            case 5 -> "penta";
            case 6 -> "hexa";
            case 7 -> "hepta";
            case 8 -> "octa";
            case 9 -> "nona";
            case 10 -> "deca";
            default -> "";
        };
    }

    // --- Formula Parsing ---
    private static Compound parseCompound(String formulaStr) {
        String formula = formulaStr.trim();
        // Remove any leading coefficient (e.g., "2 H2O" -> "H2O")
        if (!formula.isEmpty() && Character.isDigit(formula.charAt(0))) {
            int idx = 0;
            while (idx < formula.length() && Character.isDigit(formula.charAt(idx))) {
                idx++;
            }
            if (idx < formula.length() && formula.charAt(idx) == ' ') idx++;
            formula = formula.substring(idx);
        }
        // Extract and parse charge if present (e.g., SO4^2- or SO4(2-) or SO4-)
        int charge = 0;
        Pattern chargePattern = Pattern.compile("(.*?)(?:\\^(\\d+)?([+-])|\\((\\d+)([+-])\\)|([+-]))$");
        Matcher matcher = chargePattern.matcher(formula);
        String coreFormula = formula;
        if (matcher.matches()) {
            coreFormula = matcher.group(1);
            String number = null;
            String sign = null;
            if (matcher.group(3) != null) {
                number = matcher.group(2);
                sign = matcher.group(3);
            } else if (matcher.group(5) != null) {
                number = matcher.group(4);
                sign = matcher.group(5);
            } else if (matcher.group(6) != null) {
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

    // Parse a formula into element composition (supports nested parentheses)
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
                // Find matching closing parenthesis
                int depth = 1;
                int j = i + 1;
                while (j < formula.length() && depth > 0) {
                    if (formula.charAt(j) == '(') depth++;
                    if (formula.charAt(j) == ')') depth--;
                    j++;
                }
                if (depth != 0) throw new IllegalArgumentException("Unmatched parentheses in formula");
                String subformula = formula.substring(i+1, j-1);
                // Check for numeric multiplier after ')'
                int k = j;
                StringBuilder numBuilder = new StringBuilder();
                while (k < formula.length() && Character.isDigit(formula.charAt(k))) {
                    numBuilder.append(formula.charAt(k));
                    k++;
                }
                int count = numBuilder.length() > 0 ? Integer.parseInt(numBuilder.toString()) : 1;
                // Recurse into the parenthesized subformula
                parseFormulaRecursive(subformula, multiplier * count, comp);
                i = k;
            } else if (Character.isUpperCase(ch)) {
                // Parse element symbol starting with uppercase
                StringBuilder elem = new StringBuilder();
                elem.append(ch);
                i++;
                // Append any lowercase letters (for symbols like Fe, Na, etc.)
                while (i < formula.length() && Character.isLowerCase(formula.charAt(i))) {
                    elem.append(formula.charAt(i));
                    i++;
                }
                // Parse any numeric subscript after the element
                StringBuilder numBuilder = new StringBuilder();
                while (i < formula.length() && Character.isDigit(formula.charAt(i))) {
                    numBuilder.append(formula.charAt(i));
                    i++;
                }
                int count = numBuilder.length() > 0 ? Integer.parseInt(numBuilder.toString()) : 1;
                String element = elem.toString();
                comp.put(element, comp.getOrDefault(element, 0) + count * multiplier);
            } else {
                // Skip unexpected characters (such as spaces or stray plus signs)
                i++;
            }
        }
    }

    // Compute oxidation numbers for each element in a compound (heuristic rules)
    private static Map<String,Integer> assignOxidationNumbers(Compound comp) {
        Map<String,Integer> oxMap = new HashMap<>();
        Map<String,Integer> compMap = comp.composition;
        int netCharge = comp.charge;
        if (compMap.isEmpty()) return oxMap;
        if (compMap.size() == 1) {
            // Single element compound (element in elemental form or monatomic ion)
            String elem = compMap.keySet().iterator().next();
            if (netCharge == 0) {
                // Elemental form (e.g., O2, Fe)
                oxMap.put(elem, 0);
            } else {
                // Monatomic ion: oxidation = charge / count
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
                // Hydrogen: +1 except when with metals (then -1)
                boolean withMetal = false;
                for (String e2 : compMap.keySet()) {
                    if (!e2.equals("H") && isMetal(e2)) { withMetal = true; break; }
                }
                ox = withMetal ? -1 : +1;
                oxMap.put(elem, ox);
                sumKnown += ox * count;
            } else if (elem.equals("O")) {
                // Oxygen: assume -2 (adjust later for peroxides)
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
                // Elements typically forming negative ions (halogens, chalcogens, etc.)
                if ((elem.equals("Cl")||elem.equals("Br")||elem.equals("I")) && compMap.containsKey("O")) {
                    // Halogen present with oxygen (likely positive oxidation state) -> leave as unknown
                    unknowns.add(elem);
                } else if ((elem.equals("S")||elem.equals("Se")||elem.equals("Te")) && compMap.containsKey("O")) {
                    // Sulfur/selenium/tellurium in an oxyanion -> unknown (varied oxidation states)
                    unknowns.add(elem);
                } else if ((elem.equals("N")||elem.equals("P")) && compMap.containsKey("O")) {
                    // Nitrogen/phosphorus in presence of oxygen -> unknown
                    unknowns.add(elem);
                } else {
                    // Otherwise assign typical negative oxidation state
                    if (elem.equals("Cl")||elem.equals("Br")||elem.equals("I")) ox = -1;
                    else if (elem.equals("S")||elem.equals("Se")||elem.equals("Te")) ox = -2;
                    else if (elem.equals("N")||elem.equals("P")) ox = -3;
                    else ox = -1;
                    oxMap.put(elem, ox);
                    sumKnown += ox * count;
                }
            } else {
                // Unknown (likely transition metal or element with variable states)
                unknowns.add(elem);
            }
        }
        if (!unknowns.isEmpty()) {
            if (unknowns.size() == 1) {
                String elem = unknowns.get(0);
                int count = compMap.get(elem);
                // Assign whatever oxidation number balances net charge
                int ox = (netCharge - sumKnown) / count;
                oxMap.put(elem, ox);
            } else {
                // Multiple unknown oxidation states (e.g., organic molecules) -> assign 0 as placeholder
                for (String elem : unknowns) {
                    oxMap.put(elem, 0);
                }
            }
        }
        // Adjust for peroxides: if O appears to be -2 but total doesn't match, treat O as -1
        if (compMap.containsKey("O")) {
            int Ocount = compMap.get("O");
            if (Ocount == 2 && oxMap.getOrDefault("O", -2) == -2) {
                // Check if using -2 for O leads to incorrect total charge
                int totalCalc = 0;
                for (String elem : compMap.keySet()) {
                    totalCalc += oxMap.getOrDefault(elem, 0) * compMap.get(elem);
                }
                if (totalCalc != netCharge) {
                    // Assume it's a peroxide: set O to -1
                    oxMap.put("O", -1);
                    // Recompute one unknown element's oxidation state if applicable
                    for (String elem : unknowns) {
                        if (!elem.equals("O")) {
                            int count = compMap.get(elem);
                            int ox = (netCharge - (sumKnown + Ocount)) / count;
                            oxMap.put(elem, ox);
                        }
                    }
                }
            }
        }
        return oxMap;
    }

    // --- Helper Methods for Redox Balancing ---
    private static String formatHalfReaction(HalfReaction half) {
        // Construct a string like "Reactant + extras -> Product + extras (+ e- if any)" for half-reaction
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
    private static void addSpecies(Map<String,Integer> map, String species, int count) {
        if (species == null || species.isEmpty() || count == 0) return;
        map.put(species, map.getOrDefault(species, 0) + count);
    }
    private static void cancelSpecies(Map<String,Integer> left, Map<String,Integer> right, String species) {
        int leftCount = left.getOrDefault(species, 0);
        int rightCount = right.getOrDefault(species, 0);
        int cancel = Math.min(leftCount, rightCount);
        if (cancel > 0) {
            left.put(species, leftCount - cancel);
            right.put(species, rightCount - cancel);
        }
    }
    private static String formatEquation(Map<String,Integer> left, Map<String,Integer> right) {
        // Format maps of species into "coef Species + ... -> coef Species + ..."
        StringBuilder sb = new StringBuilder();
        boolean firstTerm = true;
        for (Map.Entry<String,Integer> entry : left.entrySet()) {
            int count = entry.getValue();
            String species = entry.getKey();
            if (count == 0) continue;
            if (!firstTerm) sb.append(" + ");
            if (count > 1) sb.append(count).append(" ");
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
            if (count > 1) sb.append(count).append(" ");
            sb.append(species);
            firstTerm = false;
        }
        return sb.toString();
    }
    private static int extraElementCount(List<Compound> extras, String element) {
        // Count occurrences of an element in a list of extra compounds
        int total = 0;
        for (Compound comp : extras) {
            total += comp.composition.getOrDefault(element, 0);
        }
        return total;
    }
    private static int extraCharge(List<Compound> extras) {
        // Sum charges of a list of extra compounds
        int total = 0;
        for (Compound comp : extras) {
            total += comp.charge;
        }
        return total;
    }
    private static int totalCharge(Compound comp) {
        return comp.charge;
    }
    private static int lcm(int a, int b) {
        if (a == 0 || b == 0) return 0;
        int gcd = gcd(a, b);
        return Math.abs(a / gcd * b);
    }
    private static int gcd(int a, int b) {
        return BigInteger.valueOf(a).gcd(BigInteger.valueOf(b)).intValue();
    }

    // --- Utility Helpers ---
    private static boolean isMetal(String element) {
        // Simplified: treat everything not explicitly a known nonmetal as metal
        // Known nonmetals and metalloids set:
        Set<String> nonmetals = new HashSet<>(Arrays.asList(
            "H","He","B","C","N","O","F","Ne","Si","P","S","Cl","Ar",
            "As","Se","Br","Kr","Te","I","Xe","At","Rn","Og"
        ));
        return !nonmetals.contains(element);
    }
    private static boolean isGroup1(String element) {
        return Arrays.asList("Li","Na","K","Rb","Cs","Fr").contains(element);
    }
    private static boolean isGroup2(String element) {
        return Arrays.asList("Be","Mg","Ca","Sr","Ba","Ra").contains(element);
    }
    private static boolean isElementIn(String symbol, String... list) {
        for (String s : list) {
            if (s.equals(symbol)) return true;
        }
        return false;
    }

    // Array of element symbols in atomic number order (populated on first use)
    private static String[] ELEMENT_SYMBOLS = null;

    // Load atomic weights for elements 1-118 into a LinkedHashMap to preserve order
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
        weights.put("Ar", 39.948);
        weights.put("K", 39.098);
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
}
