package com.chemcalc.ui;

import com.chemcalc.data.BundledPeriodicTableRepository;
import com.chemcalc.domain.BalanceConservationCheck;
import com.chemcalc.domain.BalanceMatrixSnapshot;
import com.chemcalc.domain.BalanceNormalizationTrace;
import com.chemcalc.domain.BalanceRowOperation;
import com.chemcalc.domain.BalancedReaction;
import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.ExplanationStep;
import com.chemcalc.domain.MolarMassComponent;
import com.chemcalc.domain.MolarMassResult;
import com.chemcalc.domain.NomenclatureResult;
import com.chemcalc.domain.PercentCompositionEntry;
import com.chemcalc.domain.Species;
import com.chemcalc.domain.StructureCandidate;
import com.chemcalc.domain.StructureInputFormat;
import com.chemcalc.domain.StructureRequest;
import com.chemcalc.domain.StructureResolutionStatus;
import com.chemcalc.domain.StructureResult;
import com.chemcalc.domain.ValidationMessage;
import com.chemcalc.engine.BalancerService;
import com.chemcalc.engine.DefaultFormulaParser;
import com.chemcalc.engine.DefaultReactionParser;
import com.chemcalc.engine.DefaultStoichiometryService;
import com.chemcalc.engine.FormulaParser;
import com.chemcalc.engine.LinearEquationBalancer;
import com.chemcalc.engine.NomenclatureService;
import com.chemcalc.engine.ReactionParser;
import com.chemcalc.engine.RuleBasedNomenclatureService;
import com.chemcalc.engine.StoichiometryService;
import com.chemcalc.engine.StructureService;
import com.chemcalc.structurecdk.CdkStructureService;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.function.Consumer;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.beans.property.SimpleObjectProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Worker;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.Separator;
import javafx.scene.control.SplitPane;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.TitledPane;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.web.WebEngine;
import javafx.scene.web.WebView;
import javafx.stage.Stage;
import netscape.javascript.JSObject;

public final class ChemistrySuperApp extends Application {
    private final BundledPeriodicTableRepository repository = new BundledPeriodicTableRepository();
    private final FormulaParser formulaParser = new DefaultFormulaParser();
    private final ReactionParser reactionParser = new DefaultReactionParser(formulaParser);
    private final BalancerService balancerService = new LinearEquationBalancer();
    private final StoichiometryService stoichiometryService = new DefaultStoichiometryService(repository);
    private final NomenclatureService nomenclatureService = new RuleBasedNomenclatureService(repository);
    private final StructureService structureService = new CdkStructureService();
    private ChromiumStructureViewer chromiumStructureViewer;

    @Override
    public void start(Stage stage) {
        TabPane tabs = new TabPane(
                createWorkbenchTab(),
                createPeriodicTableTab(),
                createVisualizationTab());

        Scene scene = new Scene(tabs, 1520, 960);
        stage.setTitle("Chemistry Super App");
        stage.setScene(scene);
        stage.show();
    }

    @Override
    public void stop() {
        if (chromiumStructureViewer != null) {
            chromiumStructureViewer.dispose();
        }
    }

    private Tab createWorkbenchTab() {
        TextField equationField = new TextField("KMnO4 + HCl -> KCl + MnCl2 + Cl2 + H2O");
        equationField.setPromptText("Enter an equation to balance");
        HBox.setHgrow(equationField, Priority.ALWAYS);

        Button balanceButton = new Button("Balance");
        Label balancedEquationLabel = new Label("Balanced equation will appear here.");
        balancedEquationLabel.setWrapText(true);
        balancedEquationLabel.setStyle("-fx-font-size: 20px; -fx-font-weight: bold;");

        Label statusLabel = new Label("Ready to solve.");
        statusLabel.setWrapText(true);
        statusLabel.setStyle("-fx-text-fill: #4b5563;");

        TextArea explanationArea = createOutputArea();
        explanationArea.setPrefRowCount(7);
        explanationArea.setStyle("-fx-font-family: Consolas, 'Courier New', monospace;");

        TableView<BalanceConservationCheck> conservationTable = createConservationTable();
        conservationTable.setPlaceholder(new Label("Balance a reaction to view conservation totals."));
        conservationTable.setPrefHeight(220);

        TextArea mathTraceArea = createOutputArea();
        mathTraceArea.setStyle("-fx-font-family: Consolas, 'Courier New', monospace;");
        TitledPane tracePane = new TitledPane("Exact Matrix Trace", mathTraceArea);
        tracePane.setExpanded(false);

        balanceButton.setOnAction(event -> {
            try {
                BalancedReaction balanced = balancerService.balance(reactionParser.parse(equationField.getText()));
                balancedEquationLabel.setText(balanced.format());
                statusLabel.setText(renderBalanceStatus(balanced));
                explanationArea.setText(renderChemistryExplanation(balanced));
                conservationTable.setItems(FXCollections.observableArrayList(balanced.trace().conservationChecks()));
                mathTraceArea.setText(renderMathTrace(balanced));
            } catch (Exception exception) {
                balancedEquationLabel.setText("Unable to balance equation");
                statusLabel.setText(exception.getMessage());
                explanationArea.clear();
                conservationTable.getItems().clear();
                mathTraceArea.clear();
            }
        });

        VBox balancePanel = new VBox(
                12,
                sectionLabel("Reaction Workbench"),
                inputRow(equationField, balanceButton),
                card("Balanced Equation", balancedEquationLabel, statusLabel),
                card("Chemistry Explanation", explanationArea),
                card("Conservation Totals", conservationTable),
                tracePane);
        balancePanel.setPadding(new Insets(16));
        VBox.setVgrow(conservationTable, Priority.NEVER);
        VBox.setVgrow(tracePane, Priority.ALWAYS);

        TextField formulaField = new TextField("CuSO4\u00b75H2O");
        formulaField.setPromptText("Enter a formula");
        HBox.setHgrow(formulaField, Priority.ALWAYS);

        Button analyzeButton = new Button("Analyze");
        TextArea formulaOutput = createOutputArea();

        analyzeButton.setOnAction(event -> {
            try {
                formulaOutput.clear();
                Species species = formulaParser.parseSpecies(formulaField.getText());
                MolarMassResult molarMassResult = stoichiometryService.molarMass(species);
                NomenclatureResult nomenclatureResult = nomenclatureService.name(species);

                formulaOutput.appendText("Canonical Formula\n");
                formulaOutput.appendText(species.displayFormula() + "\n\n");
                formulaOutput.appendText("Composition\n");
                species.composition().forEach((symbol, count) -> formulaOutput.appendText("- " + symbol + ": " + count + "\n"));
                formulaOutput.appendText("\nMolar Mass\n");
                formulaOutput.appendText(String.format("%.4f g/mol%n", molarMassResult.totalMass()));
                for (MolarMassComponent component : molarMassResult.components()) {
                    formulaOutput.appendText(String.format("- %s x%d -> %.4f g/mol%n",
                            component.symbol(), component.count(), component.contribution()));
                }
                formulaOutput.appendText("\nPercent Composition\n");
                for (PercentCompositionEntry entry : stoichiometryService.percentComposition(species).entries()) {
                    formulaOutput.appendText(String.format("- %s: %.2f%%%n", entry.symbol(), entry.percent()));
                }
                formulaOutput.appendText("\nName\n");
                formulaOutput.appendText(nomenclatureResult.supported() ? nomenclatureResult.displayName() : "Unsupported");
                if (!nomenclatureResult.messages().isEmpty()) {
                    formulaOutput.appendText("\n\nWarnings\n");
                    nomenclatureResult.messages().forEach(message -> formulaOutput.appendText("- " + message.message() + "\n"));
                }
            } catch (Exception exception) {
                formulaOutput.setText("Unable to analyze formula:\n" + exception.getMessage());
            }
        });

        VBox analysisPanel = new VBox(
                12,
                sectionLabel("Formula Analysis"),
                inputRow(formulaField, analyzeButton),
                card("Analysis Output", formulaOutput));
        analysisPanel.setPadding(new Insets(16));

        SplitPane splitPane = new SplitPane(balancePanel, analysisPanel);
        splitPane.setDividerPositions(0.58);

        Tab tab = new Tab("Workbench", splitPane);
        tab.setClosable(false);
        return tab;
    }

    private Tab createPeriodicTableTab() {
        List<ElementRecord> allElements = repository.findAll();

        TextField searchField = new TextField();
        searchField.setPromptText("Search by symbol, name, or atomic number");
        HBox.setHgrow(searchField, Priority.ALWAYS);

        ComboBox<String> categoryFilter = createFilterCombo("All categories",
                distinctValues(allElements, ElementRecord::category));
        ComboBox<String> blockFilter = createFilterCombo("All blocks",
                distinctValues(allElements, ElementRecord::block));
        ComboBox<String> periodFilter = createFilterCombo("All periods",
                allElements.stream()
                        .map(element -> Integer.toString(element.period()))
                        .distinct()
                        .sorted(Comparator.comparingInt(Integer::parseInt))
                        .toList());
        ComboBox<String> stateFilter = createFilterCombo("All states",
                allElements.stream()
                        .map(this::standardStateLabel)
                        .distinct()
                        .sorted()
                        .toList());

        Label resultCountLabel = new Label("118 elements visible");
        resultCountLabel.setStyle("-fx-text-fill: #4b5563;");

        WebView tableView = new WebView();
        tableView.setPrefWidth(980);
        WebEngine engine = tableView.getEngine();
        engine.loadContent(PeriodicTableHtmlRenderer.renderShell());

        PeriodicTableDetailPane detailPane = new PeriodicTableDetailPane();
        SimpleObjectProperty<ElementRecord> selectedElement = new SimpleObjectProperty<>();
        ObservableList<ElementRecord> filteredElements = FXCollections.observableArrayList(allElements);

        PeriodicTableBridge bridge = new PeriodicTableBridge(symbol -> {
            ElementRecord match = filteredElements.stream()
                    .filter(element -> element.symbol().equalsIgnoreCase(symbol))
                    .findFirst()
                    .orElse(null);
            Platform.runLater(() -> selectedElement.set(match));
        });

        final boolean[] pageReady = {false};
        engine.getLoadWorker().stateProperty().addListener((observable, oldState, newState) -> {
            if (newState != Worker.State.SUCCEEDED) {
                return;
            }
            JSObject window = (JSObject) engine.executeScript("window");
            window.setMember("periodicBridge", bridge);
            pageReady[0] = true;
            renderPeriodicTable(engine, filteredElements, selectedElement.get());
        });

        Runnable refresh = () -> {
            List<ElementRecord> matches = filterElements(
                    allElements,
                    searchField.getText(),
                    categoryFilter.getValue(),
                    blockFilter.getValue(),
                    periodFilter.getValue(),
                    stateFilter.getValue());
            filteredElements.setAll(matches);

            ElementRecord currentSelection = selectedElement.get();
            if (currentSelection != null && matches.stream().noneMatch(element -> element.symbol().equals(currentSelection.symbol()))) {
                selectedElement.set(null);
            }

            resultCountLabel.setText(matches.size() + (matches.size() == 1 ? " element visible" : " elements visible"));
            detailPane.show(selectedElement.get(), matches.size());
            if (pageReady[0]) {
                renderPeriodicTable(engine, filteredElements, selectedElement.get());
            }
        };

        selectedElement.addListener((observable, oldValue, newValue) -> {
            detailPane.show(newValue, filteredElements.size());
            if (pageReady[0]) {
                renderPeriodicTable(engine, filteredElements, newValue);
            }
        });

        searchField.textProperty().addListener((observable, oldValue, newValue) -> refresh.run());
        categoryFilter.valueProperty().addListener((observable, oldValue, newValue) -> refresh.run());
        blockFilter.valueProperty().addListener((observable, oldValue, newValue) -> refresh.run());
        periodFilter.valueProperty().addListener((observable, oldValue, newValue) -> refresh.run());
        stateFilter.valueProperty().addListener((observable, oldValue, newValue) -> refresh.run());

        HBox filters = new HBox(10, searchField, categoryFilter, blockFilter, periodFilter, stateFilter);
        filters.setPadding(new Insets(12, 12, 0, 12));

        VBox leftPane = new VBox(10, filters, resultCountLabel, tableView);
        leftPane.setPadding(new Insets(0, 0, 12, 12));
        VBox.setVgrow(tableView, Priority.ALWAYS);

        SplitPane splitPane = new SplitPane(leftPane, detailPane.root());
        splitPane.setDividerPositions(0.74);

        BorderPane root = new BorderPane(splitPane);
        detailPane.show(null, filteredElements.size());
        refresh.run();

        Tab tab = new Tab("Periodic Table", root);
        tab.setClosable(false);
        return tab;
    }

    private Tab createVisualizationTab() {
        ComboBox<StructureInputFormat> inputFormat = new ComboBox<>(FXCollections.observableArrayList(StructureInputFormat.values()));
        inputFormat.setValue(StructureInputFormat.FORMULA);
        inputFormat.setPrefWidth(140);

        ComboBox<String> stylePreset = new ComboBox<>(FXCollections.observableArrayList("ball-stick", "stick", "sphere", "line"));
        stylePreset.setValue("ball-stick");
        stylePreset.setPrefWidth(130);

        TextField singleLineInput = new TextField(defaultInputFor(StructureInputFormat.FORMULA));
        HBox.setHgrow(singleLineInput, Priority.ALWAYS);

        TextArea multiLineInput = new TextArea();
        multiLineInput.setPromptText("Paste a Molfile or SDF block here");
        multiLineInput.setWrapText(false);
        multiLineInput.setVisible(false);
        multiLineInput.setManaged(false);
        multiLineInput.setPrefRowCount(8);

        Button loadButton = new Button("Load Structure");
        Button resetButton = new Button("Reset");

        ComboBox<StructureCandidate> candidateBox = new ComboBox<>();
        candidateBox.setPromptText("Select a candidate");
        candidateBox.setVisible(false);
        candidateBox.setManaged(false);
        candidateBox.setMaxWidth(Double.MAX_VALUE);

        Label statusLabel = new Label("Ready to resolve a structure.");
        statusLabel.setWrapText(true);
        statusLabel.setStyle("-fx-text-fill: #4b5563;");

        TextArea messageArea = createOutputArea();
        messageArea.setPrefRowCount(6);
        messageArea.setStyle("-fx-font-family: Consolas, 'Courier New', monospace;");

        StructureMetadataPane metadataPane = new StructureMetadataPane();
        SimpleObjectProperty<StructureResult> currentResult = new SimpleObjectProperty<>();
        final boolean[] candidateUpdate = {false};

        chromiumStructureViewer = new ChromiumStructureViewer(
                (atomIndex, label) -> metadataPane.updateSelectedAtom(atomIndex, label));

        candidateBox.setCellFactory(listView -> new javafx.scene.control.ListCell<>() {
            @Override
            protected void updateItem(StructureCandidate item, boolean empty) {
                super.updateItem(item, empty);
                setText(empty || item == null ? null : item.preferredName() + " [" + item.canonicalFormula() + "]");
            }
        });
        candidateBox.setButtonCell(new javafx.scene.control.ListCell<>() {
            @Override
            protected void updateItem(StructureCandidate item, boolean empty) {
                super.updateItem(item, empty);
                setText(empty || item == null ? null : item.preferredName() + " [" + item.canonicalFormula() + "]");
            }
        });

        Runnable loadStructure = () -> {
            StructureInputFormat format = inputFormat.getValue();
            String input = format == StructureInputFormat.MOLFILE ? multiLineInput.getText() : singleLineInput.getText();
            String selectedCandidateKey = candidateBox.isVisible() && candidateBox.getValue() != null
                    ? candidateBox.getValue().id()
                    : "";

            StructureResult result = structureService.resolve(new StructureRequest(format, input, selectedCandidateKey));
            currentResult.set(result);

            candidateUpdate[0] = true;
            candidateBox.setItems(FXCollections.observableArrayList(result.candidates()));
            boolean showCandidates = !result.candidates().isEmpty();
            candidateBox.setVisible(showCandidates);
            candidateBox.setManaged(showCandidates);
            if (showCandidates && !selectedCandidateKey.isBlank()) {
                result.candidates().stream()
                        .filter(candidate -> candidate.id().equalsIgnoreCase(selectedCandidateKey))
                        .findFirst()
                        .ifPresent(candidateBox::setValue);
            } else if (!showCandidates) {
                candidateBox.setValue(null);
            }
            candidateUpdate[0] = false;

            statusLabel.setText(renderStructureStatus(result));
            messageArea.setText(renderStructureMessages(result));
            metadataPane.show(result);
            renderStructureView(chromiumStructureViewer, result, stylePreset.getValue());
        };

        inputFormat.valueProperty().addListener((observable, oldValue, newValue) -> {
            boolean molfileMode = newValue == StructureInputFormat.MOLFILE;
            singleLineInput.setVisible(!molfileMode);
            singleLineInput.setManaged(!molfileMode);
            multiLineInput.setVisible(molfileMode);
            multiLineInput.setManaged(molfileMode);
            if (molfileMode) {
                if (multiLineInput.getText().isBlank()) {
                    multiLineInput.setText(defaultInputFor(newValue));
                }
            } else {
                if (singleLineInput.getText().isBlank()) {
                    singleLineInput.setText(defaultInputFor(newValue));
                }
            }
            candidateBox.setItems(FXCollections.observableArrayList());
            candidateBox.setVisible(false);
            candidateBox.setManaged(false);
            statusLabel.setText("Input mode changed to " + newValue + ".");
            metadataPane.resetSelectedAtom();
        });

        stylePreset.valueProperty().addListener((observable, oldValue, newValue) -> {
            chromiumStructureViewer.setStylePreset(newValue);
        });

        candidateBox.valueProperty().addListener((observable, oldValue, newValue) -> {
            if (candidateUpdate[0] || newValue == null) {
                return;
            }
            loadStructure.run();
        });

        loadButton.setOnAction(event -> loadStructure.run());
        resetButton.setOnAction(event -> {
            inputFormat.setValue(StructureInputFormat.FORMULA);
            singleLineInput.setText(defaultInputFor(StructureInputFormat.FORMULA));
            multiLineInput.clear();
            candidateBox.setItems(FXCollections.observableArrayList());
            candidateBox.setValue(null);
            candidateBox.setVisible(false);
            candidateBox.setManaged(false);
            statusLabel.setText("Reset to the default formula example.");
            metadataPane.show(null);
            currentResult.set(null);
            renderStructureView(chromiumStructureViewer, null, stylePreset.getValue());
        });

        VBox controls = new VBox(
                10,
                sectionLabel("Visualization"),
                new HBox(10, inputFormat, stylePreset),
                new HBox(8, singleLineInput, loadButton, resetButton),
                multiLineInput,
                candidateBox,
                card("Status", statusLabel),
                card("Messages", messageArea),
                metadataPane.root());
        controls.setPadding(new Insets(16));
        controls.setPrefWidth(440);
        VBox.setVgrow(messageArea, Priority.NEVER);
        VBox.setVgrow(metadataPane.root(), Priority.ALWAYS);

        BorderPane viewerPane = new BorderPane(chromiumStructureViewer.node());
        viewerPane.setPadding(new Insets(16, 16, 16, 0));

        SplitPane splitPane = new SplitPane(viewerPane, controls);
        splitPane.setDividerPositions(0.68);

        loadStructure.run();

        Tab tab = new Tab("Visualization", splitPane);
        tab.setClosable(false);
        return tab;
    }

    private TableView<BalanceConservationCheck> createConservationTable() {
        TableView<BalanceConservationCheck> table = new TableView<>();
        table.setColumnResizePolicy(features -> TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN.call(features));

        TableColumn<BalanceConservationCheck, String> quantityColumn = new TableColumn<>("Quantity");
        quantityColumn.setCellValueFactory(cell -> new ReadOnlyStringWrapper(cell.getValue().label()));

        TableColumn<BalanceConservationCheck, Number> reactantColumn = new TableColumn<>("Reactants");
        reactantColumn.setCellValueFactory(cell -> new ReadOnlyObjectWrapper<>(cell.getValue().reactantTotal()));

        TableColumn<BalanceConservationCheck, Number> productColumn = new TableColumn<>("Products");
        productColumn.setCellValueFactory(cell -> new ReadOnlyObjectWrapper<>(cell.getValue().productTotal()));

        TableColumn<BalanceConservationCheck, String> statusColumn = new TableColumn<>("Status");
        statusColumn.setCellValueFactory(cell -> new ReadOnlyStringWrapper(
                cell.getValue().reactantTotal() == cell.getValue().productTotal() ? "Balanced" : "Mismatch"));

        table.getColumns().setAll(quantityColumn, reactantColumn, productColumn, statusColumn);
        return table;
    }

    private ComboBox<String> createFilterCombo(String allLabel, List<String> values) {
        ComboBox<String> comboBox = new ComboBox<>();
        comboBox.setItems(FXCollections.observableArrayList(allLabel));
        comboBox.getItems().addAll(values);
        comboBox.setValue(allLabel);
        comboBox.setPrefWidth(145);
        return comboBox;
    }

    private List<String> distinctValues(List<ElementRecord> elements, java.util.function.Function<ElementRecord, String> extractor) {
        return elements.stream()
                .map(extractor)
                .filter(value -> value != null && !value.isBlank())
                .distinct()
                .sorted()
                .toList();
    }

    private List<ElementRecord> filterElements(
            List<ElementRecord> elements,
            String query,
            String category,
            String block,
            String period,
            String standardState) {
        String normalizedQuery = query == null ? "" : query.trim().toLowerCase(Locale.ROOT);
        return elements.stream()
                .filter(element -> normalizedQuery.isBlank()
                        || element.symbol().toLowerCase(Locale.ROOT).contains(normalizedQuery)
                        || element.name().toLowerCase(Locale.ROOT).contains(normalizedQuery)
                        || Integer.toString(element.atomicNumber()).equals(normalizedQuery))
                .filter(element -> isAllFilter(category) || element.category().equals(category))
                .filter(element -> isAllFilter(block) || element.block().equals(block))
                .filter(element -> isAllFilter(period) || Integer.toString(element.period()).equals(period))
                .filter(element -> isAllFilter(standardState) || standardStateLabel(element).equals(standardState))
                .toList();
    }

    private boolean isAllFilter(String value) {
        return value == null || value.startsWith("All ");
    }

    private String standardStateLabel(ElementRecord element) {
        return element.standardState().isBlank() ? "Unknown" : element.standardState();
    }

    private void renderPeriodicTable(WebEngine engine, List<ElementRecord> elements, ElementRecord selectedElement) {
        engine.executeScript(PeriodicTableHtmlRenderer.renderUpdateScript(
                elements,
                selectedElement == null ? null : selectedElement.symbol()));
    }

    private Label sectionLabel(String text) {
        Label label = new Label(text);
        label.setStyle("-fx-font-size: 22px; -fx-font-weight: bold;");
        return label;
    }

    private HBox inputRow(Region input, Button button) {
        HBox row = new HBox(8, input, button);
        row.setFillHeight(true);
        return row;
    }

    private VBox card(String title, Region content) {
        Label header = new Label(title);
        header.setStyle("-fx-font-size: 15px; -fx-font-weight: bold;");
        VBox card = new VBox(8, header, content);
        card.setPadding(new Insets(12));
        card.setStyle("-fx-background-color: #ffffff; -fx-background-radius: 16; -fx-border-color: #dbe3ec; -fx-border-radius: 16;");
        return card;
    }

    private VBox card(String title, Region primary, Region secondary) {
        Label header = new Label(title);
        header.setStyle("-fx-font-size: 15px; -fx-font-weight: bold;");
        VBox card = new VBox(8, header, primary, secondary);
        card.setPadding(new Insets(12));
        card.setStyle("-fx-background-color: #ffffff; -fx-background-radius: 16; -fx-border-color: #dbe3ec; -fx-border-radius: 16;");
        return card;
    }

    private TextArea createOutputArea() {
        TextArea area = new TextArea();
        area.setWrapText(true);
        area.setEditable(false);
        VBox.setVgrow(area, Priority.ALWAYS);
        return area;
    }

    private void renderStructureView(ChromiumStructureViewer viewer, StructureResult result, String stylePreset) {
        if (result == null) {
            viewer.renderScript(StructureViewerHtmlRenderer.renderMessageScript(
                    "No structure loaded",
                    "Load a formula, SMILES, Molfile/SDF, or InChI to populate the viewer."), stylePreset);
            return;
        }
        if (result.hasRenderablePayload()) {
            viewer.renderScript(StructureViewerHtmlRenderer.renderStructureScript(result), stylePreset);
            return;
        }
        String detail = result.messages().isEmpty()
                ? "No renderable structure payload is available."
                : result.messages().get(0).message();
        viewer.renderScript(StructureViewerHtmlRenderer.renderMessageScript(
                result.status().name(),
                detail), stylePreset);
    }

    private String renderStructureStatus(StructureResult result) {
        return switch (result.status()) {
            case RESOLVED -> "Resolved " + result.preferredName() + " with "
                    + result.atomCount() + " atoms and " + result.bondCount() + " bonds.";
            case AMBIGUOUS -> "Multiple known structures match " + result.canonicalFormula() + ".";
            case UNSUPPORTED -> "No offline structure matches " + result.canonicalFormula()
                    + "; enter SMILES, Molfile/SDF, or InChI for an explicit structure.";
            case INVALID -> "The structure input could not be parsed.";
        };
    }

    private String renderStructureMessages(StructureResult result) {
        if (result.messages().isEmpty()) {
            return "No warnings.";
        }
        StringBuilder builder = new StringBuilder();
        for (ValidationMessage message : result.messages()) {
            builder.append(message.severity()).append(" ")
                    .append(message.code()).append(": ")
                    .append(message.message()).append("\n");
        }
        return builder.toString().trim();
    }

    private String defaultInputFor(StructureInputFormat format) {
        return switch (format) {
            case FORMULA -> "H2O";
            case SMILES -> "CCO";
            case MOLFILE -> "";
            case INCHI -> "InChI=1S/H2O/h1H2";
        };
    }

    private String renderBalanceStatus(BalancedReaction balanced) {
        int balancedChecks = (int) balanced.trace().conservationChecks().stream()
                .filter(check -> check.reactantTotal() == check.productTotal())
                .count();
        return "Balanced " + balancedChecks + " conserved quantities across "
                + balanced.trace().rowOperations().size() + " matrix operations.";
    }

    private String renderChemistryExplanation(BalancedReaction balanced) {
        StringBuilder builder = new StringBuilder();
        builder.append("Balanced equation:\n");
        builder.append(balanced.format()).append("\n\n");

        builder.append("Why these coefficients work:\n");
        for (ExplanationStep step : balanced.explanation()) {
            builder.append("- ").append(step.title()).append(": ").append(step.detail()).append("\n");
        }

        builder.append("\nConserved quantities:\n");
        balanced.trace().conservationChecks().forEach(check -> builder.append("- ")
                .append(check.label())
                .append(": ")
                .append(check.reactantTotal())
                .append(" on reactants, ")
                .append(check.productTotal())
                .append(" on products\n"));
        return builder.toString();
    }

    private String renderMathTrace(BalancedReaction balanced) {
        StringBuilder builder = new StringBuilder();
        builder.append("Initial Conservation Matrix\n");
        builder.append(formatMatrix(balanced.trace().initialMatrix())).append("\n");

        int stepNumber = 1;
        for (BalanceRowOperation operation : balanced.trace().rowOperations()) {
            builder.append("Step ").append(stepNumber++).append(": ").append(operation.label()).append("\n");
            builder.append(operation.detail()).append("\n");
            builder.append(formatMatrix(operation.snapshot())).append("\n");
        }

        builder.append("Nullspace Solution\n");
        builder.append("Pivot variables: ").append(String.join(", ", balanced.trace().solution().pivotVariables())).append("\n");
        builder.append("Free variables: ").append(String.join(", ", balanced.trace().solution().freeVariables())).append("\n");
        builder.append("Vector: ").append(String.join(", ", balanced.trace().solution().nullspaceVector())).append("\n");
        builder.append(balanced.trace().solution().detail()).append("\n\n");

        BalanceNormalizationTrace normalization = balanced.trace().normalization();
        builder.append("Normalization\n");
        builder.append("LCM of denominators: ").append(normalization.leastCommonMultiple()).append("\n");
        builder.append("Scaled vector: ").append(String.join(", ", normalization.scaledVector())).append("\n");
        builder.append("Sign flip applied: ").append(normalization.signFlipped() ? "yes" : "no").append("\n");
        builder.append("Positive vector: ").append(String.join(", ", normalization.signNormalizedVector())).append("\n");
        builder.append("GCD reduction: ").append(normalization.greatestCommonDivisor()).append("\n");
        builder.append("Final coefficients: ").append(normalization.coefficients()).append("\n");
        return builder.toString();
    }

    private String formatMatrix(BalanceMatrixSnapshot snapshot) {
        if (snapshot.rowLabels().isEmpty() || snapshot.columnLabels().isEmpty()) {
            return "(empty matrix)\n";
        }

        int[] widths = new int[snapshot.columnLabels().size() + 1];
        widths[0] = snapshot.rowLabels().stream().mapToInt(String::length).max().orElse(6);
        for (int column = 0; column < snapshot.columnLabels().size(); column++) {
            int valueWidth = 1;
            for (List<String> row : snapshot.values()) {
                valueWidth = Math.max(valueWidth, row.get(column).length());
            }
            widths[column + 1] = Math.max(snapshot.columnLabels().get(column).length(), valueWidth);
        }

        StringBuilder builder = new StringBuilder();
        builder.append(padRight("", widths[0])).append(" | ");
        for (int column = 0; column < snapshot.columnLabels().size(); column++) {
            if (column > 0) {
                builder.append(" ");
            }
            builder.append(padLeft(snapshot.columnLabels().get(column), widths[column + 1]));
        }
        builder.append("\n");

        for (int row = 0; row < snapshot.rowLabels().size(); row++) {
            builder.append(padRight(snapshot.rowLabels().get(row), widths[0])).append(" | ");
            for (int column = 0; column < snapshot.columnLabels().size(); column++) {
                if (column > 0) {
                    builder.append(" ");
                }
                builder.append(padLeft(snapshot.values().get(row).get(column), widths[column + 1]));
            }
            builder.append("\n");
        }
        return builder.toString();
    }

    private String padLeft(String value, int width) {
        return " ".repeat(Math.max(0, width - value.length())) + value;
    }

    private String padRight(String value, int width) {
        return value + " ".repeat(Math.max(0, width - value.length()));
    }

    public static void main(String[] args) {
        launch(args);
    }

    private static final class PeriodicTableDetailPane {
        private final Label title = new Label("Periodic Table Explorer");
        private final Label subtitle = new Label("Select an element to inspect its metadata.");
        private final Label symbol = new Label("--");
        private final Label name = new Label("No element selected");
        private final GridPane facts = new GridPane();
        private final Label source = new Label("Source will appear here.");
        private final VBox root;
        private final Label visibleCount = new Label();

        private PeriodicTableDetailPane() {
            title.setStyle("-fx-font-size: 24px; -fx-font-weight: bold;");
            subtitle.setWrapText(true);
            subtitle.setStyle("-fx-text-fill: #4b5563;");

            symbol.setStyle("-fx-font-size: 48px; -fx-font-weight: bold;");
            name.setStyle("-fx-font-size: 18px;");

            facts.setHgap(12);
            facts.setVgap(8);

            source.setWrapText(true);
            source.setStyle("-fx-text-fill: #4b5563;");

            VBox hero = new VBox(2, visibleCount, symbol, name);
            hero.setPadding(new Insets(12));
            hero.setStyle("-fx-background-color: linear-gradient(to bottom right, #f7fafc, #e2e8f0); -fx-background-radius: 18;");

            root = new VBox(
                    14,
                    title,
                    subtitle,
                    hero,
                    new Separator(Orientation.HORIZONTAL),
                    facts,
                    new Separator(Orientation.HORIZONTAL),
                    source);
            root.setPadding(new Insets(18));
            root.setStyle("-fx-background-color: #ffffff; -fx-border-color: #dbe3ec; -fx-border-radius: 0; -fx-background-radius: 0;");
            VBox.setVgrow(facts, Priority.ALWAYS);
        }

        private VBox root() {
            return root;
        }

        private void show(ElementRecord element, int matchCount) {
            visibleCount.setText(matchCount + (matchCount == 1 ? " matching element" : " matching elements"));
            facts.getChildren().clear();

            if (element == null) {
                symbol.setText("--");
                name.setText("No element selected");
                subtitle.setText("Use the filters or click a tile in the WebView explorer.");
                source.setText("Select a tile to view atomic number, block, category, state, and provenance.");
                addFactRow("Category", "None");
                addFactRow("State", "None");
                addFactRow("Block", "None");
                addFactRow("Period", "None");
                return;
            }

            symbol.setText(element.symbol());
            name.setText(element.name() + " (#" + element.atomicNumber() + ")");
            subtitle.setText("Metadata updates live as the WebView selection changes.");
            addFactRow("Category", element.category());
            addFactRow("Standard state", element.standardState().isBlank() ? "Unknown" : element.standardState());
            addFactRow("Block", element.block());
            addFactRow("Group", element.group() == null ? "Unknown" : Integer.toString(element.group()));
            addFactRow("Period", Integer.toString(element.period()));
            addFactRow("Atomic weight", element.property("standardAtomicWeight")
                    .map(value -> value.value() + " " + value.units())
                    .orElse("Unavailable"));
            addFactRow("Year discovered", element.yearDiscovered().isBlank() ? "Unavailable" : element.yearDiscovered());
            source.setText(element.property("standardAtomicWeight")
                    .map(value -> "Atomic weight source: " + value.source().title())
                    .orElse("Atomic weight source: unavailable"));
        }

        private void addFactRow(String labelText, String valueText) {
            int rowIndex = facts.getChildren().size() / 2;
            Label label = new Label(labelText);
            label.setStyle("-fx-font-weight: bold;");
            Label value = new Label(valueText);
            value.setWrapText(true);
            facts.addRow(rowIndex, label, value);
        }
    }

    private static final class StructureMetadataPane {
        private final Label title = new Label("Structure Metadata");
        private final Label preferredName = new Label("No structure loaded");
        private final Label canonicalFormula = new Label("--");
        private final Label source = new Label("Source unavailable");
        private final Label atomCount = new Label("--");
        private final Label bondCount = new Label("--");
        private final Label charge = new Label("--");
        private final Label selectedAtom = new Label("No atom selected");
        private final VBox root;

        private StructureMetadataPane() {
            title.setStyle("-fx-font-size: 22px; -fx-font-weight: bold;");
            preferredName.setStyle("-fx-font-size: 18px; -fx-font-weight: bold;");
            canonicalFormula.setStyle("-fx-font-size: 16px;");

            GridPane facts = new GridPane();
            facts.setHgap(12);
            facts.setVgap(8);

            facts.addRow(0, factLabel("Source"), source);
            facts.addRow(1, factLabel("Atoms"), atomCount);
            facts.addRow(2, factLabel("Bonds"), bondCount);
            facts.addRow(3, factLabel("Formal charge"), charge);
            facts.addRow(4, factLabel("Selected atom"), selectedAtom);

            root = new VBox(10, title, preferredName, canonicalFormula, facts);
            root.setPadding(new Insets(12));
            root.setStyle("-fx-background-color: #ffffff; -fx-background-radius: 16; -fx-border-color: #dbe3ec; -fx-border-radius: 16;");
        }

        private VBox root() {
            return root;
        }

        private void show(StructureResult result) {
            if (result == null) {
                preferredName.setText("No structure loaded");
                canonicalFormula.setText("--");
                source.setText("Source unavailable");
                atomCount.setText("--");
                bondCount.setText("--");
                charge.setText("--");
                selectedAtom.setText("No atom selected");
                return;
            }
            preferredName.setText(result.preferredName().isBlank() ? result.status().name() : result.preferredName());
            canonicalFormula.setText(result.canonicalFormula().isBlank() ? "--" : result.canonicalFormula());
            source.setText(result.source().title().isBlank()
                    ? "Source unavailable"
                    : result.source().title() + " via " + result.inputFormat().name());
            atomCount.setText(result.atomCount() <= 0 ? "--" : Integer.toString(result.atomCount()));
            bondCount.setText(result.bondCount() <= 0 ? "--" : Integer.toString(result.bondCount()));
            charge.setText(Integer.toString(result.formalCharge()));
            selectedAtom.setText("No atom selected");
        }

        private void updateSelectedAtom(String atomIndex, String label) {
            selectedAtom.setText(label + " (index " + atomIndex + ")");
        }

        private void resetSelectedAtom() {
            selectedAtom.setText("No atom selected");
        }

        private Label factLabel(String text) {
            Label label = new Label(text);
            label.setStyle("-fx-font-weight: bold;");
            return label;
        }
    }

    public static final class PeriodicTableBridge {
        private final Consumer<String> onSelect;

        public PeriodicTableBridge(Consumer<String> onSelect) {
            this.onSelect = onSelect;
        }

        public void selectElement(String symbol) {
            onSelect.accept(symbol);
        }
    }
}
