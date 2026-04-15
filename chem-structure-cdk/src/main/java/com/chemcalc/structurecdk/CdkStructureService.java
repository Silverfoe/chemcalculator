package com.chemcalc.structurecdk;

import com.chemcalc.data.BundledCompoundStructureRepository;
import com.chemcalc.domain.CompoundStructureRecord;
import com.chemcalc.domain.CompoundStructureRepository;
import com.chemcalc.domain.CompoundStructureKind;
import com.chemcalc.domain.SourceReference;
import com.chemcalc.domain.Species;
import com.chemcalc.domain.StructureCandidate;
import com.chemcalc.domain.StructureInputFormat;
import com.chemcalc.domain.StructureRequest;
import com.chemcalc.domain.StructureResolutionStatus;
import com.chemcalc.domain.StructureResult;
import com.chemcalc.domain.StructureViewerPayload;
import com.chemcalc.domain.ValidationMessage;
import com.chemcalc.engine.DefaultFormulaParser;
import com.chemcalc.engine.FormulaParser;
import com.chemcalc.engine.StructureService;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.inchi.InChIToStructure;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public final class CdkStructureService implements StructureService {
    private static final SourceReference DIRECT_INPUT_SOURCE =
            new SourceReference("DIRECT_INPUT", "Direct Input", "", "local");

    private final FormulaParser formulaParser;
    private final CompoundStructureRepository compoundStructureRepository;
    private final SmilesParser smilesParser;

    public CdkStructureService() {
        this(new DefaultFormulaParser(), new BundledCompoundStructureRepository());
    }

    public CdkStructureService(FormulaParser formulaParser, CompoundStructureRepository compoundStructureRepository) {
        this.formulaParser = Objects.requireNonNull(formulaParser, "formulaParser");
        this.compoundStructureRepository = Objects.requireNonNull(compoundStructureRepository, "compoundStructureRepository");
        this.smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    }

    @Override
    public boolean available() {
        return true;
    }

    @Override
    public StructureResult resolve(StructureRequest request) {
        if (request.input().isBlank()) {
            return invalid(request.format(), "Input must not be empty.");
        }
        try {
            return switch (request.format()) {
                case FORMULA -> resolveFormula(request);
                case SMILES -> resolveSmiles(request);
                case MOLFILE -> resolveMolfile(request);
                case INCHI -> resolveInchi(request);
            };
        } catch (Exception exception) {
            return new StructureResult(
                    StructureResolutionStatus.INVALID,
                    request.format(),
                    "",
                    "",
                    SourceReference.unknown(),
                    0,
                    0,
                    0,
                    null,
                    List.of(),
                    List.of(ValidationMessage.error("STRUCTURE_RESOLUTION_FAILED", exception.getMessage())));
        }
    }

    private StructureResult resolveFormula(StructureRequest request) {
        Species species = formulaParser.parseSpecies(request.input());
        String normalizedFormula = toHillFormula(species.composition());
        List<CompoundStructureRecord> matches = compoundStructureRepository.findByFormulaAndCharge(normalizedFormula, species.charge());
        List<StructureCandidate> candidates = matches.stream()
                .map(this::toCandidate)
                .sorted(Comparator.comparing(StructureCandidate::preferredName, String.CASE_INSENSITIVE_ORDER))
                .toList();
        if (matches.isEmpty()) {
            return new StructureResult(
                    StructureResolutionStatus.UNSUPPORTED,
                    StructureInputFormat.FORMULA,
                    "",
                    normalizedFormula,
                    SourceReference.unknown(),
                    0,
                    0,
                    species.charge(),
                    null,
                    List.of(),
                    List.of(ValidationMessage.warning(
                            "FORMULA_NOT_INDEXED",
                            "No offline structure is indexed for " + normalizedFormula + withChargeSuffix(species.charge())
                                    + ". Enter SMILES, Molfile/SDF, or InChI for an explicit structure.")));
        }

        CompoundStructureRecord selected = selectCandidate(matches, request.selectedCandidateKey());
        if (matches.size() > 1 && selected == null) {
            return new StructureResult(
                    StructureResolutionStatus.AMBIGUOUS,
                    StructureInputFormat.FORMULA,
                    "",
                    normalizedFormula,
                    matches.get(0).source(),
                    0,
                    0,
                    species.charge(),
                    null,
                    candidates,
                    List.of(ValidationMessage.warning(
                            "FORMULA_AMBIGUOUS",
                            normalizedFormula + withChargeSuffix(species.charge())
                                    + " matches multiple known structures. Select a candidate to render.")));
        }

        CompoundStructureRecord resolved = selected == null ? matches.get(0) : selected;
        List<ValidationMessage> messages = new ArrayList<>();
        if (matches.size() > 1) {
            messages.add(ValidationMessage.info(
                    "FORMULA_CANDIDATE_SELECTED",
                    "Resolved formula through the selected candidate " + resolved.preferredName() + "."));
        }
        addStructureKindMessage(messages, resolved.structureKind(), resolved.preferredName());
        return resolveCanonicalSmiles(
                resolved.canonicalSmiles(),
                StructureInputFormat.FORMULA,
                resolved.preferredName(),
                resolved.canonicalFormula(),
                resolved.charge(),
                resolved.source(),
                candidates,
                messages);
    }

    private StructureResult resolveSmiles(StructureRequest request) throws CDKException, IOException {
        IAtomContainer container = smilesParser.parseSmiles(request.input());
        return buildResolvedResult(
                container,
                StructureInputFormat.SMILES,
                "Direct SMILES Input",
                "",
                Integer.MIN_VALUE,
                DIRECT_INPUT_SOURCE,
                List.of(),
                List.of());
    }

    private StructureResult resolveMolfile(StructureRequest request) throws CDKException, IOException {
        IAtomContainer container = readMolfile(firstMolBlock(request.input()));
        return buildResolvedResult(
                container,
                StructureInputFormat.MOLFILE,
                "Direct Molfile Input",
                "",
                Integer.MIN_VALUE,
                DIRECT_INPUT_SOURCE,
                List.of(),
                List.of());
    }

    private StructureResult resolveInchi(StructureRequest request) throws Exception {
        InChIToStructure converter = InChIGeneratorFactory.getInstance()
                .getInChIToStructure(request.input(), DefaultChemObjectBuilder.getInstance());
        IAtomContainer container = converter.getAtomContainer();
        return buildResolvedResult(
                container,
                StructureInputFormat.INCHI,
                "Direct InChI Input",
                "",
                Integer.MIN_VALUE,
                DIRECT_INPUT_SOURCE,
                List.of(),
                List.of());
    }

    private StructureResult resolveCanonicalSmiles(
            String canonicalSmiles,
            StructureInputFormat inputFormat,
            String preferredName,
            String canonicalFormula,
            int formalCharge,
            SourceReference source,
            List<StructureCandidate> candidates,
            List<ValidationMessage> baseMessages) {
        try {
            IAtomContainer container = smilesParser.parseSmiles(canonicalSmiles);
            return buildResolvedResult(
                    container,
                    inputFormat,
                    preferredName,
                    canonicalFormula,
                    formalCharge,
                    source,
                    candidates,
                    baseMessages);
        } catch (Exception exception) {
            return new StructureResult(
                    StructureResolutionStatus.INVALID,
                    inputFormat,
                    preferredName,
                    canonicalFormula,
                    source,
                    0,
                    0,
                    formalCharge,
                    null,
                    candidates,
                    List.of(ValidationMessage.error(
                            "CANONICAL_STRUCTURE_FAILED",
                            "Failed to parse the indexed structure for " + preferredName + ": " + exception.getMessage())));
        }
    }

    private StructureResult buildResolvedResult(
            IAtomContainer input,
            StructureInputFormat inputFormat,
            String preferredName,
            String canonicalFormula,
            int formalChargeOverride,
            SourceReference source,
            List<StructureCandidate> candidates,
            List<ValidationMessage> baseMessages) throws CDKException, IOException {
        List<ValidationMessage> messages = new ArrayList<>(baseMessages);
        IAtomContainer container = input;
        normalize(container, messages);

        boolean has3dCoordinates = attempt3dCoordinates(container, messages);
        if (!has3dCoordinates) {
            generate2dCoordinates(container, messages);
        }

        String resolvedFormula = canonicalFormula.isBlank() ? formulaFromContainer(container) : canonicalFormula;
        int resolvedCharge = formalChargeOverride == Integer.MIN_VALUE ? totalFormalCharge(container) : formalChargeOverride;
        String molfile = writeMolfile(container);
        String svg = depictSvg(container, messages);

        return new StructureResult(
                StructureResolutionStatus.RESOLVED,
                inputFormat,
                preferredName,
                resolvedFormula,
                source,
                container.getAtomCount(),
                container.getBondCount(),
                resolvedCharge,
                new StructureViewerPayload("mol", molfile, has3dCoordinates, svg),
                candidates,
                messages);
    }

    private void normalize(IAtomContainer container, List<ValidationMessage> messages) {
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
        } catch (Exception exception) {
            messages.add(ValidationMessage.warning(
                    "ATOM_TYPING_SKIPPED",
                    "Atom typing could not be completed: " + exception.getMessage()));
        }

        try {
            CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(container);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
        } catch (Exception exception) {
            messages.add(ValidationMessage.warning(
                    "HYDROGEN_NORMALIZATION_SKIPPED",
                    "Hydrogen normalization could not be completed: " + exception.getMessage()));
        }

        try {
            Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
            aromaticity.apply(container);
        } catch (Exception exception) {
            messages.add(ValidationMessage.warning(
                    "AROMATICITY_SKIPPED",
                    "Aromaticity perception could not be completed: " + exception.getMessage()));
        }
    }

    private boolean attempt3dCoordinates(IAtomContainer container, List<ValidationMessage> messages) {
        try {
            IAtomContainer with3d = ModelBuilder3D.getInstance(DefaultChemObjectBuilder.getInstance())
                    .generate3DCoordinates(container, false);
            copyAtomCoordinates(with3d, container);
            return GeometryUtil.has3DCoordinates(container);
        } catch (Exception exception) {
            return false;
        }
    }

    private void generate2dCoordinates(IAtomContainer container, List<ValidationMessage> messages) {
        try {
            StructureDiagramGenerator generator = new StructureDiagramGenerator();
            generator.setMolecule(container);
            generator.generateCoordinates();
            IAtomContainer laidOut = generator.getMolecule();
            copyAtomCoordinates(laidOut, container);
        } catch (Exception exception) {
            messages.add(ValidationMessage.warning(
                    "TWO_D_COORDINATES_UNAVAILABLE",
                    "2D coordinate generation failed: " + exception.getMessage()));
        }
    }

    private IAtomContainer readMolfile(String molfile) throws CDKException, IOException {
        if (molfile.contains("V3000")) {
            try (MDLV3000Reader reader = new MDLV3000Reader(new StringReader(molfile))) {
                return reader.read(DefaultChemObjectBuilder.getInstance().newAtomContainer());
            }
        }
        try (MDLV2000Reader reader = new MDLV2000Reader(new StringReader(molfile))) {
            return reader.read(DefaultChemObjectBuilder.getInstance().newAtomContainer());
        }
    }

    private String writeMolfile(IAtomContainer container) throws CDKException, IOException {
        StringWriter writer = new StringWriter();
        try (MDLV2000Writer molWriter = new MDLV2000Writer(writer)) {
            molWriter.write(container);
        }
        return writer.toString();
    }

    private String depictSvg(IAtomContainer container, List<ValidationMessage> messages) {
        try {
            return new DepictionGenerator()
                    .withSize(760, 560)
                    .withFillToFit()
                    .depict(container)
                    .toSvgStr();
        } catch (Exception exception) {
            messages.add(ValidationMessage.warning(
                    "SVG_FALLBACK_UNAVAILABLE",
                    "CDK SVG depiction could not be generated: " + exception.getMessage()));
            return "";
        }
    }

    private void copyAtomCoordinates(IAtomContainer from, IAtomContainer to) {
        for (int index = 0; index < Math.min(from.getAtomCount(), to.getAtomCount()); index++) {
            IAtom sourceAtom = from.getAtom(index);
            IAtom targetAtom = to.getAtom(index);
            targetAtom.setPoint2d(sourceAtom.getPoint2d());
            targetAtom.setPoint3d(sourceAtom.getPoint3d());
        }
    }

    private CompoundStructureRecord selectCandidate(List<CompoundStructureRecord> matches, String selectedCandidateKey) {
        if (selectedCandidateKey == null || selectedCandidateKey.isBlank()) {
            return matches.size() == 1 ? matches.get(0) : null;
        }
        return matches.stream()
                .filter(candidate -> candidate.id().equalsIgnoreCase(selectedCandidateKey))
                .findFirst()
                .orElse(null);
    }

    private StructureCandidate toCandidate(CompoundStructureRecord record) {
        return new StructureCandidate(
                record.id(),
                record.preferredName(),
                record.canonicalFormula(),
                record.charge(),
                record.canonicalSmiles(),
                record.source());
    }

    private void addStructureKindMessage(
            List<ValidationMessage> messages,
            CompoundStructureKind structureKind,
            String preferredName) {
        switch (structureKind) {
            case REPRESENTATIVE_FORMULA_UNIT -> messages.add(ValidationMessage.info(
                    "REPRESENTATIVE_FORMULA_UNIT",
                    "Rendering a representative formula unit for " + preferredName
                            + " because ionic and extended solids do not have a single unique molecular geometry."));
            case ELEMENTAL_ATOM -> messages.add(ValidationMessage.info(
                    "ELEMENTAL_ATOM",
                    "Rendering an isolated atom representation."));
            case ELEMENTAL_MOLECULE -> messages.add(ValidationMessage.info(
                    "ELEMENTAL_MOLECULE",
                    "Rendering the standard elemental molecule for this formula."));
            case DISCRETE_MOLECULE -> {
            }
        }
    }

    private String firstMolBlock(String input) {
        int delimiter = input.indexOf("$$$$");
        return delimiter >= 0 ? input.substring(0, delimiter).trim() : input.trim();
    }

    private String formulaFromContainer(IAtomContainer container) {
        try {
            return MolecularFormulaManipulator.getString(MolecularFormulaManipulator.getMolecularFormula(container));
        } catch (Exception exception) {
            return "";
        }
    }

    private int totalFormalCharge(IAtomContainer container) {
        int total = 0;
        for (IAtom atom : container.atoms()) {
            Integer formalCharge = atom.getFormalCharge();
            if (formalCharge != null) {
                total += formalCharge;
            }
        }
        return total;
    }

    private String toHillFormula(Map<String, Integer> composition) {
        if (composition.isEmpty()) {
            return "";
        }
        Map<String, Integer> ordered = new LinkedHashMap<>();
        if (composition.containsKey("C")) {
            ordered.put("C", composition.get("C"));
            if (composition.containsKey("H")) {
                ordered.put("H", composition.get("H"));
            }
        }
        composition.keySet().stream()
                .filter(symbol -> !ordered.containsKey(symbol))
                .sorted()
                .forEach(symbol -> ordered.put(symbol, composition.get(symbol)));

        StringBuilder builder = new StringBuilder();
        ordered.forEach((symbol, count) -> {
            builder.append(symbol);
            if (count != 1) {
                builder.append(count);
            }
        });
        return builder.toString();
    }

    private String withChargeSuffix(int charge) {
        if (charge == 0) {
            return "";
        }
        int magnitude = Math.abs(charge);
        return magnitude + (charge > 0 ? "+" : "-");
    }

    private StructureResult invalid(StructureInputFormat inputFormat, String message) {
        return new StructureResult(
                StructureResolutionStatus.INVALID,
                inputFormat,
                "",
                "",
                SourceReference.unknown(),
                0,
                0,
                0,
                null,
                List.of(),
                List.of(ValidationMessage.error("STRUCTURE_INPUT_INVALID", message)));
    }
}
