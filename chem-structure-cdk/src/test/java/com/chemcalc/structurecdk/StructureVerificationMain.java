package com.chemcalc.structurecdk;

import com.chemcalc.domain.StructureInputFormat;
import com.chemcalc.domain.StructureRequest;
import com.chemcalc.domain.StructureResolutionStatus;
import com.chemcalc.domain.StructureResult;

public final class StructureVerificationMain {
    private StructureVerificationMain() {
    }

    public static void main(String[] args) {
        CdkStructureService service = new CdkStructureService();

        verifyFormulaResolved(service);
        verifyFormulaNormalized(service);
        verifyFormulaElement(service);
        verifyElementalMolecule(service);
        verifyFormulaAmbiguous(service);
        verifyFormulaRepresentative(service);
        verifyFormulaUnsupported(service);
        verifySmiles(service);
        verifyMolfile(service);
        verifyInchi(service);

        System.out.println("CDK structure verification passed.");
    }

    private static void verifyFormulaResolved(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "H2O"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected formula H2O to resolve");
        require(result.hasRenderablePayload(), "Expected formula H2O to produce a viewer payload");
        require("Water".equals(result.preferredName()), "Expected H2O to resolve to water");
    }

    private static void verifyFormulaNormalized(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "HCl"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected HCl to resolve");
        require("Hydrogen chloride".equals(result.preferredName()), "Expected HCl to resolve to hydrogen chloride");
        require("ClH".equals(result.canonicalFormula()), "Expected HCl to normalize to ClH");
    }

    private static void verifyFormulaElement(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "Na"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected Na to resolve as an element");
        require("Sodium".equals(result.preferredName()), "Expected Na to resolve to sodium");
        require(result.hasRenderablePayload(), "Expected Na to produce a viewer payload");
    }

    private static void verifyElementalMolecule(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "Cl2"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected Cl2 to resolve");
        require("Chlorine".equals(result.preferredName()), "Expected Cl2 to resolve to chlorine");
    }

    private static void verifyFormulaAmbiguous(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "C2H6O"));
        require(result.status() == StructureResolutionStatus.AMBIGUOUS, "Expected C2H6O to be ambiguous");
        require(result.candidates().size() == 2, "Expected C2H6O to expose two bundled candidates");
    }

    private static void verifyFormulaRepresentative(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "NaCl"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected NaCl to resolve from the offline catalog");
        require("Sodium chloride".equals(result.preferredName()), "Expected NaCl to resolve to sodium chloride");
        require(result.messages().stream().anyMatch(message -> "REPRESENTATIVE_FORMULA_UNIT".equals(message.code())),
                "Expected NaCl to explain that it is using a representative formula unit");
    }

    private static void verifyFormulaUnsupported(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.FORMULA, "ClH4NaO"));
        require(result.status() == StructureResolutionStatus.UNSUPPORTED,
                "Expected unsupported formulas outside the offline catalog to remain unsupported");
    }

    private static void verifySmiles(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.SMILES, "CCO"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected SMILES to resolve");
        require(result.hasRenderablePayload(), "Expected SMILES to produce a viewer payload");
        require(result.viewerPayload().modelData().contains("M  END"), "Expected SMILES payload to export a Molfile block");
    }

    private static void verifyMolfile(CdkStructureService service) {
        StructureResult source = service.resolve(new StructureRequest(StructureInputFormat.SMILES, "O"));
        StructureResult roundTrip = service.resolve(new StructureRequest(StructureInputFormat.MOLFILE, source.viewerPayload().modelData()));
        require(roundTrip.status() == StructureResolutionStatus.RESOLVED, "Expected Molfile to round-trip");
        require(roundTrip.hasRenderablePayload(), "Expected Molfile round-trip to produce a viewer payload");
    }

    private static void verifyInchi(CdkStructureService service) {
        StructureResult result = service.resolve(new StructureRequest(StructureInputFormat.INCHI, "InChI=1S/H2O/h1H2"));
        require(result.status() == StructureResolutionStatus.RESOLVED, "Expected InChI to resolve");
        require(result.hasRenderablePayload(), "Expected InChI to produce a viewer payload");
    }

    private static void require(boolean condition, String message) {
        if (!condition) {
            throw new IllegalStateException(message);
        }
    }
}
