package com.chemcalc.domain;

import java.util.List;
import java.util.Objects;

public record StructureResult(
        StructureResolutionStatus status,
        StructureInputFormat inputFormat,
        String preferredName,
        String canonicalFormula,
        SourceReference source,
        int atomCount,
        int bondCount,
        int formalCharge,
        StructureViewerPayload viewerPayload,
        List<StructureCandidate> candidates,
        List<ValidationMessage> messages) {

    public StructureResult {
        status = status == null ? StructureResolutionStatus.INVALID : status;
        inputFormat = inputFormat == null ? StructureInputFormat.FORMULA : inputFormat;
        preferredName = Objects.requireNonNullElse(preferredName, "").trim();
        canonicalFormula = Objects.requireNonNullElse(canonicalFormula, "").trim();
        source = source == null ? SourceReference.unknown() : source;
        candidates = List.copyOf(candidates == null ? List.of() : candidates);
        messages = List.copyOf(messages == null ? List.of() : messages);
    }

    public boolean hasRenderablePayload() {
        return viewerPayload != null && !viewerPayload.modelData().isBlank();
    }
}
