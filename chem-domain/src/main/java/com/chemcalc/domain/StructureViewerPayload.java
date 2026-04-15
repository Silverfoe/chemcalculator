package com.chemcalc.domain;

import java.util.Objects;

public record StructureViewerPayload(
        String modelFormat,
        String modelData,
        boolean has3dCoordinates,
        String svgData) {

    public StructureViewerPayload {
        modelFormat = Objects.requireNonNullElse(modelFormat, "").trim();
        modelData = Objects.requireNonNullElse(modelData, "").trim();
        svgData = Objects.requireNonNullElse(svgData, "").trim();
    }
}
