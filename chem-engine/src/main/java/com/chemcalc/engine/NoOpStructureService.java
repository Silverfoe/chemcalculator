package com.chemcalc.engine;

import com.chemcalc.domain.SourceReference;
import com.chemcalc.domain.StructureRequest;
import com.chemcalc.domain.StructureResolutionStatus;
import com.chemcalc.domain.StructureResult;
import com.chemcalc.domain.ValidationMessage;
import java.util.List;

public final class NoOpStructureService implements StructureService {
    @Override
    public boolean available() {
        return false;
    }

    @Override
    public StructureResult resolve(StructureRequest request) {
        return new StructureResult(
                StructureResolutionStatus.UNSUPPORTED,
                request.format(),
                "",
                "",
                SourceReference.unknown(),
                0,
                0,
                0,
                null,
                List.of(),
                List.of(ValidationMessage.warning(
                        "STRUCTURE_BACKEND_UNAVAILABLE",
                        "No structure backend is configured for this application.")));
    }
}
