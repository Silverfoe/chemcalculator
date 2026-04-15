package com.chemcalc.engine;

import com.chemcalc.domain.StructureRequest;
import com.chemcalc.domain.StructureResult;

public interface StructureService {
    boolean available();

    StructureResult resolve(StructureRequest request);
}
