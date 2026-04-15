package com.chemcalc.engine;

import com.chemcalc.domain.NomenclatureResult;
import com.chemcalc.domain.Species;

public interface NomenclatureService {
    NomenclatureResult name(Species species);
}
