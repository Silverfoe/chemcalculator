package com.chemcalc.engine;

import com.chemcalc.domain.MolarMassResult;
import com.chemcalc.domain.PercentCompositionResult;
import com.chemcalc.domain.Species;

public interface StoichiometryService {
    MolarMassResult molarMass(Species species);

    PercentCompositionResult percentComposition(Species species);
}
