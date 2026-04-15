package com.chemcalc.engine;

import com.chemcalc.domain.FormulaAst;
import com.chemcalc.domain.Species;

public interface FormulaParser {
    FormulaAst parseFormula(String input);

    Species parseSpecies(String input);
}
