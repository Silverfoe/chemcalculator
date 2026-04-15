package com.chemcalc.domain;

import java.util.List;

public interface CompoundStructureRepository {
    List<CompoundStructureRecord> findAll();

    List<CompoundStructureRecord> findByFormulaAndCharge(String canonicalFormula, int charge);
}
