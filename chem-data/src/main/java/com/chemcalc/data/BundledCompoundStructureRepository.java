package com.chemcalc.data;

import com.chemcalc.domain.CompoundStructureRecord;
import com.chemcalc.domain.CompoundStructureRepository;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

public final class BundledCompoundStructureRepository implements CompoundStructureRepository {
    private final List<CompoundStructureRecord> records;
    private final Map<String, List<CompoundStructureRecord>> byFormulaAndCharge;

    public BundledCompoundStructureRepository() {
        this.records = CompoundStructureCatalog.allRecords();
        Map<String, List<CompoundStructureRecord>> grouped = new LinkedHashMap<>();
        for (CompoundStructureRecord record : records) {
            grouped.computeIfAbsent(key(record.canonicalFormula(), record.charge()), ignored -> new ArrayList<>())
                    .add(record);
        }
        Map<String, List<CompoundStructureRecord>> immutable = new LinkedHashMap<>();
        grouped.forEach((groupKey, value) -> immutable.put(groupKey, List.copyOf(value)));
        this.byFormulaAndCharge = Map.copyOf(immutable);
    }

    @Override
    public List<CompoundStructureRecord> findAll() {
        return records;
    }

    @Override
    public List<CompoundStructureRecord> findByFormulaAndCharge(String canonicalFormula, int charge) {
        return byFormulaAndCharge.getOrDefault(key(canonicalFormula, charge), List.of());
    }

    private String key(String canonicalFormula, int charge) {
        return canonicalFormula.trim().toUpperCase(Locale.ROOT) + "#" + charge;
    }
}
