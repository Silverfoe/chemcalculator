package com.chemcalc.engine;

import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.MolarMassComponent;
import com.chemcalc.domain.MolarMassResult;
import com.chemcalc.domain.PercentCompositionEntry;
import com.chemcalc.domain.PercentCompositionResult;
import com.chemcalc.domain.PeriodicTableRepository;
import com.chemcalc.domain.PropertyValue;
import com.chemcalc.domain.Species;
import com.chemcalc.domain.ValidationMessage;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public final class DefaultStoichiometryService implements StoichiometryService {
    private final PeriodicTableRepository periodicTableRepository;

    public DefaultStoichiometryService(PeriodicTableRepository periodicTableRepository) {
        this.periodicTableRepository = periodicTableRepository;
    }

    @Override
    public MolarMassResult molarMass(Species species) {
        List<MolarMassComponent> components = new ArrayList<>();
        List<ValidationMessage> messages = new ArrayList<>();
        double total = 0.0;
        for (String symbol : species.composition().keySet().stream().sorted().toList()) {
            int count = species.composition().get(symbol);
            ElementRecord element = periodicTableRepository.findBySymbol(symbol).orElse(null);
            if (element == null) {
                messages.add(ValidationMessage.error("UNKNOWN_ELEMENT", "Unknown element: " + symbol));
                continue;
            }
            PropertyValue<String> weight = element.property("standardAtomicWeight").orElse(null);
            if (weight == null || !weight.isPresent()) {
                messages.add(ValidationMessage.error("MISSING_ATOMIC_WEIGHT", "Missing atomic weight for " + symbol));
                continue;
            }
            double atomicWeight = parseNumeric(weight.value());
            double contribution = atomicWeight * count;
            total += contribution;
            components.add(new MolarMassComponent(symbol, count, atomicWeight, contribution));
        }
        components.sort(Comparator.comparing(MolarMassComponent::symbol));
        return new MolarMassResult(total, components, messages);
    }

    @Override
    public PercentCompositionResult percentComposition(Species species) {
        MolarMassResult molarMassResult = molarMass(species);
        if (!molarMassResult.successful() || molarMassResult.totalMass() == 0.0) {
            return new PercentCompositionResult(List.of(), molarMassResult.messages());
        }
        List<PercentCompositionEntry> entries = new ArrayList<>();
        for (MolarMassComponent component : molarMassResult.components()) {
            entries.add(new PercentCompositionEntry(
                    component.symbol(),
                    component.contribution(),
                    component.contribution() / molarMassResult.totalMass() * 100.0));
        }
        return new PercentCompositionResult(entries, List.of());
    }

    private double parseNumeric(String value) {
        String normalized = value.trim();
        if (normalized.startsWith("[") && normalized.endsWith("]")) {
            String[] bounds = normalized.substring(1, normalized.length() - 1).split(",");
            if (bounds.length == 2) {
                return (Double.parseDouble(bounds[0].trim()) + Double.parseDouble(bounds[1].trim())) / 2.0;
            }
        }
        return Double.parseDouble(normalized);
    }
}
