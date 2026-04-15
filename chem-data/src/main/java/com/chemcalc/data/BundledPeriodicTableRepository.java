package com.chemcalc.data;

import com.chemcalc.domain.ElementRecord;
import com.chemcalc.domain.PeriodicTableRepository;
import com.chemcalc.domain.PropertyKind;
import com.chemcalc.domain.PropertyValue;
import com.chemcalc.domain.SourceReference;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Optional;

public final class BundledPeriodicTableRepository implements PeriodicTableRepository {
    private final List<ElementRecord> elements;
    private final Map<String, ElementRecord> bySymbol;
    private final Map<Integer, ElementRecord> byAtomicNumber;

    public BundledPeriodicTableRepository() {
        this.elements = loadElements();
        Map<String, ElementRecord> symbolMap = new LinkedHashMap<>();
        Map<Integer, ElementRecord> numberMap = new LinkedHashMap<>();
        for (ElementRecord element : elements) {
            symbolMap.put(element.symbol().toUpperCase(Locale.ROOT), element);
            numberMap.put(element.atomicNumber(), element);
        }
        this.bySymbol = Map.copyOf(symbolMap);
        this.byAtomicNumber = Map.copyOf(numberMap);
    }

    @Override
    public List<ElementRecord> findAll() {
        return elements;
    }

    @Override
    public Optional<ElementRecord> findByAtomicNumber(int atomicNumber) {
        return Optional.ofNullable(byAtomicNumber.get(atomicNumber));
    }

    @Override
    public Optional<ElementRecord> findBySymbol(String symbol) {
        if (symbol == null) {
            return Optional.empty();
        }
        return Optional.ofNullable(bySymbol.get(symbol.trim().toUpperCase(Locale.ROOT)));
    }

    @Override
    public List<ElementRecord> search(String query) {
        if (query == null || query.isBlank()) {
            return elements;
        }
        String normalized = query.trim().toLowerCase(Locale.ROOT);
        return elements.stream()
                .filter(element -> element.symbol().toLowerCase(Locale.ROOT).contains(normalized)
                        || element.name().toLowerCase(Locale.ROOT).contains(normalized)
                        || Integer.toString(element.atomicNumber()).equals(normalized))
                .toList();
    }

    private List<ElementRecord> loadElements() {
        SourceRegistry sourceRegistry = SourceRegistry.load();
        List<ElementRecord> loaded = new ArrayList<>();
        try (InputStream input = BundledPeriodicTableRepository.class.getResourceAsStream("/data/elements.psv")) {
            if (input == null) {
                throw new IllegalStateException("Bundled element snapshot not found");
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8))) {
                String header = reader.readLine();
                if (header == null) {
                    return List.of();
                }
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.isBlank() || line.startsWith("#")) {
                        continue;
                    }
                    String[] parts = line.split("\\|", -1);
                    if (parts.length < 7) {
                        continue;
                    }

                    int atomicNumber = Integer.parseInt(parts[0]);
                    String symbol = parts[1];
                    String name = parts[2];
                    String atomicWeight = parts[3];
                    String atomicWeightSource = parts[4];
                    String standardState = parts[5];
                    String yearDiscovered = parts[6];

                    PeriodicTableTopology.Layout layout = PeriodicTableTopology.forElement(atomicNumber, symbol);
                    Map<String, PropertyValue<String>> properties = new LinkedHashMap<>();
                    if (!atomicWeight.isBlank()) {
                        SourceReference source = sourceRegistry.lookup(atomicWeightSource);
                        properties.put(
                                "standardAtomicWeight",
                                PropertyValue.of(atomicWeight, "g/mol", PropertyKind.TEXT, source));
                    }
                    loaded.add(new ElementRecord(
                            atomicNumber,
                            symbol,
                            name,
                            layout.group(),
                            layout.period(),
                            layout.block(),
                            layout.category(),
                            standardState,
                            yearDiscovered,
                            properties,
                            List.of()));
                }
            }
        } catch (IOException exception) {
            throw new IllegalStateException("Failed to load bundled periodic table data", exception);
        }
        return List.copyOf(loaded);
    }
}
