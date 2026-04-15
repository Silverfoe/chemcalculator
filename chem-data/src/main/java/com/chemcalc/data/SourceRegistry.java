package com.chemcalc.data;

import com.chemcalc.domain.SourceReference;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

final class SourceRegistry {
    private final Map<String, SourceReference> sources;

    private SourceRegistry(Map<String, SourceReference> sources) {
        this.sources = Map.copyOf(sources);
    }

    static SourceRegistry load() {
        Map<String, SourceReference> loaded = new LinkedHashMap<>();
        try (InputStream input = SourceRegistry.class.getResourceAsStream("/data/sources.psv")) {
            if (input == null) {
                return new SourceRegistry(Map.of());
            }
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(input, StandardCharsets.UTF_8))) {
                String header = reader.readLine();
                if (header == null) {
                    return new SourceRegistry(Map.of());
                }
                String line;
                while ((line = reader.readLine()) != null) {
                    if (line.isBlank() || line.startsWith("#")) {
                        continue;
                    }
                    String[] parts = line.split("\\|", -1);
                    if (parts.length < 4) {
                        continue;
                    }
                    loaded.put(parts[0], new SourceReference(parts[0], parts[1], parts[2], parts[3]));
                }
            }
        } catch (IOException exception) {
            throw new IllegalStateException("Failed to load source registry", exception);
        }
        return new SourceRegistry(loaded);
    }

    SourceReference lookup(String id) {
        if (id == null || id.isBlank()) {
            return SourceReference.unknown();
        }
        return sources.getOrDefault(id.trim(), SourceReference.unknown());
    }
}
