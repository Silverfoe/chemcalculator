package com.chemcalc.domain;

import java.util.List;

public record PercentCompositionResult(
        List<PercentCompositionEntry> entries,
        List<ValidationMessage> messages) {

    public PercentCompositionResult {
        entries = List.copyOf(entries == null ? List.of() : entries);
        messages = List.copyOf(messages == null ? List.of() : messages);
    }
}
