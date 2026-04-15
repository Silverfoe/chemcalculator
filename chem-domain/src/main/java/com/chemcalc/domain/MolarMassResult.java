package com.chemcalc.domain;

import java.util.List;

public record MolarMassResult(
        double totalMass,
        List<MolarMassComponent> components,
        List<ValidationMessage> messages) {

    public MolarMassResult {
        components = List.copyOf(components == null ? List.of() : components);
        messages = List.copyOf(messages == null ? List.of() : messages);
    }

    public boolean successful() {
        return messages.stream().noneMatch(message -> message.severity() == ValidationSeverity.ERROR);
    }
}
