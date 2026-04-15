package com.chemcalc.domain;

import java.util.Objects;
import java.util.Optional;

public record PropertyValue<T>(
        T value,
        String units,
        PropertyKind kind,
        SourceReference source,
        boolean precise) {

    public PropertyValue {
        units = Objects.requireNonNullElse(units, "").trim();
        kind = kind == null ? PropertyKind.UNKNOWN : kind;
        source = source == null ? SourceReference.unknown() : source;
    }

    public static <T> PropertyValue<T> of(T value, String units, PropertyKind kind, SourceReference source) {
        return new PropertyValue<>(value, units, kind, source, true);
    }

    public static <T> PropertyValue<T> approximate(T value, String units, PropertyKind kind, SourceReference source) {
        return new PropertyValue<>(value, units, kind, source, false);
    }

    public Optional<T> optionalValue() {
        return Optional.ofNullable(value);
    }

    public boolean isPresent() {
        return value != null && !value.toString().isBlank();
    }
}
