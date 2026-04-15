package com.chemcalc.domain;

import java.util.Objects;

public record SourceReference(String id, String title, String revision, String url) {
    public SourceReference {
        id = Objects.requireNonNullElse(id, "").trim();
        title = Objects.requireNonNullElse(title, "").trim();
        revision = Objects.requireNonNullElse(revision, "").trim();
        url = Objects.requireNonNullElse(url, "").trim();
    }

    public static SourceReference unknown() {
        return new SourceReference("UNKNOWN", "Unknown Source", "", "");
    }
}
