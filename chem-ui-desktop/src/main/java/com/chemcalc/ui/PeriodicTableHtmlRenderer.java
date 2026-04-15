package com.chemcalc.ui;

import com.chemcalc.domain.ElementRecord;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.List;

final class PeriodicTableHtmlRenderer {
    private static final String TEMPLATE = """
            <html>
              <head>
                <meta charset="utf-8">
                <style>__CSS__</style>
              </head>
              <body>
                <div class="page">
                  <div class="hero">
                    <div>
                      <div class="eyebrow">Chemistry Super App</div>
                      <h1>Periodic Table Explorer</h1>
                      <p>Dense desktop explorer with native search and a synchronized detail panel.</p>
                    </div>
                    <div id="stats" class="stats">Loading element data...</div>
                  </div>
                  <div id="legend" class="legend"></div>
                  <div id="empty-state" class="empty-state hidden">No elements match the current filters.</div>
                  <div id="main-grid" class="grid"></div>
                  <section class="series-section">
                    <div class="section-heading">Lanthanides</div>
                    <div id="lanthanides" class="series-grid"></div>
                  </section>
                  <section class="series-section">
                    <div class="section-heading">Actinides</div>
                    <div id="actinides" class="series-grid"></div>
                  </section>
                </div>
                <script>__JS__</script>
              </body>
            </html>
            """;

    private static final String CSS = readResource("/ui/periodic-table.css");
    private static final String SCRIPT = readResource("/ui/periodic-table.js");

    private PeriodicTableHtmlRenderer() {
    }

    static String renderShell() {
        return TEMPLATE.replace("__CSS__", CSS).replace("__JS__", SCRIPT);
    }

    static String renderUpdateScript(List<ElementRecord> elements, String selectedSymbol) {
        return "window.renderTable(" + toJson(elements) + ", " + quote(selectedSymbol) + ");";
    }

    private static String toJson(List<ElementRecord> elements) {
        StringBuilder builder = new StringBuilder("[");
        for (int index = 0; index < elements.size(); index++) {
            if (index > 0) {
                builder.append(',');
            }
            ElementRecord element = elements.get(index);
            builder.append('{')
                    .append("\"atomicNumber\":").append(element.atomicNumber()).append(',')
                    .append("\"symbol\":").append(quote(element.symbol())).append(',')
                    .append("\"name\":").append(quote(element.name())).append(',')
                    .append("\"group\":").append(element.group() == null ? "null" : element.group()).append(',')
                    .append("\"period\":").append(element.period()).append(',')
                    .append("\"block\":").append(quote(element.block())).append(',')
                    .append("\"category\":").append(quote(element.category())).append(',')
                    .append("\"standardState\":").append(quote(element.standardState().isBlank() ? "Unknown" : element.standardState())).append(',')
                    .append("\"atomicWeight\":").append(quote(element.property("standardAtomicWeight")
                            .map(value -> value.value() + " " + value.units())
                            .orElse("Unavailable")))
                    .append('}');
        }
        return builder.append(']').toString();
    }

    private static String quote(String value) {
        if (value == null) {
            return "null";
        }
        return "\"" + escape(value) + "\"";
    }

    private static String escape(String value) {
        return value
                .replace("\\", "\\\\")
                .replace("\"", "\\\"")
                .replace("\r", "\\r")
                .replace("\n", "\\n");
    }

    private static String readResource(String path) {
        try (InputStream input = PeriodicTableHtmlRenderer.class.getResourceAsStream(path)) {
            if (input == null) {
                throw new IllegalStateException("Missing UI resource: " + path);
            }
            return new String(input.readAllBytes(), StandardCharsets.UTF_8);
        } catch (IOException exception) {
            throw new IllegalStateException("Failed to load UI resource: " + path, exception);
        }
    }
}
