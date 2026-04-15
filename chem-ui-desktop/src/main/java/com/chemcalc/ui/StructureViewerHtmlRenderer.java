package com.chemcalc.ui;

import com.chemcalc.domain.StructureResult;
import com.chemcalc.domain.ValidationMessage;

final class StructureViewerHtmlRenderer {
    private StructureViewerHtmlRenderer() {
    }

    static String renderStructureScript(StructureResult result) {
        StringBuilder builder = new StringBuilder("window.renderStructure(");
        builder.append("{")
                .append("\"preferredName\":").append(quote(result.preferredName())).append(',')
                .append("\"canonicalFormula\":").append(quote(result.canonicalFormula())).append(',')
                .append("\"inputFormat\":").append(quote(result.inputFormat().name())).append(',')
                .append("\"has3dCoordinates\":").append(result.viewerPayload().has3dCoordinates()).append(',')
                .append("\"modelFormat\":").append(quote(result.viewerPayload().modelFormat())).append(',')
                .append("\"modelData\":").append(quote(result.viewerPayload().modelData())).append(',')
                .append("\"svgData\":").append(quote(result.viewerPayload().svgData())).append(',')
                .append("\"messages\":").append(messagesToJson(result))
                .append("});");
        return builder.toString();
    }

    static String renderMessageScript(String title, String detail) {
        return "window.showViewerMessage(" + quote(title) + ", " + quote(detail) + ");";
    }

    private static String messagesToJson(StructureResult result) {
        StringBuilder builder = new StringBuilder("[");
        for (int index = 0; index < result.messages().size(); index++) {
            ValidationMessage message = result.messages().get(index);
            if (index > 0) {
                builder.append(',');
            }
            builder.append('{')
                    .append("\"severity\":").append(quote(message.severity().name())).append(',')
                    .append("\"code\":").append(quote(message.code())).append(',')
                    .append("\"message\":").append(quote(message.message()))
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

}
