package com.chemcalc.domain;

public enum Phase {
    SOLID("s"),
    LIQUID("l"),
    GAS("g"),
    AQUEOUS("aq"),
    UNKNOWN("");

    private final String token;

    Phase(String token) {
        this.token = token;
    }

    public String token() {
        return token;
    }

    public static Phase fromToken(String token) {
        if (token == null) {
            return UNKNOWN;
        }
        return switch (token.trim().toLowerCase()) {
            case "s" -> SOLID;
            case "l" -> LIQUID;
            case "g" -> GAS;
            case "aq" -> AQUEOUS;
            default -> UNKNOWN;
        };
    }
}
