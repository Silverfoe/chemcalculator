package com.chemcalc.domain;

public enum ReactionArrow {
    FORWARD("->"),
    REVERSIBLE("<->");

    private final String token;

    ReactionArrow(String token) {
        this.token = token;
    }

    public String token() {
        return token;
    }
}
