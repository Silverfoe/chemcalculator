package com.chemcalc.domain;

import java.util.List;

public record Reaction(
        List<Species> reactants,
        List<Species> products,
        ReactionArrow arrow,
        List<ValidationMessage> messages) {

    public Reaction {
        reactants = List.copyOf(reactants == null ? List.of() : reactants);
        products = List.copyOf(products == null ? List.of() : products);
        arrow = arrow == null ? ReactionArrow.FORWARD : arrow;
        messages = List.copyOf(messages == null ? List.of() : messages);
    }

    public String format() {
        return reactants.stream().map(Species::displayFormula).reduce((a, b) -> a + " + " + b).orElse("")
                + " " + arrow.token() + " "
                + products.stream().map(Species::displayFormula).reduce((a, b) -> a + " + " + b).orElse("");
    }
}
