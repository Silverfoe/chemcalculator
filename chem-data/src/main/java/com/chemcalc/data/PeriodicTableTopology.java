package com.chemcalc.data;

import java.util.Set;

final class PeriodicTableTopology {
    private static final Set<String> ALKALI = Set.of("Li", "Na", "K", "Rb", "Cs", "Fr");
    private static final Set<String> ALKALINE = Set.of("Be", "Mg", "Ca", "Sr", "Ba", "Ra");
    private static final Set<String> METALLOIDS = Set.of("B", "Si", "Ge", "As", "Sb", "Te", "Po");
    private static final Set<String> HALOGENS = Set.of("F", "Cl", "Br", "I", "At", "Ts");
    private static final Set<String> NOBLE_GASES = Set.of("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og");
    private static final Set<String> POST_TRANSITION = Set.of("Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi", "Nh", "Fl", "Mc", "Lv");
    private static final Set<String> OTHER_NONMETALS = Set.of("H", "C", "N", "O", "P", "S", "Se");

    record Layout(Integer group, int period, String block, String category) {
    }

    private PeriodicTableTopology() {
    }

    static Layout forElement(int atomicNumber, String symbol) {
        Integer group = null;
        int period;

        if (atomicNumber == 1) {
            group = 1;
            period = 1;
        } else if (atomicNumber == 2) {
            group = 18;
            period = 1;
        } else if (atomicNumber <= 10) {
            period = 2;
            group = atomicNumber <= 4 ? atomicNumber - 2 : atomicNumber + 8;
        } else if (atomicNumber <= 18) {
            period = 3;
            group = atomicNumber <= 12 ? atomicNumber - 10 : atomicNumber;
        } else if (atomicNumber <= 36) {
            period = 4;
            group = atomicNumber - 18;
        } else if (atomicNumber <= 54) {
            period = 5;
            group = atomicNumber - 36;
        } else if (atomicNumber <= 56) {
            period = 6;
            group = atomicNumber - 54;
        } else if (atomicNumber <= 71) {
            period = 6;
            group = 3;
        } else if (atomicNumber <= 86) {
            period = 6;
            group = atomicNumber - 68;
        } else if (atomicNumber <= 88) {
            period = 7;
            group = atomicNumber - 86;
        } else if (atomicNumber <= 103) {
            period = 7;
            group = 3;
        } else {
            period = 7;
            group = atomicNumber - 100;
        }

        String block;
        if (atomicNumber >= 57 && atomicNumber <= 71 || atomicNumber >= 89 && atomicNumber <= 103) {
            block = "f";
        } else if ("He".equals(symbol) || group != null && group <= 2) {
            block = "s";
        } else if (group != null && group <= 12) {
            block = "d";
        } else {
            block = "p";
        }

        String category;
        if (atomicNumber >= 57 && atomicNumber <= 71) {
            category = "Lanthanide";
        } else if (atomicNumber >= 89 && atomicNumber <= 103) {
            category = "Actinide";
        } else if (NOBLE_GASES.contains(symbol)) {
            category = "Noble Gas";
        } else if (HALOGENS.contains(symbol)) {
            category = "Halogen";
        } else if (ALKALI.contains(symbol)) {
            category = "Alkali Metal";
        } else if (ALKALINE.contains(symbol)) {
            category = "Alkaline Earth Metal";
        } else if (METALLOIDS.contains(symbol)) {
            category = "Metalloid";
        } else if (POST_TRANSITION.contains(symbol)) {
            category = "Post-Transition Metal";
        } else if (OTHER_NONMETALS.contains(symbol)) {
            category = "Nonmetal";
        } else if (group != null && group >= 3 && group <= 12) {
            category = "Transition Metal";
        } else {
            category = "Metal";
        }

        return new Layout(group, period, block, category);
    }
}
