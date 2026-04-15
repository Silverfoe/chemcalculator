package com.chemcalc.domain;

public sealed interface FormulaPart permits ElementNode, GroupNode {
    String toFormulaString();
}
