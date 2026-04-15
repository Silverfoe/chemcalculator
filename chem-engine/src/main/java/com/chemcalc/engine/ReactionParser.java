package com.chemcalc.engine;

import com.chemcalc.domain.Reaction;

public interface ReactionParser {
    Reaction parse(String input);
}
