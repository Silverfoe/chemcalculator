package com.chemcalc.engine;

import com.chemcalc.domain.BalancedReaction;
import com.chemcalc.domain.Reaction;

public interface BalancerService {
    BalancedReaction balance(Reaction reaction);
}
