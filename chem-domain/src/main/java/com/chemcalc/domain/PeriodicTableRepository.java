package com.chemcalc.domain;

import java.util.List;
import java.util.Optional;

public interface PeriodicTableRepository {
    List<ElementRecord> findAll();

    Optional<ElementRecord> findByAtomicNumber(int atomicNumber);

    Optional<ElementRecord> findBySymbol(String symbol);

    List<ElementRecord> search(String query);
}
