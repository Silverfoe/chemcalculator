package com.chemcalc.engine;

import java.math.BigInteger;

final class Rational implements Comparable<Rational> {
    static final Rational ZERO = new Rational(BigInteger.ZERO, BigInteger.ONE);
    static final Rational ONE = new Rational(BigInteger.ONE, BigInteger.ONE);

    private final BigInteger numerator;
    private final BigInteger denominator;

    Rational(long numerator) {
        this(BigInteger.valueOf(numerator), BigInteger.ONE);
    }

    Rational(BigInteger numerator, BigInteger denominator) {
        if (denominator.signum() == 0) {
            throw new ArithmeticException("Denominator must not be zero");
        }
        if (denominator.signum() < 0) {
            numerator = numerator.negate();
            denominator = denominator.negate();
        }
        BigInteger gcd = numerator.gcd(denominator);
        if (!gcd.equals(BigInteger.ZERO)) {
            numerator = numerator.divide(gcd);
            denominator = denominator.divide(gcd);
        }
        this.numerator = numerator;
        this.denominator = denominator;
    }

    static Rational of(long value) {
        return new Rational(value);
    }

    Rational add(Rational other) {
        return new Rational(
                numerator.multiply(other.denominator).add(other.numerator.multiply(denominator)),
                denominator.multiply(other.denominator));
    }

    Rational subtract(Rational other) {
        return new Rational(
                numerator.multiply(other.denominator).subtract(other.numerator.multiply(denominator)),
                denominator.multiply(other.denominator));
    }

    Rational multiply(Rational other) {
        return new Rational(numerator.multiply(other.numerator), denominator.multiply(other.denominator));
    }

    Rational divide(Rational other) {
        return new Rational(numerator.multiply(other.denominator), denominator.multiply(other.numerator));
    }

    Rational negate() {
        return new Rational(numerator.negate(), denominator);
    }

    boolean isZero() {
        return numerator.equals(BigInteger.ZERO);
    }

    BigInteger numerator() {
        return numerator;
    }

    BigInteger denominator() {
        return denominator;
    }

    @Override
    public int compareTo(Rational other) {
        return numerator.multiply(other.denominator).compareTo(other.numerator.multiply(denominator));
    }

    @Override
    public String toString() {
        return denominator.equals(BigInteger.ONE) ? numerator.toString() : numerator + "/" + denominator;
    }
}
