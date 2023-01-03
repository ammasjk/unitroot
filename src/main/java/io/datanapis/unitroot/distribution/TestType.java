/**
 * Author: Jayakumar Muthukumarasamy
 */
package io.datanapis.unitroot.distribution;

public enum TestType {
    TAU(1),
    Z(2);

    private final int value;

    TestType(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }
}
