/**
 * Author: Jayakumar Muthukumarasamy
 */
package io.datanapis.unitroot.distribution;

public enum RegressionType {
    NC(1),      /* No intercept or trend */
    C(2),       /* Constant */
    CT(3),      /* Constant and trend */
    CTT(4);     /* Constant, trend and trend-squared */

    private final int value;

    RegressionType(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }
}
