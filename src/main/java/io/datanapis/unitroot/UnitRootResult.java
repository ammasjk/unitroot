package io.datanapis.unitroot;

/**
 * Author: Jayakumar Muthukumarasamy
 *
 * Result of testing for unit root
 */
public class UnitRootResult {
    private final UnitRootEvaluator.Type type;
    private final int lag;
    private final double statistic;
    private final double pValueTau;
    private final double pValueZ;
    private final double[] coefficients;

    UnitRootResult(UnitRootEvaluator.Type type, int lag, double statistic, double pValueTau, double pValueZ, double[] coefficients) {
        this.type = type;
        this.lag = lag;
        this.statistic = statistic;
        this.coefficients = coefficients;
        this.pValueTau = pValueTau;
        this.pValueZ = pValueZ;
    }

    public UnitRootEvaluator.Type getType() {
        return type;
    }

    public int getLag() {
        return this.lag;
    }

    public double[] getCoefficients() {
        return this.coefficients;
    }

    public double pValueTau() {
        return this.pValueTau;
    }

    public double pValueZ() {
        return this.pValueZ;
    }

    public double lambda() {
        int index = (type == UnitRootEvaluator.Type.NO_CONSTANT) ? 0 : 1;
        return this.coefficients[index];
    }

    @Override
    public String toString() {
        if (type == UnitRootEvaluator.Type.NO_CONSTANT) {
            return String.format("(%s, b1=%f, p=(%f, %f), t=%f)",
                    type.name(), coefficients[0], pValueTau, pValueZ, statistic);
        } else {
            return String.format("(%s, b0=%f, b1=%f, p=(%f, %f), t=%f)",
                    type.name(), coefficients[0], coefficients[1], pValueTau, pValueZ, statistic);
        }
    }
}
