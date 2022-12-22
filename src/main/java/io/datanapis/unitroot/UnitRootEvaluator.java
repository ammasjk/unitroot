package io.datanapis.unitroot;

import io.datanapis.unitroot.distribution.JGMData;
import io.datanapis.unitroot.distribution.RegressionType;
import io.datanapis.unitroot.distribution.TestType;
import io.datanapis.unitroot.distribution.UrcDist;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


public class UnitRootEvaluator {
    private static final int NP = 9;

    enum Type {
        NO_CONSTANT(RegressionType.NC),
        CONSTANT(RegressionType.C),
        CONSTANT_WITH_TREND(RegressionType.CT),
        CONSTANT_WITH_TREND_SQUARED(RegressionType.CTT);

        private final RegressionType value;

        Type(RegressionType value) {
            this.value = value;
        }

        public RegressionType getValue() {
            return this.value;
        }
    }

    private final double[] x;
    private final int lags;
    private final Type type;

    public UnitRootEvaluator(double[] ts) {
        this(ts, 1, Type.NO_CONSTANT);
    }

    public UnitRootEvaluator(double[] ts, int lags, Type type) {
        this.x = new double [ts.length];
        System.arraycopy(ts, 0, x, 0, ts.length);

        this.lags = lags;
        this.type = type;
    }

    public int getLags() {
        return lags;
    }

    public Type getType() {
        return type;
    }

    private static double[] diff(double[] x, int lags) {
        int length = x.length;
        if (lags <= 0 || lags >= length)
            throw new IllegalArgumentException();

        double[] y = new double [x.length - lags];
        for (int i = 0; i < y.length; i++) {
            y[i] = x[i+lags] - x[i];
        }

        return y;
    }

    /* Columns are stored sequentially */
    private static double[] embed(double[] y, int dimension) {
        int length = y.length;
        if (dimension <= 0 || dimension >= length)
            throw new IllegalArgumentException();

        int size = length - dimension + 1;
        double[] value = new double [dimension * size];
        for (int i = 0; i < dimension; i++) {
            System.arraycopy(y, dimension - i - 1, value, i * size, size);
        }

        return value;
    }

    private static double[] column(double[] y, int column, int length) {
        double[] c = new double [length];
        System.arraycopy(y, column * length, c, 0, length);
        return c;
    }

    public double pValue(TestType testType) {
        int lags = this.lags + 1;

        double[] y = diff(x, 1);
        int n = y.length;
        double[] z = embed(y, lags);

        double[] yDiff = column(z, 0, z.length / lags);
        double[] yLag1 = new double [n - lags + 1];
        System.arraycopy(x, lags - 1, yLag1, 0, yLag1.length);

        assert yLag1.length == (z.length / lags);

        double[] values = new double [(lags + 1) * yLag1.length];
        for (int i = 0; i < yLag1.length; i++) {
            for (int j = 0; j < lags + 1; j++) {
                if (j == 0) {
                    values[i * (lags + 1)] = yDiff[i];
                } else if (j == 1) {
                    values[i * (lags + 1) + 1] = yLag1[i];
                } else {
                    values[i * (lags + 1) + j] = z[(j - 1) * yLag1.length + i];
                }
            }
        }

        OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
        int index;
        if (this.type == Type.NO_CONSTANT) {
            index = 0;
            regression.setNoIntercept(true);
        } else {
            index = 1;
            regression.setNoIntercept(false);
        }
        regression.newSampleData(values, yLag1.length, lags);

        double[] estimates = regression.estimateRegressionParameters();
        double[] stdErrors = regression.estimateRegressionParametersStandardErrors();

        double statistic = estimates[index] / stdErrors[index];
        JGMData data = JGMData.getInstance(1, testType, type.getValue());
        double pValue = UrcDist.fpval(data, statistic, 2.0, n, NP);

        return pValue;
    }
}
