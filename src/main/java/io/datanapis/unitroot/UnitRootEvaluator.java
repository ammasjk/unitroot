package io.datanapis.unitroot;

import io.datanapis.unitroot.distribution.JGMData;
import io.datanapis.unitroot.distribution.RegressionType;
import io.datanapis.unitroot.distribution.TestType;
import io.datanapis.unitroot.distribution.UrcDist;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


public class UnitRootEvaluator {
    private static final int NP = 9;

    public enum Type {
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
    private final int lag;
    private final Type type;

    public UnitRootEvaluator(double[] ts) {
        this(ts, 1, Type.NO_CONSTANT);
    }

    public UnitRootEvaluator(double[] ts, int lag, Type type) {
        this.x = new double [ts.length];
        System.arraycopy(ts, 0, x, 0, ts.length);

        this.lag = lag;
        this.type = type;
    }

    public int getLag() {
        return lag;
    }

    public Type getType() {
        return type;
    }

    /**
     * Return lagged difference of x. The first value will be x[lag] - x[0], the second value will be x[1+lag] - x[1]
     * and the i'th value will be x[i+lag] - x[i].
     *
     * @param x the values to lag
     * @param lag the lag
     * @return the lagged values. The return array will be of size x.length - lag
     */
    private static double[] diff(double[] x, int lag) {
        int length = x.length;
        if (lag <= 0 || lag >= length)
            throw new IllegalArgumentException();

        double[] y = new double [x.length - lag];
        for (int i = 0; i < y.length; i++) {
            y[i] = x[i+lag] - x[i];
        }

        return y;
    }

    /**
     * See R function embed for documentation
     *
     * @param y the values to embed
     * @param dimension the dimensions
     * @return a two-dimensional matrix consisting of 'dimension' columns each with lagged values of y starting from
     * dimension - 1 (0th column) till 0 (last column). The column values are stored sequentially. Each column is of
     * length y.length - dimension + 1 (since the 0th column will only have y.length - dimension + 1 valid values)
     */
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

    /**
     * Return values for 'column' from a two-dimensional matrix y. Each column is expected to be stored
     * sequentially within y and is of size 'length'
     *
     * @param y the two-dimensional matrix with column values stored sequentially
     * @param column the column whose value needs to be returned
     * @param length the length of each column
     * @return the values for column
     */
    private static double[] column(double[] y, int column, int length) {
        double[] c = new double [length];
        System.arraycopy(y, column * length, c, 0, length);
        return c;
    }

    /**
     * This is a port of the function unitrootTest from the R package fUnitRoot (file UnitRootTests.R)
     *
     * @return an instance of UnitRootResult containing the coefficients (0th index is the constant term if Type != NC),
     * pValue and the statistic
     */
    public UnitRootResult test() {
        int lag = this.lag + 1;

        /* Δy(t) = λy(t−1) + μ + βt + α(1)Δy(t−1) + ... + α(k)Δy(t−k) + ε(t) */

        double[] y = diff(x, 1);    /* get differences between consecutive terms */
        int n = y.length;
        double[] z = embed(y, lag);     /* See R function embed, logic from UnitrootTests.R from fUnitRoots */

        /* 0th column becomes the response variable */
        double[] yDiff = column(z, 0, z.length / lag);

        /* the lagged original values become the first predictor */
        double[] yLag1 = new double [n - lag + 1];
        System.arraycopy(x, lag - 1, yLag1, 0, yLag1.length);

        assert yLag1.length == (z.length / lag);

        int offset = 1;
        int nCols = lag + 1;   /* +1 for the response variable - see OLSMultipleLinearRegression for details */
        if (this.type == Type.CONSTANT_WITH_TREND) {
            /* Need to add a lag term as a predictor */
            nCols += 1;
            offset += 1;
        } else if (this.type == Type.CONSTANT_WITH_TREND_SQUARED) {
            /* Need to add both a lag and lag**2 as predictors */
            nCols += 2;
            offset += 2;
        }

        int tt = lag;

        double[] values = new double [nCols * yLag1.length];
        for (int i = 0; i < yLag1.length; i++) {
            for (int j = 0; j < nCols; j++) {
                switch (j) {
                    case 0:
                        /* response values for regression, i.e. difference in y values, see equation above */
                        values[i * nCols] = yDiff[i];
                        continue;
                    case 1:
                        /* first predictor, i.e. lagged y, see equation above */
                        values[i * nCols + 1] = yLag1[i];
                        continue;
                    case 2:
                        if (this.type == Type.CONSTANT_WITH_TREND) {
                            /* lag becomes a predictor */
                            values[i * nCols + j] = tt;
                            ++tt;
                            continue;
                        } else if (this.type == Type.CONSTANT_WITH_TREND_SQUARED) {
                            /* both lag and lag**2 become predictors */
                            values[i * nCols + j] = tt;
                            j++;
                            values[i * nCols + j] = tt * tt;
                            ++tt;
                            continue;
                        }
                    default:
                        /* columns, 1 through lag are used as predictors as well */
                        values[i * nCols + j] = z[(j - offset) * yLag1.length + i];
                }
            }
        }

        try {
            /* regress the response against the predictors */
            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
            int index;
            if (this.type == Type.NO_CONSTANT) {
                index = 0;
                regression.setNoIntercept(true);
            } else {
                index = 1;
                regression.setNoIntercept(false);
            }
            regression.newSampleData(values, yLag1.length, nCols - 1);

            // Compute the t-statistic of the coefficient of y(t-1) in the equation above. Compare the t-statistic
            // with the critical value as defined by the JGM surface, and return the p-value. If the p-value is
            // small enough, we can conclude that the residuals are stationary i.e. the time series are cointegrated
            double[] coefficients = regression.estimateRegressionParameters();
            double[] stdErrors = regression.estimateRegressionParametersStandardErrors();

            double pValueTau = -1.0, pValueZ = -1.0;
            double statistic = coefficients[index] / stdErrors[index];
            JGMData data;
            data = JGMData.getInstance(1, TestType.TAU, type.getValue());
            if (data != null) {
                pValueTau = UrcDist.fpval(data, statistic, 2.0, n, NP);
            }
            data = JGMData.getInstance(1, TestType.Z, type.getValue());
            if (data != null) {
                pValueZ = UrcDist.fpval(data, statistic, 2.0, n, NP);
            }

            return new UnitRootResult(this.type, this.lag, statistic, pValueTau, pValueZ, coefficients);
        } catch (Exception e) {
            throw e;
        }
    }
}
