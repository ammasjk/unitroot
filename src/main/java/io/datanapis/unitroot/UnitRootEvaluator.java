package io.datanapis.unitroot;

import io.datanapis.unitroot.distribution.JGMData;
import io.datanapis.unitroot.distribution.RegressionType;
import io.datanapis.unitroot.distribution.TestType;
import io.datanapis.unitroot.distribution.UrcDist;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import java.util.Arrays;


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
     * This is a port of the function unitrootTest from the R package fUnitRoot (file UnitRootTests.R)
     *
     * @return an instance of UnitRootResult containing the coefficients (0th index is the constant term if Type != NC),
     * pValue and the statistic
     */
    public UnitRootResult test() {
        int lag = this.lag + 1;

        /* Δy(t) = λy(t−1) + μ + βt + α(1)Δy(t−1) + ... + α(k)Δy(t−k) + ε(t) */

        double[] y = RUtils.diff(x, 1);    /* get differences between consecutive terms */
        int n = y.length;
        double[] z = RUtils.embed(y, lag);     /* See R function embed, logic from UnitrootTests.R from fUnitRoots */

        /* the lagged original values become the first predictor */
        double[] yLag1 = new double [n - lag + 1];
        System.arraycopy(x, lag - 1, yLag1, 0, yLag1.length);

        assert yLag1.length == (z.length / lag);

        int nCols = lag + 1;   /* +1 for the response variable - see OLSMultipleLinearRegression for details */
        if (this.type == Type.CONSTANT_WITH_TREND) {
            /* Need to add a lag term as a predictor */
            nCols += 1;
        } else if (this.type == Type.CONSTANT_WITH_TREND_SQUARED) {
            /* Need to add both a lag and lag**2 as predictors */
            nCols += 2;
        }

        double[] values = new double [nCols * yLag1.length];
        int j = 0;
        /* response values for regression, 0th column becomes the response variable, see equation above */
        RUtils.copyColumnTransposed(z, 0, values, j, nCols, yLag1.length); ++j;
        /* first predictor, i.e. lagged y, see equation above */
        RUtils.copyColumnTransposed(yLag1, 0, values, j, nCols, yLag1.length); ++j;
        if (this.type == Type.CONSTANT_WITH_TREND || this.type == Type.CONSTANT_WITH_TREND_SQUARED) {
            double[] tt = new double [yLag1.length];
            for (int i = 0; i < yLag1.length; i++) {
                tt[i] = lag + i;
            }
            /* lag becomes a predictor */
            RUtils.copyColumnTransposed(tt, 0, values, j, nCols, yLag1.length); ++j;
            if (this.type == Type.CONSTANT_WITH_TREND_SQUARED) {
                double[] ttSquared = new double [yLag1.length];
                for (int i = 0; i < yLag1.length; i++) {
                    ttSquared[i] = Math.pow(lag + i, 2.0);
                }
                /* lag**2 becomes a predictor as well */
                RUtils.copyColumnTransposed(ttSquared, 0, values, j, nCols, yLag1.length); ++j;
            }
        }
        /* columns, 1 through lag are used as predictors as well, 0'th column is the response variable */
        for (int i = 1; i < lag; i++) {
            RUtils.copyColumnTransposed(z, i, values, j, nCols, yLag1.length); ++j;
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
