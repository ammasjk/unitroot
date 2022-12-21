package io.datanapis.unitroot;

import io.datanapis.unitroot.data.JGMData;
import io.datanapis.unitroot.data.Probs;
import io.datanapis.unitroot.data.RegressionType;
import io.datanapis.unitroot.data.TestType;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


public class UnitRootEvaluator {
    private static final int NP = 9;

    enum Type {
        NO_INTERCEPT(RegressionType.NC),
        CONSTANT(RegressionType.C),
        CONSTANT_WITH_TREND(RegressionType.CT);

        private RegressionType value;

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
        this(ts, 1, Type.NO_INTERCEPT);
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
        if (this.type == Type.NO_INTERCEPT) {
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
        double pValue = fpval(data, statistic, 2.0, n, NP);

        return pValue;
    }

    public static double criticalValue(double level, int nobs, int niv, TestType testType, RegressionType type) {
        JGMData data = JGMData.getInstance(niv, testType, type);
        double criticalValue = fcrit(data, level, 2.0, nobs, NP);
        return criticalValue;
    }

    private static double fcrit(JGMData data, double level, double precrt, int nobs, int np) {
        if (level < 0.0001 || level > 0.9999)
            throw new IllegalArgumentException("level must be between 0.0001 and 0.9999");

        double[] crits = new double[JGMData.NROWS];
        double[] yvect;
        double[] xmat;
        double[] omega;

        double[] beta = new double[data.getNvar()];

        double diffMin = Double.MAX_VALUE;
        int imin = Probs.PROBS.length;
        for (int i = 0; i < Probs.PROBS.length; i++) {
            double diff = Math.abs(level - Probs.PROBS[i]);
            if (diff < diffMin) {
                diffMin = diff;
                imin = i;
                if (diffMin < 1d - 6)
                    break;
            }
        }

        int nph = np / 2;
        int npTop = JGMData.NROWS - 1 - nph;
        if (imin > nph && imin < npTop) {
            // imin is not too close to the end. Use np points around stat.
            yvect = new double [np];
            xmat = new double [np * 4];
            for (int i = 0; i < np; i++) {
                int ic = imin - nph + i;

                System.arraycopy(data.getBeta(), ic * data.getNvar(), beta, 0, data.getNvar());
                crits[ic] = eval(beta, data.getModel(), data.getNreg(), nobs);

                yvect[i] = crits[ic];
                initialize(xmat, i*4, Probs.CNORM[ic]);
            }

            // form omega matrix
            omega = new double [np * np];
            for (int i = 0; i < np; i++) {
                for (int j = i; j < np; j++) {
                    int ic = imin - nph + i;
                    int jc = imin - nph + j;

                    double top = Probs.PROBS[ic] * (1.0 - Probs.PROBS[jc]);
                    double bot = Probs.PROBS[jc] * (1.0 - Probs.PROBS[ic]);

                    omega[i*np + j] = data.getWght()[ic] * data.getWght()[jc] * Math.sqrt(top/bot);
                }
            }
            fill(omega, np);

            return criticalValueEstimate(xmat, yvect, omega, np, precrt, level);
        } else {
            /* imin is close to one of the ends. Use points from imin +/- nph to end. */
            int np1;
            if (imin < np) {
                np1 = imin + nph;
                if (np1 < 5)
                    np1 = 5;

                yvect = new double [np1];
                xmat = new double [np1 * 4];

                for (int i = 0; i < np1; i++) {
                    System.arraycopy(data.getBeta(), i * data.getNvar(), beta, 0, data.getNvar());
                    crits[i] = eval(beta, data.getModel(), data.getNreg(), nobs);

                    yvect[i] = crits[i];
                    initialize(xmat, i*4, Probs.CNORM[i]);
                }
            } else {
                np1 = JGMData.NROWS - 1 - imin + nph;
                if (np1 < 5)
                    np1 = 5;

                yvect = new double [np1];
                xmat = new double [np1 * 4];

                for (int i = 0; i < np1; i++) {
                    int ic = JGMData.NROWS - 1 - i;

                    System.arraycopy(data.getBeta(), ic * data.getNvar(), beta, 0, data.getNvar());
                    crits[ic] = eval(beta, data.getModel(), data.getNreg(), nobs);

                    yvect[i] = crits[ic];
                    initialize(xmat, i*4, Probs.CNORM[ic]);
                }
            }

            // form omega matrix
            omega = new double [np1 * np1];
            for (int i = 0; i < np1; i++) {
                for (int j = i; j < np1; j++) {
                    if (imin < np) {
                        double top = Probs.PROBS[i] * (1.0 - Probs.PROBS[j]);
                        double bot = Probs.PROBS[j] * (1.0 - Probs.PROBS[i]);
                        omega[i*np1 + j] = data.getWght()[i] * data.getWght()[j] * Math.sqrt(top/bot);
                    } else {
                        /* This is to avoid numerical singularities at the upper end */
                        omega[i*np1 + j] = 0.0;
                        if (i == j)
                            omega[i*np1 + j] = 1.0;
                    }
                }
            }
            fill(omega, np1);

            return criticalValueEstimate(xmat, yvect, omega, np1, precrt, level);
        }
    }

    public static double pValue(double statistic, int nobs, int niv, TestType testType, RegressionType type) {
        JGMData data = JGMData.getInstance(niv, testType, type);
        double pValue = fpval(data, statistic, 2.0, nobs, NP);
        return pValue;
    }

    private static double fpval(JGMData data, double stat, double precrt, int nobs, int np) {
        double[] crits = new double [JGMData.NROWS];
        double[] yvect;
        double[] xmat;
        double[] omega;

        double[] beta = new double [data.getNvar()];

        // first, compute all the estimated critical values
        for (int i = 0; i < JGMData.NROWS; i++) {
            System.arraycopy(data.getBeta(), i * data.getNvar(), beta, 0, data.getNvar());
            crits[i] = eval(beta, data.getModel(), data.getNreg(), nobs);
        }

        // find critical value closest to test statistic
        double diffMin = Double.MAX_VALUE;
        int imin = JGMData.NROWS;
        for (int i = 0; i < JGMData.NROWS; i++) {
            double diff = Math.abs(stat - crits[i]);
            if (diff < diffMin) {
                diffMin = diff;
                imin = i;
            }
        }

        int nph = np / 2;
        int npTop = JGMData.NROWS - 1 - nph;
        if (imin > nph && imin < npTop) {
            // imin is not too close to the end. Use np points around stat.
            yvect = new double [np];
            xmat = new double [np * 4];
            for (int i = 0; i < np; i++) {
                int ic = imin - nph + i;
                yvect[i] = Probs.CNORM[ic];
                initialize(xmat, i*4, crits[ic]);
            }

            // form omega matrix
            omega = new double [np * np];
            for (int i = 0; i < np; i++) {
                for (int j = i; j < np; j++) {
                    int ic = imin - nph + i;
                    int jc = imin - nph + j;

                    double top = Probs.PROBS[ic] * (1.0 - Probs.PROBS[jc]);
                    double bot = Probs.PROBS[jc] * (1.0 - Probs.PROBS[ic]);

                    omega[i*np + j] = data.getWght()[ic] * data.getWght()[jc] * Math.sqrt(top/bot);
                }
            }
            fill(omega, np);

            return pValueEstimate(xmat, yvect, omega, np, precrt, stat);
        } else {
            /* imin is close to one of the ends. Use points from imin +/- nph to end. */
            int np1;
            if (imin < np) {
                np1 = imin + nph;
                if (np1 < 5)
                    np1 = 5;

                yvect = new double [np1];
                xmat = new double [np1 * 4];

                for (int i = 0; i < np1; i++) {
                    yvect[i] = Probs.CNORM[i];
                    initialize(xmat, i*4, crits[i]);
                }
            } else {
                np1 = JGMData.NROWS - 1 - imin + nph;
                if (np1 < 5)
                    np1 = 5;

                yvect = new double [np1];
                xmat = new double [np1 * 4];

                for (int i = 0; i < np1; i++) {
                    int ic = JGMData.NROWS - 1 - i;
                    yvect[i] = Probs.CNORM[ic];
                    initialize(xmat, i*4, crits[ic]);
                }
            }

            // form omega matrix
            omega = new double [np1 * np1];
            for (int i = 0; i < np1; i++) {
                for (int j = i; j < np1; j++) {
                    if (imin < np) {
                        double top = Probs.PROBS[i] * (1.0 - Probs.PROBS[j]);
                        double bot = Probs.PROBS[j] * (1.0 - Probs.PROBS[i]);
                        omega[i*np1 + j] = data.getWght()[i] * data.getWght()[j] * Math.sqrt(top/bot);
                    } else {
                        /* This is to avoid numerical singularities at the upper end */
                        omega[i*np1 + j] = 0.0;
                        if (i == j)
                            omega[i*np1 + j] = 1.0;
                    }
                }
            }
            fill(omega, np1);

            return pValueEstimate(xmat, yvect, omega, np1, precrt, stat);
        }
    }

    private static void initialize(double[] xmat, int start, double value) {
        xmat[start] = 1.0;
        xmat[start + 1] = value;
        xmat[start + 2] = xmat[start + 1] * value;
        xmat[start + 3] = xmat[start + 2] * value;
    }

    private static double[] reshape(double[] xmatInput, int rows) {
        double[] xmatOutput = new double [rows * 3];
        for (int i = 0; i < rows; i++) {
            System.arraycopy(xmatInput, i * 4, xmatOutput, i * 3, 3);
        }

        return xmatOutput;
    }

    private static void fill(double[] omega, int np) {
        for (int i = 0; i < np; i++) {
            for (int j = i; j < np; j++) {
                omega[j*np + i] = omega[i*np + j];
            }
        }
    }

    private static double criticalValueEstimate(double[] xmat, double[] yvect, double[] omega, int np, double precrt, double size) {
        GLS gls = new GLS(xmat, yvect, omega, np, 4, true);
        gls.gls();

        // check to see if gamma(4) is needed
        double[] gamma = gls.getBeta();
        double sd4 = Math.sqrt((gls.getSsrt()/(np-4)) * gls.getXomx()[3*4 + 3]);
        double tTest = Math.abs(gamma[3]) / sd4;

        double criticalValue;
        if (tTest > precrt) {
            double anorm = innorz(size);
            criticalValue = gamma[0] + gamma[1]*anorm + gamma[2]*Math.pow(anorm,2) + gamma[3]*Math.pow(anorm,3);
        } else {
            /* adjust xmat from np x 4 to np x 3 */
            xmat = reshape(xmat, np);
            gls = new GLS(xmat, yvect, omega, np, 3, false);
            gls.gls();
            gamma = gls.getBeta();

            double anorm = innorz(size);
            criticalValue = gamma[0] + gamma[1]*anorm + gamma[2]*Math.pow(anorm,2);
        }

        return criticalValue;
    }

    private static double pValueEstimate(double[] xmat, double[] yvect, double[] omega, int np, double precrt, double stat) {
        GLS gls = new GLS(xmat, yvect, omega, np, 4, true);
        gls.gls();

        double pValue;

        // check to see if gamma(4) is needed
        double[] gamma = gls.getBeta();
        double sd4 = Math.sqrt((gls.getSsrt()/(np-4)) * gls.getXomx()[3*4 + 3]);
        double tTest = Math.abs(gamma[3]) / sd4;
        if (tTest > precrt) {
            double crfit = gamma[0] + gamma[1]*stat + gamma[2]*Math.pow(stat,2) + gamma[3]*Math.pow(stat,3);
            pValue = ddnor(crfit);
        } else {
            /* adjust xmat from np x 4 to np x 3 */
            xmat = reshape(xmat, np);
            gls = new GLS(xmat, yvect, omega, np, 3, false);
            gls.gls();
            gamma = gls.getBeta();
            double crfit = gamma[0] + gamma[1]*stat + gamma[2]*Math.pow(stat,2);
            pValue = ddnor(crfit);
        }

        return pValue;
    }

    private static class GLS {
        private final double[] xmat;
        private final double[] yvect;
        private final double[] omega;
        private final int nobs;
        private final int nvar;
        private final boolean invert;

        private final double[] xomx;
        private final double[] xomy;
        private final double[] beta;
        private final double[] fits;
        private final double[] resid;
        private double ssr;
        private double ssrt;

        private GLS(double[] xmat, double[] yvect, double[] omega, int nobs, int nvar, boolean invert) {
            this.xmat = xmat;       /* nobs x nvar */
            this.yvect = yvect;     /* nobs */
            this.omega = omega;     /* nobs x nobs */
            this.nobs = nobs;       /* usually around 9, but varies near the ends */
            this.nvar = nvar;       /* usually 3 or 4 */
            this.invert = invert;

            this.xomx = new double [nvar * nvar];
            this.xomy = new double [nvar];
            this.beta = new double [nvar];
            this.fits = new double [nobs];
            this.resid = new double [nobs];
        }

        public double[] getXomx() {
            return xomx;
        }

        public double[] getXomy() {
            return xomy;
        }

        public double[] getBeta() {
            return beta;
        }

        public double[] getFits() {
            return fits;
        }

        public double[] getResid() {
            return resid;
        }

        public double getSsr() {
            return ssr;
        }

        public double getSsrt() {
            return ssrt;
        }

        /**
         * Copyright (c) James G. MacKinnon, 1995
         * Subroutine to do GLS estimation the obvious way
         * Use only when sample size is small (nobs <= 50)
         * 1995-1-3
         */
        public void gls() {
            if (invert)
                cholx(omega, nobs);

            // form xomx matrix and xomy vector
            for (int j = 0; j < nvar; j++) {
                xomy[j] = 0.0;
                for (int l = j; l < nvar; l++) {
                    xomx[j*nvar + l] = 0.0;
                }
            }

            for (int i = 0; i < nobs; i++) {
                for (int k = 0; k < nobs; k++) {
                    for (int j = 0; j < nvar; j++) {
                        /* xomy(j) = xomy(j) + xmat(i,j)*omega(k,i)*yvect(k) */
                        xomy[j] = xomy[j] + xmat[i*nvar + j] * omega[k*nobs + i] * yvect[k];
                        for (int l = j; l < nvar; l++) {
                            /* xomx(j,l) = xomx(j,l) + xmat(i,j)*omega(k,i)*xmat(k,l) */
                            xomx[j*nvar + l] = xomx[j*nvar + l] + xmat[i*nvar + j] * omega[k*nobs + i] * xmat[k*nvar + l];
                        }
                    }
                }
            }

            for (int j = 0; j < nvar; j++) {
                for (int l = j; l < nvar; l++) {
                    /* xomx(l,j) = xomx(j,l) */
                    xomx[l*nvar + j] = xomx[j*nvar + l];
                }
            }

            // invert xomx matrix
            cholx(xomx, nvar);

            // now form estimates of beta.
            for (int i = 0; i < nvar; i++) {
                beta[i] = 0.0;
                for (int j = 0; j < nvar; j++) {
                    /* beta(i) = beta(i) + xomx(i,j) * xomy(j) */
                    beta[i] = beta[i] + xomx[i*nvar + j] * xomy[j];
                }
            }

            // find ssr, fitted values, and residuals
            ssr = 0.0;
            for (int i = 0; i < nobs; i++) {
                fits[i] = 0.0;
                for (int j = 0; j < nvar; j++) {
                    /* fits(i) = fits(i) + xmat(i,j)*beta(j) */
                    fits[i] = fits[i] + xmat[i*nvar + j] * beta[j];
                }
                /* resid(i) = yvect(i) - fits(i) */
                resid[i] = yvect[i] - fits[i];
                /* ssr = ssr + resid(i)**2 */
                ssr = ssr + resid[i] * resid[i];
            }

            // find ssr from transformed regression
            for (int i = 0; i < nobs; i++) {
                for (int k = 0; k < nobs; k++) {
                    /* ssrt = ssrt + resid(i)*omega(k,i)*resid(k) */
                    ssrt = ssrt + resid[i] * omega[k*nobs + i] * resid[k];
                }
            }
        }

        /**
         * Copyright (c) James G. MacKinnon, 1993
         * This routine uses the Cholesky decomposition to invert a real symmetric matrix.
         *
         * The original Fortran implementation has been adapted to Java. Specifically, arrays are
         * allocated to the required size and values are returned rather than passed out using out-parameters.
         * Two-dimensional arrays are also represented as a single contiguous array
         *
         * @param amat the matrix to be inverted
         * @param n the dimensions of the matrix i.e., n x n
         */
        private static int cholx(double[] amat, int n) {
            double ooa = 1.0;

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    if (i > 0) {
                        for (int k = 0; k < i; k++) {
                            /* amat(i,j) = amat(i,j) - amat(k,i)*amat(k,j) */
                            amat[i*n + j] = amat[i*n + j] - amat[k*n + i] * amat[k*n + j];
                        }
                    } else if (amat[i*n + i] <= 0.0) {
                        return i;
                    }

                    if (i == j) {
                        /* amat(i,i) = dsqrt(amat(i,i)) */
                        amat[i*n + i] = Math.sqrt(amat[i*n + i]);
                    } else {
                        if (j == i+1)
                            ooa = 1.0/amat[i*n + i];
                        amat[i*n + j] = amat[i*n + j] * ooa;
                    }
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    ooa = 1.0/amat[j*n + j];
                    double t;
                    if (i >= j) {
                        t = 1.0;
                    } else {
                        t = 0.0;
                        for (int k = i; k < j; k++) {
                            t = t - amat[i*n + k] * amat[k*n + j];
                        }
                    }
                    amat[i*n + j] = t * ooa;
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    double t = 0.0;
                    for (int k = j; k < n; k++) {
                        /* t = t + amat(i,k)*amat(j,k) */
                        t = t + amat[i*n + k] * amat[j*n + k];
                    }
                    amat[i*n + j] = t;
                    amat[j*n + i] = t;
                }
            }

            return 0;
        }
    }

    /**
     * Copyright (c) James G. MacKinnon, 1995
     * Routine to evaluate response surface for specified betas and sample size.
     *
     * Adapted from Fortran to Java.
     *
     * @param beta - coefficients from the Urc* data classes, beta[4]
     * @param model - defines the cardinality of beta (i.e. either 3 or 4), from the Urc* data classes
     * @param nreg - number of variables (?)
     * @param nobs - number of observations
     * @return the critical value (?)
     */
    private static double eval(double[] beta, int model, int nreg, int nobs) {
        if (nobs == 0) {
            return beta[0];
        }

        if (model == 2) {
            double onobs = 1.0/nobs;
            /* beta(1) + beta(2)*onobs + beta(3)*onobs**2 */
            return beta[0] + beta[1] * onobs + beta[2] * Math.pow(onobs, 2);
        } else if (model == 3) {
            double onobs = 1.0/nobs;
            /* beta(1) + beta(2)*onobs + beta(3)*onobs**2 + beta(4)*onobs**3 */
            return beta[0] + beta[1] * onobs + beta[2] * Math.pow(onobs, 2) + beta[3] * Math.pow(onobs, 3);
        } else if (model == 4) {
            double onobs = 1.0/(nobs - nreg);
            /* beta(1) + beta(2)*onobs + beta(3)*onobs**2 */
            return beta[0] + beta[1] * onobs + beta[2] * Math.pow(onobs, 2);
        } else if (model == 5) {
            double onobs = 1.0/(nobs - nreg);
            /* beta(1) + beta(2)*onobs + beta(3)*onobs**2 + beta(4)*onobs**3 */
            return beta[0] + beta[1] * onobs + beta[2] * Math.pow(onobs, 2) + beta[3] * Math.pow(onobs, 3);
        }

        throw new IllegalArgumentException(String.format("Invalid value for model [%d]", model));
    }

    /**
     * Copyright (c) James G. MacKinnon, 1993
     * Routine to evaluate cumulative normal distribution
     * Written originally in late 1970's
     * Modified 1993 to avoid changing the argument
     *
     * This subroutine uses Cody's method to evaluate the cumulative
     * normal distribution. It is probably accurate to 19 or 20
     * significant digits. It was written in 1977, based on the Cody
     * article referred to in the documentation for IMSL subroutine mdnor.
     *
     * @param ystar
     */
    private static double ddnor(double ystar) {
        final double[] p = {
                -6.58749161529837803157e-04,
                -1.60837851487422766278e-02,
                -1.25781726111229246204e-01,
                -3.60344899949804439429e-01,
                -3.05326634961232344035e-01,
                -1.63153871373020978498e-02
        };
        final double[] q = {
                2.33520497626869185443e-03,
                6.05183413124413191178e-02,
                5.27905102951428412248e-01,
                1.87295284992346047209e00,
                2.56852019228982242072e00
        };
        final double[] a = {
                1.23033935479799725272e03,
                2.05107837782607146532e03,
                1.71204761263407058314e03,
                8.81952221241769090411e02,
                2.98635138197400131132e02,
                6.61191906371416294775e01,
                8.88314979438837594118e00,
                5.64188496988670089180e-01,
                2.15311535474403846343e-08
        };
        final double[] b = {
                1.23033935480374942043e03,
                3.43936767414372163696e03,
                4.36261909014324715820e03,
                3.29079923573345962678e03,
                1.62138957456669018874e03,
                5.37181101862009857509e02,
                1.17693950891312499305e02,
                1.57449261107098347253e01
        };
        final double[] c = {
                3.209377589138469472562e03,
                3.774852376853020208137e02,
                1.138641541510501556495e02,
                3.161123743870565596947e00,
                1.857777061846031526730e-01
        };
        final double[] d = {
                2.844236833439170622273e03,
                1.282616526077372275645e03,
                2.440246379344441733056e02,
                2.360129095234412093499e01
        };
        final double orpi = .5641895835477562869483e0;
        final double root2 = .70710678118654752440083e0;

        int isw = 1;

        if (ystar < -16.0)
            ystar = -16.0;
        else if (ystar > 16.0)
            ystar = 16.0;

        double x = -ystar * root2;

        if (!(x > 0.0 || x < 0.0)) {
            /* x == 0.0 */
            return 0.5;
        }

        if (x < 0.0) {
            x = -x;
            isw = -1;
        }

        if (x > 4.0) {
            double x2 = x*x;
            double xm2 = 1.0/x2;
            double xm4 = xm2*xm2;
            double xm6 = xm4*xm2;
            double xm8 = xm4*xm4;
            double xm10 = xm6*xm4;
            double top = p[0] + p[1]*xm2 + p[2]*xm4 + p[3]*xm6 + p[4]*xm8 + p[5]*xm10;
            double bot = q[0] + q[1]*xm2 + q[2]*xm4 + q[3]*xm6 + q[4]*xm8 + xm10;
            double crap = orpi + top/(bot*x2);
            double erfc = Math.exp(-x2) * crap/x;

            if (isw == -1)
                erfc = 2.0 - erfc;

            return erfc * .5;
        } else if (x < 0.477) {
            double x2 = x*x;
            double x4 = x2*x2;
            double x6 = x4*x2;
            double x8 = x4*x4;
            double top = c[0] + c[1]*x2 + c[2]*x4 + c[3]*x6 + c[4]*x8;
            double bot = d[0] + d[1]*x2 + d[2]*x4 + d[3]*x6 + x8;
            double erf = x*top/bot;

            erf = erf * isw;
            double erfc = 1.0 - erf;

            return erfc * .5;
        } else {
            double x2 = x*x;
            double x3 = x2*x;
            double x4 = x2*x2;
            double x5 = x3*x2;
            double x6 = x3*x3;
            double x7 = x3*x4;
            double x8 = x4*x4;
            double top = a[0] + a[1]*x + a[2]*x2 + a[3]*x3 + a[4]*x4 + a[5]*x5 + a[6]*x6 + a[7]*x7 + a[8]*x8;
            double bot = b[0] + b[1]*x + b[2]*x2 + b[3]*x3 + b[4]*x4 + b[5]*x5 + b[6]*x6 + b[7]*x7 + x8;
            double erfc = Math.exp(-x2)*top/bot;

            if (isw == -1)
                erfc = 2.0 - erfc;

            return erfc * .5;
        }
    }

    /**
     * Copyright (c) James G. MacKinnon, 1995
     * Inverse normal routine that adjusts crude result twice.
     * It seems to be accurate to about 14 digits.
     * Crude result is taken from Abramowitz & Stegun (1968)
     * It should have abs. error < 4.5 * 10^-4
     *
     * @param prob
     */
    private static double innorz(double prob) {
        double c0 = 2.515517e0;
        double c1 = 0.802853e0;
        double c2 = 0.010328e0;
        double d1 = 1.432788e0;
        double d2 = 0.189269e0;
        double d3 = 0.001308e0;
        double constant = .398942280401432678e0;

        if (prob < 0.0 || prob > 1.0) {
            throw new IllegalArgumentException(String.format("prob [%.2f] must be in the range [0.0, 1.0]", prob));
        }

        double pr = prob;
        if (prob > 0.5)
            pr = 1.0 - prob;

        double t = Math.sqrt(-2 * Math.log(pr));
        double t2 = Math.pow(t, 2);
        double t3 = Math.pow(t, 3);
        double anorm = t - (c0 + c1*t + c2*t2) / (1 + d1*t + d2*t2 + d3*t3);

        //
        // now correct crude result by direct method
        //
        double prob2 = ddnor(anorm);
        double pr2 = 1.0 - prob2;
        t = Math.sqrt(-2 * Math.log(pr2));
        t2 = Math.pow(t, 2);
        t3 = Math.pow(t, 3);
        double anorm2 = t - (c0 + c1*t + c2*t2) / (1 + d1*t + d2*t2 + d3*t3);
        anorm = anorm + anorm - anorm2;
        if (prob < 0.5)
            anorm = -anorm;

        //
        // now correct better result, using Taylor series approximation
        //
        prob2 = ddnor(anorm);
        double error = prob2 - prob;
        double dens = constant * Math.exp(-0.5 * Math.pow(anorm, 2));
        return anorm - error/dens;
    }
}
