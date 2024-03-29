package io.datanapis.unitroot.distribution;

import java.io.IOException;

/**
 * Author: Jayakumar Muthukumarasamy
 *
 * Class containing distribution data for a given set of parameters.
 *
 * James G. MacKinnon, "Numerical distribution functions for unit root and cointegration tests,"
 * Journal of Applied Econometrics, 11, 1996, 601-618.
 */
public class MackinnonData {
    public static final int NROWS = 221;

    private final String name;
    private final int nz;
    private final int nreg;
    private final int model;
    private final int minSize;
    private final double[] beta;
    private final double[] wght;

    protected MackinnonData(String name, int nz, int nreg, int model, int minSize, double[] beta) {
        this.name = name;
        this.nz = nz;
        this.nreg = nreg;
        this.model = model;
        this.minSize = minSize;

        int nvar = (model == 2 || model == 4) ? 3 : 4;
        assert beta.length == NROWS * (nvar + 1) : "Size mismatch";

        this.beta = new double [nvar * NROWS];
        this.wght = new double [NROWS];

        for (int i = 0; i < NROWS; i++) {
            System.arraycopy(beta, i * (nvar + 1), this.beta, i * nvar, nvar);
            wght[i] = beta[i * (nvar + 1) + nvar];
        }
    }

    public String getName() {
        return name;
    }

    public int getNz() {
        return nz;
    }

    public int getNreg() {
        return nreg;
    }

    public int getModel() {
        return model;
    }

    public int getNvar() {
        return (model == 2 || model == 4) ? 3 : 4;
    }

    public int getMinSize() {
        return minSize;
    }

    public double[] getBeta() {
        return beta;
    }

    public double[] getWght() {
        return wght;
    }

    public static MackinnonData getInstance(int niv, TestType tt, RegressionType rt) {
        try {
            DataReader dataReader = new DataReader();
            MackinnonData data = dataReader.read(niv, tt, rt);
            return data;
        } catch (IOException e) {
            return null;
        }
    }
}
