package io.datanapis.unitroot;

/**
 * Author: Jayakumar Muthukumarasamy
 *
 * R compatible methods
 */
public class RUtils {
    /**
     * Return lagged difference of x. The first value will be x[lag] - x[0], the second value will be x[1+lag] - x[1]
     * and the i'th value will be x[i+lag] - x[i].
     *
     * @param x the values to lag
     * @param lag the lag
     * @return the lagged values. The return array will be of size x.length - lag
     */
    public static double[] diff(double[] x, int lag) {
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
    public static double[] embed(double[] y, int dimension) {
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
    public static double[] column(double[] y, int column, int length) {
        double[] c = new double [length];
        System.arraycopy(y, column * length, c, 0, length);
        return c;
    }

    public static void copyColumn(double[] src, int i, double[] dest, int j, int length) {
        System.arraycopy(src, i * length, dest, j * length, length);
    }

    public static void copyColumnTransposed(double[] src, int i, double[] dest, int j, int nCols, int length) {
        for (int k = 0; k < length; k++) {
            dest[k * nCols + j] = src[i * length + k];
        }
    }
}
