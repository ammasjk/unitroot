package io.datanapis.unitroot;

import io.datanapis.unitroot.data.RegressionType;
import io.datanapis.unitroot.data.TestType;
import org.junit.Test;

import java.util.function.Function;

public class UnitRootTest {
    /*
     * This file contains the annual data used in Section 7 of the paper.
     * The five columns following the year contain, from left to right,
     *
     * pi_c (inflation rate based on the CPI; see footnote 1 of the paper)
     * pi_y (inflation rate based on the GDP deflator; see footnote 1)
     * r_s  (short term interest rate, from CANSIM B14001)
     * r_m  (medium term interest rate, from CANSIM B14010)
     * r_l  (long term interest rate, from CANSIM B14003)
     */
    public static class JGMTestData {
        private final int year;
        private final double piC;
        private final double piY;
        private final double rS;
        private final double rM;
        private final double rL;

        public JGMTestData(int year, double piC, double piY, double rS, double rM, double rL) {
            this.year = year;
            this.piC = piC;
            this.piY = piY;
            this.rS = rS;
            this.rM = rM;
            this.rL = rL;
        }

        public double piC() {
            return this.piC;
        }

        public double piY() {
            return this.piY;
        }

        public double rS() {
            return this.rS;
        }

        public double rM() {
            return this.rM;
        }

        public double rL() {
            return this.rL;
        }
    }

    private final static JGMTestData[] jgmTestData = new JGMTestData[] {
            new JGMTestData(1952, 2.46347, 4.49336, 1.06583, 3.24167, 3.53917),
            new JGMTestData(1953, -0.93133, -0.24884, 1.69083, 3.44667, 3.76917),
            new JGMTestData(1954, 0.66060, 1.44679, 1.46583, 2.66833, 3.25500),
            new JGMTestData(1955, 0.077402, 0.76162, 1.55333, 2.79083, 3.18917),
            new JGMTestData(1956, 1.42176, 3.04330, 2.90250, 3.75750, 3.60583),
            new JGMTestData(1957, 3.15460, 2.31477, 3.77750, 4.56500, 4.12500),
            new JGMTestData(1958, 2.48276, 1.36357, 2.29250, 3.46917, 4.11500),
            new JGMTestData(1959, 1.18302, 2.07219, 4.80500, 4.93833, 5.04917),
            new JGMTestData(1960, 1.23963, 1.21386, 3.32417, 4.51917, 5.18917),
            new JGMTestData(1961, 1.01561, 0.46074, 2.83417, 4.37500, 5.05833),
            new JGMTestData(1962, 1.10879, 1.37372, 4.01250, 4.60000, 5.10083),
            new JGMTestData(1963, 1.70829, 2.04184, 3.57250, 4.48000, 5.09083),
            new JGMTestData(1964, 1.81273, 2.58770, 3.74250, 4.72083, 5.18417),
            new JGMTestData(1965, 2.49684, 3.31031, 3.96750, 4.89917, 5.20167),
            new JGMTestData(1966, 3.60085, 4.72365, 4.99833, 5.54833, 5.68500),
            new JGMTestData(1967, 3.50590, 4.14178, 4.58667, 5.64000, 5.89833),
            new JGMTestData(1968, 3.96972, 3.57138, 6.24500, 6.67750, 6.72500),
            new JGMTestData(1969, 4.43129, 4.41562, 7.14667, 7.66333, 7.55833),
            new JGMTestData(1970, 3.27988, 4.55531, 6.11750, 7.10500, 7.96500),
            new JGMTestData(1971, 2.81035, 3.10383, 3.61667, 5.55500, 6.94500),
            new JGMTestData(1972, 4.69832, 5.45023, 3.54583, 6.25667, 7.22500),
            new JGMTestData(1973, 7.33239, 8.41061, 5.38750, 6.98417, 7.55000),
            new JGMTestData(1974, 10.36205, 13.43793, 7.78417, 8.12333, 8.87000),
            new JGMTestData(1975, 10.20320, 9.43372, 7.36583, 7.71917, 8.99333),
            new JGMTestData(1976, 7.22462, 8.38706, 8.89167, 8.35000, 9.22917),
            new JGMTestData(1977, 7.71755, 6.09521, 7.35250, 7.90250, 8.69333),
            new JGMTestData(1978, 8.54662, 5.85012, 8.58417, 8.99917, 9.23583),
            new JGMTestData(1979, 8.75210, 9.52905, 11.57417, 10.41917, 10.17917),
            new JGMTestData(1980, 9.68625, 10.09435, 12.68083, 12.37083, 12.33833),
            new JGMTestData(1981, 11.73003, 10.27568, 17.77750, 15.67583, 14.98917),
            new JGMTestData(1982, 10.25870, 8.38642, 13.82750, 14.00083, 14.38333),
            new JGMTestData(1983, 5.65329, 4.91071, 9.32333, 10.61250, 11.76500),
            new JGMTestData(1984, 4.24734, 3.08629, 11.09750, 11.90750, 12.73833),
            new JGMTestData(1985, 3.87488, 2.53609, 9.45750, 10.38583, 11.11417),
            new JGMTestData(1986, 4.08913, 2.35438, 8.99083, 9.21333, 9.54417),
            new JGMTestData(1987, 4.27228, 4.60797, 8.16917, 9.41583, 9.92667),
            new JGMTestData(1988, 3.94387, 4.52556, 9.41583, 9.77167, 10.22750),
            new JGMTestData(1989, 4.87426, 4.72576, 12.01583, 10.20333, 9.92167),
            new JGMTestData(1990, 4.65466, 3.10150, 12.80500, 11.19250, 10.81167),
            new JGMTestData(1991, 5.46325, 2.86143, 8.83008, 9.16250, 9.80667),
            new JGMTestData(1992, 1.49461, 1.22810, 6.50875, 7.43167, 8.77167),
            new JGMTestData(1993, 1.82462, 1.04732, 4.92675, 6.45833, 7.87667),
            new JGMTestData(1994, 0.18511, 0.60929, 5.41675, 7.78667, 8.58000),
    };

    private static double[] data(Function<JGMTestData,Double> function) {
        double[] rS = new double [jgmTestData.length];

        for (int i = 0; i < jgmTestData.length; i++) {
            rS[i] = function.apply(jgmTestData[i]);
        }

        return rS;
    }

    private final static double[] residuals = {
            -11.5513833, -12.8690876, -12.5625016, -12.1318672, -7.6358242, -10.6873806, -8.3849936, -10.2063080, -10.7940123,
            -11.9678644, -8.1320178, -4.1110958, -4.9068512, -7.2117759, -5.7826065, -4.9517302, -7.1836791, -7.8458242,
            -9.6652355, -12.3912784, -12.9561574, -12.4333779, -14.1759612, -12.5133779, -10.4616253, -9.3428028, -8.0231496,
            -5.8271524, -6.9470018, -4.2899319, -5.3944641, -3.6223190, -1.8310045, 2.2769416, 4.5610357, 5.2267910,
            0.6094337, 0.5817750, 1.6624551, 2.8992830, -0.4952948, 0.9677267, -0.9429534, 0.2185573, -0.8595394,
            -3.1321684, -1.6144185, -1.1969561, 0.7468960, 2.3046003, -0.3282250, -0.0856873, -1.0233459, 2.1904013,
            3.7100088, 5.5604469, 6.5557641, 5.7787399, 7.0615651, 7.9507345, 8.2907345, 8.3233178, 10.0010220,
            9.9699038, 11.0726833, 10.0630759, 8.7173205, 6.2262023, 5.0421539, 3.6989362, 3.8600088, 1.6688905,
            -0.6875905, -3.9414883, -4.8842222, -3.7239803, -4.6595394, -5.0938297, -4.5633915, -4.9042679, -6.1519721,
            -6.3551442, -5.1300232, -2.5006120, 0.6754309, 0.7603100, 2.1275761, -0.9896900, 1.4660653, -1.5780743,
            -1.7531953, -3.1115340, -4.3396764, -7.5046467, -7.5410365, -8.4914746, -6.5586494, -6.3918672, -5.2769425,
            -3.6234828, -3.9866093, -2.8295394, -2.5708082, -1.2705663, -1.5529077, -0.9456417, 4.7914738, 5.8081055,
            5.9942534, 3.3335733, 2.7720626, 2.0769416, 5.2965948, 6.9910813, 7.3063985, 7.8400544, 10.6960517,
            4.9071836, 4.8612319, 6.9776674, 5.7038152, 5.4942077, 5.1414282, 5.3844040, 5.9105975, 5.1801594,
            3.4590411, 3.3105062, 4.7631352, 5.3226970, 7.5867910, 11.2220033, 10.6734684, 8.5663985, 7.2988905,
            4.7654309, 6.6454766, 7.5686943, 9.2253260, 0.4525601, 2.4302643, 5.4294793, 7.4163072, 6.8671379,
            9.4683018, 12.2507345, 13.7397075, 9.8619576, 8.4317157, 8.5183018, 10.2672748, 9.7421539, 6.2814282,
            4.9694793, 5.8635277, 9.9242534, 4.0668504, 1.5423182, 3.0228020, -0.6407033, -4.7933779, -3.4578188,
            -1.0465043, -3.5440123, -1.1871387, -3.4520634, -1.5697220, -6.0762487
    };


    @Test
    public void testOne() {
        UnitRootEvaluator evaluator;
        evaluator = new UnitRootEvaluator(data(JGMTestData::rS), 1, UnitRootEvaluator.Type.CONSTANT);

        System.out.printf("critical value (tau) = %f\n", UnitRootEvaluator.criticalValue(0.95, 0, 1, TestType.TAU, RegressionType.C));
        System.out.printf("critical value (z) = %f\n", UnitRootEvaluator.criticalValue(0.95, 0, 1, TestType.Z, RegressionType.C));

        System.out.printf("p-value (tau) = %f\n", evaluator.pValue(TestType.TAU));
        System.out.printf("p-value (z)   = %f\n", evaluator.pValue(TestType.Z));

        System.out.printf("p-value (tau) = %f\n", UnitRootEvaluator.pValue(-1.775, 0, 1, TestType.TAU, RegressionType.C));
        System.out.printf("p-value (z)   = %f\n", UnitRootEvaluator.pValue(-6.339, 0, 1, TestType.Z,   RegressionType.C));

        System.out.printf("p-value (tau) = %f\n", UnitRootEvaluator.pValue(-1.775, 37, 1, TestType.TAU, RegressionType.C));
        System.out.printf("p-value (z)   = %f\n", UnitRootEvaluator.pValue(-6.339, 37, 1, TestType.Z,   RegressionType.C));

        System.out.printf("p-value (tau) = %f\n", UnitRootEvaluator.pValue(-2.000, 0, 1, TestType.TAU, RegressionType.C));
        System.out.printf("p-value (z)   = %f\n", UnitRootEvaluator.pValue(-7.842, 0, 1, TestType.Z,   RegressionType.C));

        System.out.printf("p-value (tau) = %f\n", UnitRootEvaluator.pValue(-2.000, 37, 1, TestType.TAU, RegressionType.C));
        System.out.printf("p-value (z)   = %f\n", UnitRootEvaluator.pValue(-7.842, 37, 1, TestType.Z,   RegressionType.C));
    }
}
