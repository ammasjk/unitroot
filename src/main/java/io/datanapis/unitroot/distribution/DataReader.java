package io.datanapis.unitroot.distribution;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Author: Jayakumar Muthukumarasamy
 *
 * Read and return response surface coefficients produced by James G MacKinnon
 * James G. MacKinnon, "Numerical distribution functions for unit root and cointegration tests,"
 *                     Journal of Applied Econometrics, 11, 1996, 601-618.
 */
class DataReader {
    Pattern DATA_START_PATTERN = Pattern.compile("^([a-z]+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)");

    public MackinnonData read(int niv, TestType tt, RegressionType rt) throws IOException {
        String desiredTag = getTag(niv, tt, rt);
        String name = name(niv);
        InputStream inputStream = getFileStream(name);
        if (inputStream == null)
            return null;

        try (InputStreamReader reader = new InputStreamReader(inputStream);
             BufferedReader bufferedReader = new BufferedReader(reader)) {

            String line;
            while ((line = bufferedReader.readLine()) != null) {
                /* Skip copyright line */
                if (line.startsWith("Copyright"))
                    continue;

                Matcher matcher = DATA_START_PATTERN.matcher(line);
                if (matcher.find()) {
                    String tag = matcher.group(1);
                    if (!desiredTag.equals(tag)) {
                        for (int i = 0; i < MackinnonData.NROWS; i++) {
                            bufferedReader.readLine();
                        }
                    } else {
                        int nz = Integer.parseInt(matcher.group(2));
                        int nreg = Integer.parseInt(matcher.group(3));
                        int model = Integer.parseInt(matcher.group(4));
                        int minSize = Integer.parseInt(matcher.group(5));

                        int count = count(model);
                        double[] data = new double [MackinnonData.NROWS * count];

                        for (int i = 0; i < MackinnonData.NROWS; i++) {
                            line = bufferedReader.readLine();
                            if (line == null) {
                                throw new IOException("EOF when EOF not expected!");
                            }

                            line = line.trim();
                            line = line.replace('D', 'E');

                            String[] parts = line.split("\\s+");
                            assert parts.length == count;
                            for (int j = 0; j < count; j++) {
                                data[i * count + j] = Double.parseDouble(parts[j]);
                            }
                        }

                        inputStream.close();
                        return new MackinnonData(tag, nz, nreg, model, minSize, data);
                    }
                } else {
                    inputStream.close();
                    throw new RuntimeException("Invalid line [" + line + "] in data file [" + name + "]");
                }
            }
        }

        inputStream.close();
        return null;
    }

    private static String choose(TestType tt, RegressionType rt,
                                 String nc, String c, String ct, String ctt,
                                 String anc, String ac, String act, String actt) {
        if (tt == TestType.TAU && rt == RegressionType.NC)
            return nc;
        else if (tt == TestType.TAU && rt == RegressionType.C)
            return c;
        else if (tt == TestType.TAU && rt == RegressionType.CT)
            return ct;
        else if (tt == TestType.TAU && rt == RegressionType.CTT)
            return ctt;
        else if (tt == TestType.Z && rt == RegressionType.NC)
            return anc;
        else if (tt == TestType.Z && rt == RegressionType.C)
            return ac;
        else if (tt == TestType.Z && rt == RegressionType.CT)
            return act;
        else
            return actt;
    }

    public static String getTag(int niv, TestType tt, RegressionType rt) {
        switch (niv) {
            case 1:
                return choose(tt, rt, "dfnc", "dfc", "dfct", "dfctt",
                        "dfanc", "dfac", "dfact", "dfactt");
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
            case 10:
            case 11:
            case 12:
                return choose(tt, rt, "conc", "coc", "coct", "coctt",
                        "coanc", "coac", "coact", "coactt");
            default:
                throw new RuntimeException("Invalid combination");
        }
    }

    private static String name(int niv) {
        return String.format("/urc-%d.tab", niv);
    }

    private static InputStream getDataStream(String name) {
        return DataReader.class.getResourceAsStream(name);
    }

    private static InputStream getFileStream(String name) throws IOException {
        return new FileInputStream("src/main/resources" + name);
    }

    private static String devName(int niv) {
        return String.format("src/main/resources/urc-%d.tab", niv);
    }

    private static int count(int model) {
        if (model == 2 || model == 4)
            return 4;
        else
            return 5;
    }
}
