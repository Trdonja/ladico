import org.jetbrains.annotations.Contract;

import java.util.Arrays;

/**
 * Functions for continuous interpolation of grid data.
 */
public class Interpolation {

    /**
     * Calculates coefficients of cubic B-spline, which interpolates values on 1D grid.
     * @param values grid values to be interpolated with cubic B-spline
     * @return coefficients of cubic B-spline, which interpolates given values
     */
    @Contract(pure = true)
    static double[] cubicBSplineTransform(double[] values) {
        final int n = values.length;
        double[] a = new double[n];
        double[] b = Arrays.stream(values).map(x -> 6 * x).toArray();
        double[] c = new double[n + 2];
        // Forward filtering
        a[0] = 0;
        b[0] = b[0] / 2;
        for (int i = 1; i < n - 1; i++) {
            double factor = -1 / a[i - 1];
            a[i] = 4 + factor;
            b[i] += factor * b[i - 1];
        }
        double factor = -2 / a[n - 2];
        a[n - 1] = 4 + factor;
        b[n - 1] += factor*b[n - 2];
        // Backwards filtering
        c[n] = b[n - 1] / a[n - 1];
        for (int i = n - 2; i > 0; i--) {
            c[i + 1] = (b[i] - c[i + 2]) / a[i];
        }
        c[0] = c[2];
        c[n + 1] = c[n - 1];

        return c;
    }

    /**
     * Calculates coefficients of cubic B-spline, which interpolates values on 2D grid.
     * @param values grid values to be interpolated with cubic B-spline
     * @return coefficients of cubic B-spline, which interpolates given values
     */
    @Contract(pure = true)
    static double[][] cubicBSplineTransform2D(double[][] values) {
        final int m = values.length;
        final int n = values[0].length;
        // Apply 1D transform row-wise
        double[][] intermediate = new double[m][n + 2];
        for (int i = 0; i < m; i++) {
            intermediate[i] = cubicBSplineTransform(values[i]);
        }
        // Apply 1D transform column-wise
        double[][] result = new double[m + 2][n + 2];
        for (int j = 0; j < n + 2; j++) {
            double[] column = new double[m];
            final int finalJ = j;
            Arrays.setAll(column, i -> intermediate[i][finalJ]);
            double[] c = cubicBSplineTransform(column);
            for(int i = 0; i < m + 2; i++) {
                result[i][j] = c[i];
            }
        }
        return result;
    }

    /**
     * Calculates coefficients of cubic B-splines, which interpolates values on 3D grid.
     * @param values grid values to be interpolated with cubic B-spline
     * @return coefficients of cubic B-spline, which interpolates given values
     */
    @Contract(pure = true)
    static double[][][] cubicBSplineTransform3D(double[][][] values) {
        final int m = values.length;
        final int n = values[0].length;
        final int p = values[0][0].length;
        // Apply 2D transform by X-slices
        double[][][] intermediate = new double[m][n + 2][p + 2];
        for (int i = 0; i < m; i++) {
            intermediate[i] = cubicBSplineTransform2D(values[i]);
        }
        // Apply 1D transform in X-direction
        double[][][] result = new double[m + 2][n + 2][p + 2];
        for (int j = 0; j < n + 2; j++) {
            for (int k = 0; k < p + 2; k++) {
                double[] column = new double[m];
                final int finalJ = j;
                final int finalK = k;
                Arrays.setAll(column, i -> intermediate[i][finalJ][finalK]);
                double[] c = cubicBSplineTransform(column);
                for (int i = 0; i < m + 2; i++) {
                    result[i][j][k] = c[i];
                }
            }
        }
        return result;
    }

    /**
     * DeBoor's algorithm for calculation of value of cubic B-spline at some point from domain.
     * @param c coefficients of B-spline. The length of this array must be exactly 4.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing).
     * @return value of cubic B-spline at x
     */
    @Contract(pure = true)
    static double deBoor3(final double[] c, final double p) {
        // coefficients c^[1](x)
        c[0] = (2*c[1] + c[0] + p*(c[1] - c[0]))/3;
        c[1] = (c[2] + 2*c[1] + p*(c[2] - c[1]))/3;
        c[2] = c[2] + p*(c[3] - c[2])/3;
        // coefficients c^[2](x)
        c[0] = (c[1] + c[0] + p*(c[1] - c[0]))/2;
        c[1] = c[1] + p*(c[2] - c[1])/2;
        // coefficient c^[3](x) = s(x)
        return c[0] + p*(c[1] - c[0]);
    }

    /**
     * DeBoor's algorithm for calculation of value of cubic B-spline and its derivative at some point from domain.
     * @param c coefficients of B-spline. The length of this array must be exactly 4.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing).
     * @param spacing domain grid spacing (used in derivative calculation).
     * @return value and derivative of cubic B-spline at x.
     */
    @Contract(pure = true)
    static double[] deBoor3withDerivative(final double[] c, final double p, final double spacing) {
        double[] d = {c[1] - c[0], c[2] - c[1], c[3] - c[2]};
        // coefficients c^[1](x)
        c[0] = (2*c[1] + c[0] + p*d[0])/3;
        c[1] = (c[2] + 2*c[1] + p*d[1])/3;
        c[2] = c[2] + p*d[2]/3;
        // coefficients c^[2](x)
        c[0] = (c[1] + c[0] + p*(c[1] - c[0]))/2;
        c[1] = c[1] + p*(c[2] - c[1])/2;
        // coefficient c^[3](x) = s(x)
        c[0] = c[0] + p*(c[1] - c[0]);

        // Calculating derivative s'(x):
        // coefficients d^[1](x)
        d[0] = (d[1] + d[0] + p*(d[1] - d[0]))/2/spacing;
        d[1] = (d[1] + p*(d[2] - d[1])/2)/spacing;
        // coefficient d^[2](x) = s'(x)
        d[0] = d[0] + p*(d[1] - d[0]);

        return new double[]{c[0], d[0]};
    }

    /**
     * DeBoor's algorithm for calculation of value of quadric B-spline at some point from domain.
     * @param c coefficients of B-spline. The length of this array must be exactly 3.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing).
     * @return value of quadric B-spline at x.
     */
    @Contract(pure = true)
    static double deBoor2(double[] c, double p) {
        // coefficients c^[1](x)
        c[0] = (c[1] + c[0] + p*(c[1] - c[0]))/2;
        c[1] = c[1] + p*(c[2] - c[1])/2;
        // coefficient c^[2](x) = s(x)
        return c[0] + p*(c[1] - c[0]);
    }

    /**
     * Linear interpolation at some point from domain.
     * @param c neighbouring values of x. The length of this array must be exactly 2.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing).
     * @return value of linear B-spline at x.
     */
    @Contract(pure = true)
    static double deBoor1(double[] c, double p) {
        // coefficient c^[1](x) = s(x)
        return c[0] + p*(c[1] - c[0]);
    }

    /**
     * Calculates values of all basis cubic spline functions at given p.
     * @param p number from interval [0, 1].
     * @return values of all four basis cubic spline functions at p (length is 4).
     */
    @Contract(pure = true)
    static double[] coxDeBoor3(double p) {
        double emp = 1 - p;
        double[] b = new double[4];
        // first column
        b[2] = emp;
        b[3] = p;
        // second column
        b[1] = emp*b[2]/2;
        b[2] = b[3] + (b[2] + p*(b[2] - b[3]))/2;
        b[3] = p*b[3]/2;
        // third column
        b[0] = emp*b[1]/3;
        b[1] = (2*(b[1] + b[2]) + p*(b[1] - b[2]))/3;
        b[2] = b[3] + (b[2] + p*(b[2] - b[3]))/3;
        b[3] = p*b[3]/3;
        return b;
    }

    @Contract(pure = true)
    static double[] coxDeBoor3WithDerivatives(double p) {
        double emp = 1 - p;
        double[] b = new double[8]; // first 4 elements are values, last 4 elements are derivatives
        // first column
        b[2] = emp;
        b[3] = p;
        // second column
        b[1] = emp*b[2]/2;
        b[2] = b[3] + (b[2] + p*(b[2] - b[3]))/2;
        b[3] = p*b[3]/2;
        // derivatives (need to be scaled after this function ends, if interval width is not 1)
        b[4] = b[0] - b[1];
        b[5] = b[1] - b[2];
        b[6] = b[2] - b[3];
        b[7] = b[3];
        // third column
        b[0] = emp*b[1]/3;
        b[1] = (2*(b[1] + b[2]) + p*(b[1] - b[2]))/3;
        b[2] = b[3] + (b[2] + p*(b[2] - b[3]))/3;
        b[3] = p*b[3]/3;
        return b;
    }
}
