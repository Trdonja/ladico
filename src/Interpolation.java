import org.jetbrains.annotations.Contract;

import java.util.Arrays;

/**
 * Functions for continuous interpolation of grid data.
 */
public class Interpolation {

    /**
     * Calculates coefficients of cubic B-spline, which interpolate values on 1D grid.
     * @param values grid values to be interpolated with cubic B-spline
     * @return coefficients of cubic B-spline, which interpolates given values
     */
    @Contract(pure = true)
    static double[] cubicBSplineTransform(double[] values) {
        final int n = values.length;
        double[] a = new double[n];
        double[] b = Arrays.stream(values).map(x -> 6*x).toArray();
        double[] c = new double[n+2];
        // Forward filtering
        a[0] = 0;
        b[0] = b[0]/2;
        for (int i = 1; i < n-1; i++) {
            double factor = -1/a[i-1];
            a[i] = 4 + factor;
            b[i] += factor*b[i-1];
        }
        double factor = -2/a[n-2];
        a[n-1] = 4 + factor;
        b[n-1] += factor*b[n-2];
        // Backwards filtering
        c[n] = b[n-1]/a[n-1];
        for (int i = n-2; i > 0; i--) {
            c[i+1] = (b[i] - c[i+2])/a[i];
        }
        c[0] = c[2];
        c[n+1] = c[n-1];

        return c;
    }

    /**
     * Calculates coefficients of cubic B-spline, which interpolate values on 2D grid.
     * @param values grid values to be interpolated with cubic B-spline
     * @return coefficients of cubic B-spline, which interpolates given values
     */
    @Contract(pure = true)
    static double[][] cubicBSplineTransform2D(double[][] values) {
        final int m = values.length;
        final int n = values[0].length;
        // Apply transform row-wise
        double[][] intermediate = new double[m][n+2];
        for (int i = 0; i < m; i++) {
            intermediate[i] = cubicBSplineTransform(values[i]);
        }
        // Apply transform column-wise
        double[][] result = new double[m+2][n+2];
        for (int j = 0; j < n+2; j++) {
            double[] column = new double[m];
            int finalJ = j;
            Arrays.setAll(column, i -> intermediate[i][finalJ]);
            double[] c = cubicBSplineTransform(column);
            for(int i = 0; i < m+2; i++) {
                result[i][j] = c[i];
            }
        }
        return result;
    }

    /**
     * DeBoor's algorithm for calculation of value of cubic B-spline at some point from domain.
     * @param c coefficients of B-spline. The length of this array must be exactly 4.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing)
     * @return value of cubic B-spline at x
     */
    @Contract(pure = true)
    static double deBoor3(double[] c, double p) {
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
     * DeBoor's algorithm for calculation of value of quadric B-spline at some point from domain.
     * @param c coefficients of B-spline. The length of this array must be exactly 3.
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing)
     * @return value of quadric B-spline at x
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
     * @param p ratio of domain point x, p = x - floor(x/gridSpacing)
     * @return value of linear B-spline at x
     */
    @Contract(pure = true)
    static double deBoor1(double[] c, double p) {
        // coefficient c^[1](x) = s(x)
        return c[0] + p*(c[1] - c[0]);
    }

    //private double cubicInterpolate(double x, coefficients)
}
