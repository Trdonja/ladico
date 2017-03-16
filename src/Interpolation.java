import java.util.Arrays;

/**
 * Created by domen on 16.3.2017.
 */
public class Interpolation {

    public static double[] cubicBSplineTransform(double[] values){
        final int n = values.length;
        double[] a = new double[n];
        double[] b = Arrays.stream(values).map(x -> 6*x).toArray();
        double[] c = new double[n+2];
        // Forward filtering
        a[0] = 0;
        b[0] = b[0]/2;
        for(int i = 1; i < n-1; i++){
            double factor = -1/a[i-1];
            a[i] = 4 + factor;
            b[i] += factor*b[i-1];
        }
        double factor = -2/a[n-2];
        a[n-1] = 4 + factor;
        b[n-1] += factor*b[n-2];
        // Backwards filtering
        c[n] = b[n-1]/a[n-1];
        for(int i = n-2; i > 0; i--){
            c[i+1] = (b[i] - c[i+2])/a[i];
        }
        c[0] = c[2];
        c[n+1] = c[n-1];

        return c;
    }

    public static double[][] cubicBSplineTransform2D(double[][] values){
        final int m = values.length;
        final int n = values[0].length;
        // Apply transform row-wise
        double[][] intermediate = new double[m][n+2];
        for(int i = 0; i < m; i++){
            intermediate[i] = cubicBSplineTransform(values[i]);
        }
        // Apply transform column-wise
        double[][] result = new double[m+2][n+2];
        for(int j = 0; j < n+2; j++){
            double[] column = new double[m];
            int finalJ = j;
            Arrays.setAll(column, i -> intermediate[i][finalJ]);
            double[] c = cubicBSplineTransform(column);
            for(int i = 0; i < m+2; i++){
                result[i][j] = c[i];
            }
        }
        return result;
    }
}
