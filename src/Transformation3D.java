import java.util.Arrays;

/**
 * Parameterized transformation from R^3 to R^3.
 */
public abstract class Transformation3D {

    private final double[] parameters;

    public Transformation3D(double[] parameters) {
        this.parameters = Arrays.copyOf(parameters, parameters.length);
    }

    /**
     * Applies this transformation to given point and returns the result.
     * @param point 3-dimensional point to be mapped.
     * @return result of mapping.
     */
    public abstract Point3D map(Point3D point);

    /**
     * Applies this transformation to given point and at the same time calculates gradient of this map with
     * respect to transformation parameter vector, calculated at given parameters and given point.
     * @param point 3-dimensional point, where map and gradient should be calculated.
     * @return resulting point and gradient.
     */
    public abstract Pair<Point3D, double[][]> pointAndGradientAt(Point3D point);

    public int dimensionOfParameterSpace() {
        return parameters.length;
    }
}
