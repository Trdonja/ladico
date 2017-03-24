import java.util.Arrays;

/**
 * Parameterized transformation from R^2 to R^2.
 */
public abstract class Transformation2D {

    private final double[] parameters;

    public Transformation2D(double[] parameters) {
        this.parameters = Arrays.copyOf(parameters, parameters.length);
    }

    /**
     * Applies this transformation to given point and returns the result.
     * @param point 2-dimensional point to be mapped.
     * @return result of mapping.
     */
    public abstract Point2D map(Point2D point);

    /**
     * Applies this transformation to given point and at the same time calculates gradient of this map with
     * respect to transformation parameter vector, calculated at given parameters and given point.
     * @param point 3-dimensional point, where map and gradient should be calculated.
     * @return resulting point and gradient.
     */
    public abstract Pair<Point2D, double[][]> pointAndGradientAt(Point2D point);

    public int dimensionOfParameterSpace() {
        return parameters.length;
    }
}
