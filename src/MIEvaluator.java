import java.util.Random;

/**
 * Class for evaluation of mutual information.
 */
public class MIEvaluator {

    private final GridFunction3D reference;
    private final TricubicInterpolator target;
    private final Transformation3D transformation;

    public MIEvaluator(GridFunction3D reference, GridFunction3D target, Transformation3D transformation) {
        this.reference = reference;
        this.target = new TricubicInterpolator(target);
        this.transformation = transformation;
    }

    private MutualInformation createMutualInformation(int numberOfSamples) {
        Random random = new Random(System.currentTimeMillis());
        MutualInformation mi = new MutualInformation(32, 32, 1, 1,
                transformation.dimensionOfParameterSpace(), MutualInformation.ParzenWindow.LINEAR_BSPLINE_WINDOW);
        for (int i = 0; i < numberOfSamples; i++) {
            // calculate random coordinates on reference grid
            final int x = random.nextInt(reference.sizeX);
            final int y = random.nextInt(reference.sizeY);
            final int z = random.nextInt(reference.sizeZ);
            final Point3D u = new Point3D(x * reference.spacingX, y * reference.spacingY, z * reference.spacingZ);
            final Pair<Point3D, double[][]> pointGradientPair = transformation.pointAndGradientAt(u);
            final double[][] gradG = pointGradientPair.right;
            try {
                final double[] valueDerivativePair = target.valueAndDerivativeAt(pointGradientPair.left.x,
                        pointGradientPair.left.y, pointGradientPair.left.z);
                final double[] gradFT = new double[transformation.dimensionOfParameterSpace()];
                for (int j = 0; j < gradFT.length; j++) {
                    gradFT[j] = valueDerivativePair[1] * gradG[j][0] + valueDerivativePair[2] * gradG[j][1] +
                            valueDerivativePair[3] * gradG[j][2];
                }
                mi.include(reference.values[x][y][z], valueDerivativePair[0], gradFT);
            } catch (OutOfGridException e) {
                // TODO: Increase counter of points, which are out of grid.
            }
        }
        return mi;
    }
}
