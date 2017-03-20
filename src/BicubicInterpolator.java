import org.jetbrains.annotations.Contract;

import java.util.Arrays;

/**
 * Cubic interpolator for two-dimensional grid function.
 */
public class BicubicInterpolator {

    private final GridFunction2D gridFunction; // reference to grid function, which we interpolate
    private final double[][] coefficients; // spline coefficients

    @Contract("null -> fail")
    public BicubicInterpolator(GridFunction2D gridFunction){
        if (gridFunction == null) {
            throw new IllegalArgumentException("Grid function should not be NULL.");
        }
        this.gridFunction = gridFunction;
        coefficients = Interpolation.cubicBSplineTransform2D(gridFunction.values);
    }

    public double valueAt(double x, double y) throws OutOfGridException {
        final Localizator loc = localize(x, y);
        // Interpolation along second dimension
        double[] c = new double[4]; // here will be values of interpolations along second dimension
        for (int i = 0; i < 4; i++) {
            double[] cInter = Arrays.copyOfRange(coefficients[loc.iFrom + i], loc.jFrom, loc.jTo);
            c[i] = Interpolation.deBoor3(cInter, loc.py);
        }
        // Interpolation along first dimension
        return Interpolation.deBoor3(c, loc.px);
    }

    public double[] valueAndDerivativeAt(double x, double y) throws OutOfGridException {
        final Localizator loc = localize(x, y);
        // Interpolation along second dimension
        double[] c = new double[4]; // here will be values of interpolations along second dimension
        double[] dY = new double[4]; // here will be derivatives of interpolations along second dimension
        for (int i = 0; i < 4; i++) {
            double[] cInter = Arrays.copyOfRange(coefficients[loc.iFrom + i], loc.jFrom, loc.jTo);
            // TODO: Possible optimization: value and derivative in single deBoor algorithm
            c[i] = Interpolation.deBoor3(cInter, loc.py);
            double[] dInter = {(cInter[1] - cInter[0])/gridFunction.spacingY,
                               (cInter[2] - cInter[1])/gridFunction.spacingY,
                               (cInter[3] - cInter[2])/gridFunction.spacingY};
            dY[i] = Interpolation.deBoor2(dInter, loc.py);
        }
        // Interpolation along first dimension
        double[] result = new double[3]; // space for value and both partial derivatives
        double[] dX = {(c[1] - c[0])/gridFunction.spacingX,
                       (c[2] - c[1])/gridFunction.spacingX,
                       (c[3] - c[2])/gridFunction.spacingX};
        result[0] = Interpolation.deBoor3(c, loc.px);
        result[1] = Interpolation.deBoor2(dX, loc.px);
        result[2] = Interpolation.deBoor3(dY, loc.px);
        return result;
    }

    private class Localizator {

        final int iFrom;
        final double px;
        final int jFrom;
        final int jTo;
        final double py;

        public Localizator(int iFrom, double px, int jFrom, double py) {
            this.iFrom = iFrom;
            this.px = px;
            this.jFrom = jFrom;
            this.jTo = jFrom + 4;
            this.py = py;
        }

    }

    private Localizator localize(double x, double y) throws OutOfGridException {
        double r = x/gridFunction.spacingX;
        int iFrom = (int) r;
        double px = r - Math.floor(r);
        if (iFrom < 0 || iFrom >= gridFunction.sizeX) {
            throw new OutOfGridException();
        } else if (iFrom == gridFunction.sizeX - 1) {
            if (px == 0) {
                iFrom--;
                px = gridFunction.spacingX;
            } else {
                throw new OutOfGridException();
            }
        }
        r = y/gridFunction.spacingY;
        int jFrom = (int) r;
        double py = r - Math.floor(r);
        if (jFrom < 0 || jFrom >= gridFunction.sizeY) {
            throw new OutOfGridException();
        } else if (jFrom == gridFunction.sizeY - 1) {
            if (py == 0) {
                jFrom--;
                py = gridFunction.spacingY;
            } else {
                throw new OutOfGridException();
            }
        }
        return new Localizator(iFrom, px, jFrom, py);
    }
}
