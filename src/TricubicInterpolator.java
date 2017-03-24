import org.jetbrains.annotations.Contract;

import java.util.Arrays;

/**
 * Cubic interpolator for two-dimensional grid function.
 */
public class TricubicInterpolator {

    private final GridFunction3D gridFunction; // reference to grid function, which we interpolate
    private final double[][][] coefficients; // spline coefficients

    @Contract("null -> fail")
    public TricubicInterpolator(GridFunction3D gridFunction) {
        if (gridFunction == null) {
            throw new IllegalArgumentException("Grid function should not be NULL.");
        }
        this.gridFunction = gridFunction;
        coefficients = Interpolation.cubicBSplineTransform3D(gridFunction.values);
    }

    public double valueAt(double x, double y, double z) throws OutOfGridException {
        final Localizator loc = localize(x, y, z);
        final double[][] c2 = new double[4][4]; // here will be values of interpolations along third dimension
        final int kTo = loc.kFrom + 4;
        for (int i = 0; i < 4; i++) {
            final int iLocal = loc.iFrom + i;
            for (int j = 0; j < 4; j++) {
                double[] c3 = Arrays.copyOfRange(coefficients[iLocal][loc.jFrom + j], loc.kFrom, kTo);
                c2[i][j] = Interpolation.deBoor3(c3, loc.pz);
            }
        }
        final double[] c1 = new double[4];
        for (int i = 0; i < 4; i++) {
            c1[i] = Interpolation.deBoor3(c2[i], loc.py);
        }
        return Interpolation.deBoor3(c1, loc.px);
    }

    public double[] valueAndDerivativeAt(double x, double y, double z) throws OutOfGridException {
        final Localizator loc = localize(x, y, z);
        final double c2[][] = new double[4][4]; // values of interpolation in Z-direction
        final double dz2[][] = new double[4][4]; // derivatives of interpolation in Z-direction
        final double c1[] = new double[4]; // values of interpolation in Z and Y-direction
        final double dy1[] = new double[4];
        final double dz1[] = new double[4];
        final int kTo = loc.kFrom + 4;
        for (int i = 0; i < 4; i++) {
            final int iLocal = loc.iFrom + i;
            for (int j = 0; j < 4; j++) {
                final double[] c3 = Arrays.copyOfRange(coefficients[iLocal][loc.jFrom + j], loc.kFrom, kTo);
                double[] valueDerivativePair = Interpolation.deBoor3withDerivative(c3, loc.pz, gridFunction.spacingZ);
                c2[i][j] = valueDerivativePair[0];
                dz2[i][j] = valueDerivativePair[1];
            }
            double[] valueDerivativePair = Interpolation.deBoor3withDerivative(c2[i], loc.py, gridFunction.spacingY);
            c1[i] = valueDerivativePair[0];
            dy1[i] = valueDerivativePair[1];
            dz1[i] = Interpolation.deBoor3(dz2[i], loc.py);
        }
        final double[] result = new double[4]; // storage for result; first cell for value, last 3 for derivatives
        double[] valueDerivativePair = Interpolation.deBoor3withDerivative(c1, loc.px, gridFunction.spacingX);
        result[0] = valueDerivativePair[0]; // value of spline at (x,y,z)
        result[1] = valueDerivativePair[1]; // partial derivative of spline with regard to x at (x,y,z)
        result[2] = Interpolation.deBoor3(dy1, loc.px); // partial derivative of spline with regard to y at (x,y,z)
        result[3] = Interpolation.deBoor3(dz1, loc.px); // partial derivative of spline with regard to z at (x,y,z)
        return result; // value and 3 gradient components
    }


    private class Localizator {

        final int iFrom;
        final double px;
        final int jFrom;
        final double py;
        final int kFrom;
        final double pz;

        public Localizator(int iFrom, double px, int jFrom, double py, int kFrom, double pz) {
            this.iFrom = iFrom;
            this.px = px;
            this.jFrom = jFrom;
            this.py = py;
            this.kFrom = kFrom;
            this.pz = pz;
        }
    }

    private Localizator localize(double x, double y, double z) throws OutOfGridException {
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
        r = z/gridFunction.spacingZ;
        int kFrom = (int) r;
        double pz = r - Math.floor(r);
        if (kFrom < 0 || kFrom >= gridFunction.sizeZ) {
            throw new OutOfGridException();
        } else if (kFrom == gridFunction.sizeZ - 1) {
            if (pz == 0) {
                kFrom--;
                pz = gridFunction.spacingZ;
            } else {
                throw new OutOfGridException();
            }
        }
        return new Localizator(iFrom, px, jFrom, py, kFrom, pz);
    }

}
