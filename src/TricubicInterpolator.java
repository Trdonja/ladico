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
        double[][] c2 = new double[4][4]; // here will be values of interpolations along third dimension
        final int kTo = loc.kFrom + 4;
        for (int i = 0; i < 4; i++) {
            final int iLocal = loc.iFrom + i;
            for (int j = 0; j < 4; j++) {
                double[] cInter = Arrays.copyOfRange(coefficients[iLocal][loc.jFrom + j], loc.kFrom, kTo);
                c2[i][j] = Interpolation.deBoor3(cInter, loc.pz);
            }
        }
        double[] c1 = new double[4];
        for (int i = 0; i < 4; i++) {
            c1[i] = Interpolation.deBoor3(c2[i], loc.py);
        }
        return Interpolation.deBoor3(c1, loc.px);
    }

    public double[] valueAndDerivativeAt(double x, double y, double z) throws OutOfGridException {
        // TODO: Implement this method
        return new double[4]; // value and 3 gradient components
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
