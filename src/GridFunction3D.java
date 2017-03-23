/**
 * Real valued function on three-dimensional grid.
 */
public class GridFunction3D {

    final int sizeX;
    final int sizeY;
    final int sizeZ;
    final double spacingX;
    final double spacingY;
    final double spacingZ;
    final double[][][] values;

    /**
     * Order, in which values should be inserted into grid function
     */
    public enum Ordering {
        IJK, // Values are inserted first in Z-direction, then in Y-direction and then in X-direction
        KJI  // Values are inserted first in X-direction, then in Y-direction and then in Z-direction
        // Additional values may come in the future.
    }

    public GridFunction3D(int sizeX, int sizeY, int sizeZ,
                          double spacingX, double spacingY, double spacingZ,
                          double[] values, Ordering ordering) {
        if (sizeX <= 0 || sizeY <= 0 || sizeZ <= 0) {
            throw new IllegalArgumentException("Grid size must be positive.");
        }
        if (spacingX <= 0 || spacingY <= 0 || spacingZ <= 0) {
            throw new IllegalArgumentException("Grid spacing must be positive.");
        }
        if (values == null) {
            throw new IllegalArgumentException("Values are not provided.");
        }
        if (values.length != sizeX*sizeY*sizeZ) {
            throw new IllegalArgumentException("The number of values does not match grid size.");
        }
        this.sizeX = sizeX;
        this.sizeY = sizeY;
        this.sizeZ = sizeZ;
        this.spacingX = spacingX;
        this.spacingY = spacingY;
        this.spacingZ = spacingZ;
        this.values = new double[sizeX][sizeY][sizeZ];
        int index = 0;
        if (ordering == Ordering.IJK) {
            for (int i = 0; i < sizeX; i++) {
                for (int j = 0; j < sizeY; j++) {
                    for (int k = 0; k < sizeZ; k++) {
                        this.values[i][j][k] = values[index];
                        index++;
                    }
                }
            }
        } else if (ordering == Ordering.KJI) {
            for (int k = 0; k < sizeZ; k++) {
                for (int j = 0; j < sizeY; j++) {
                    for (int i = 0; i < sizeX; i++) {
                        this.values[i][j][k] = values[index];
                        index++;
                    }
                }
            }
        }
    }

}
