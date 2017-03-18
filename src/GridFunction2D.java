/**
 * Real valued function on two-dimensional grid.
 */
public class GridFunction2D {

    final int sizeX;
    final int sizeY;
    final double spacingX;
    final double spacingY;
    final double[][] values;

    /**
     * Constructs grid function from given data.
     * @param sizeX size of the first dimension.
     * @param sizeY size of the second dimension.
     * @param spacingX grid spacing along first dimension.
     * @param spacingY grid spacing along second dimension.
     * @param values values of function in grid points.
     * @param columnOrder if true, values are rewritten in 2D array in column order, otherwise in row order.
     */
    public GridFunction2D(int sizeX, int sizeY, int spacingX, int spacingY, double[] values, boolean columnOrder) {
        if (sizeX <= 0 || sizeY <= 0) {
            throw new IllegalArgumentException("Grid size must be positive.");
        }
        if (spacingX <= 0 || spacingY <= 0) {
            throw new IllegalArgumentException("Grid spacing must be positive.");
        }
        if (values == null) {
            throw new IllegalArgumentException("Values are not provided.");
        }
        if (values.length != sizeX*sizeY) {
            throw new IllegalArgumentException("The number of values does not match grid size.");
        }
        this.sizeX = sizeX;
        this.sizeY = sizeY;
        this.spacingX = spacingX;
        this.spacingY = spacingY;
        this.values = new double[sizeX][sizeY];
        int index = 0;
        if (!columnOrder) {
            for (int i = 0; i < sizeX; i++) {
                for (int j = 0; j < sizeY; j++) {
                    this.values[i][j] = values[index];
                    index++;
                }
            }
        } else {
            for (int j = 0; j < sizeY; j++) {
                for (int i = 0; i < sizeX; i++) {
                    this.values[i][j] = values[index];
                    index++;
                }
            }
        }
    }

}
