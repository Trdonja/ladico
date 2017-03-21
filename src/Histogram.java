import java.util.PrimitiveIterator;

/**
 * One dimensional histogram, smoothed with cubic B-spline Parzen window. All histogram bins have an equal length.
 */
public class Histogram {

    final double[] bin;
    final double binWidth;
    final double[] frequency;

    /**
     * Constructs histogram by sorting given values from given interval into equal sized bins.
     * It is assumed that values lie on interval [lowerBound, upperBound].
     * @param lowerBound the lower bound of histogram domain (where the first bin is located).
     * @param upperBound the upper bound of histogram domain (where the last bin is located).
     * @param numberOfBins number of histogram bins.
     * @param values values to be sorted into histogram.
     */
    public Histogram(double lowerBound, double upperBound, int numberOfBins, PrimitiveIterator.OfDouble values) {
        if (lowerBound >= upperBound) {
            throw new IllegalArgumentException("Lower bound of histogram must be lower than upper bound.");
        }
        if (numberOfBins < 2) {
            throw new IllegalArgumentException("Number of bins must be greater than 1.");
        }
        binWidth = (upperBound - lowerBound)/(numberOfBins - 1);
        bin = new double[numberOfBins];
        int last = numberOfBins - 1;
        double binValue = lowerBound;
        for (int i = 0; i < last; i++) {
            bin[i] = binValue;
            binValue += binWidth;
        }
        bin[last] = upperBound; // explicitly, because of possible rounding errors in loop sumation
        frequency = new double[numberOfBins];
        // binning values into histogram (calculate frequencies):
        while (values.hasNext()) {
            double x = values.nextDouble();
            double r = (x - lowerBound)/binWidth;
            int iFrom = (int) r; // bins to be updated are iFrom .. iFrom + 3 (except for iFrom = 0 or last - 1)
            double p = r - Math.floor(r);
            double[] v = Interpolation.coxDeBoor3(p); // calculate incremental values for relavant bins
            if (iFrom == 0) { // shouldn't be less than 0, because values[:] >= lowerBound
                // skip v[0]
                for (int i = 1; i < 4; i++) {
                    frequency[iFrom - 1 + i] += v[i];
                }
            } else if (iFrom >= last - 1) { // should be ==last-1 or (==last and r==0), because values[:] <= upperBound
                // skip v[end]
                // FIXME if iFrom == last, then in foor loop its frequency[iFrom - 2 + i]
                for (int i = 0; i < 3; i++) {
                    frequency[iFrom - 1 + i] += v[i];
                }
            } else {
                for (int i = 0; i < 4; i++) {
                    frequency[iFrom - 1 + i] += v[i];
                }
            }
        }

    }

}
