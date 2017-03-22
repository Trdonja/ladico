import java.util.Arrays;
import java.util.Optional;

/**
 * Gradient of mutual information with respect to transformation parameter vector.
 */
public class MutualInformationGradient {

    /* TODO:
    * include histogram hR
    * rename class to MutualInformation
    * create methods value() and gradient() instead of get()
    * */
    // histogram arrays
    private final double[] hT; // histogram for target image, smoothed with cubic Parzen window
    private final double[][] h; // joint histogram for reference and target image
    /* Quasi gradient of joint probability with respect to transformation parameter vector,
       ∂h/∂μ. Needs to be multiplied with proper constant to be true probability density function. */
    private final double[][][] dp;
    // number of bins and width of bins used in histograms
    private final int nKappa; // number of bins in histogram for reference image
    private final int nIota; // number of bins in histogram for target image
    private final double epsR; // bin interval length in histogram for reference image
    private final double epsT; // bin interval length in histogram for target image
    // indicates which Parezn window is used for reference image intensities histogram
    private final ParzenWindow parzenWindow;
    // number of parameters of geometric transformation (dimension of transformation model)
    private final int parameters;
    // number of sample points used for evaluation of MI gradient
    private int cardinality;
    // values of MI gradient
    private boolean evaluated; // indicator, if MI gradient is evaluated or not
    private final double[] gradient; // gradient of MI gradient, if it is evaluated

    /**
     * Enumeration for discrimination of different Parzen windows used for reference image histogram.
     */
    public enum ParzenWindow {
        CUBIC_BSPLINE_WINDOW, // identifies cubic B-spline as Parzen window for ref. image
        LINEAR_BSPLINE_WINDOW, // identifies linear B-spline as Parzen window for ref. image
        NO_WINDOW // identifies absence of Parzen window for ref. image
        // For target image, cubic B-spline Parzen window is always used, because it has continuous derivative
    }

    /**
     * Initializes object for mutual information gradient calculation.
     * @param binsR number of bins used in histogram for reference image.
     * @param binsT number of bins used in histogram for target image.
     * @param maxfR maximum intensity gradient of reference image (minimum is assumed to be 0). I trust you!
     * @param maxfT maximum intensity gradient of target image (minimum is assumed to be 0). I trust you!
     * @param parameters number of geometric transformation parameters (degree of freedom of transformation model)
     * @param parzenWindow Parzen window used for reference image histogram
     */
    public MutualInformationGradient(int binsR, int binsT, double maxfR, double maxfT, int parameters,
                                     ParzenWindow parzenWindow) {
        if (parzenWindow == ParzenWindow.CUBIC_BSPLINE_WINDOW) {
            nKappa = binsR + 2; // we add two more bins, one at the beginning and one at the end
        } else { // parzenWindow is LINEAR_BSPLINE_WINDOW or NO_WINDOW
            nKappa = binsR;
        }
        this.parzenWindow = parzenWindow;
        nIota = binsT + 2; // we add two more bins because of use of cubic Parzen window for target image histogram
        epsR = maxfR/(binsR - 1);
        epsT = maxfT/(binsT - 1);

        hT = new double[nIota];
        h = new double[nKappa][nIota];
        dp = new double[parameters][nKappa][nIota];

        this.parameters = parameters;
        cardinality = 0;

        evaluated = false;
        gradient = new double[parameters];
    }

    /**
     * Include functional values of new sample point into internal components of gradient of mutual information.
     * This method serves as a builder for mutual information gradient. After functional values of all sample points
     * have been included, method evaluate() should be called to finalize the computation of gradient of mutual
     * information. After that, method get() gives the final gradient.
     * @param fR reference image intensity at sample point.
     * @param fT target image intensity at sample point.
     * @param gradfT gradient of image function (with respect to transformation parameter vector) at sample point.
     *               Its length must be equal to length of parameter vector.
     */
    public void include(double fR, double fT, double[] gradfT) {
        double r = fT/epsT;
        int jFrom = (int) r; // to know which bins should be updated
        if (jFrom == nIota - 3) { // if functional value equals to maximum, fix starting index:
            jFrom --;
        }
        double x = r - Math.floor(r);
        double[] contributionsFromT = Interpolation.coxDeBoor3WithDerivatives(x);
        for (int i = 4; i < 8; i++) { // to get actual derivatives, raw derivative values must be divided by bin width
            contributionsFromT[i] = contributionsFromT[i] / epsT;
        }

        r = fR/epsR;
        int kFrom = (int) r; // to know which bins should be updated
        x = r - Math.floor(r);
        if (parzenWindow == ParzenWindow.CUBIC_BSPLINE_WINDOW) {
            if (kFrom == nKappa - 3) { // if functional value equals to maximum, fix starting index:
                kFrom--;
            }
            double[] contributionsFromR = Interpolation.coxDeBoor3(x);
            // updating histogram h and ∂h/∂μ
            for (int k = 0; k < 4; k++) {
                int kNow = kFrom + k;
                // TODO: hR[kNow] += contributionsFromR[k];
                for (int j = 0; j < 4; j++) {
                    int jNow = jFrom + j;
                    h[kNow][jNow] += contributionsFromR[k] * contributionsFromT[j]; // update h
                    // update gradient of h
                    double prod = contributionsFromR[k] * contributionsFromT[4 + j];
                    for (int i = 0; i < parameters; i++) {
                        dp[i][kNow][jNow] += prod * gradfT[i];
                    }
                }
            }
        } else if (parzenWindow == ParzenWindow.LINEAR_BSPLINE_WINDOW) {
            if (kFrom == nKappa - 1) { // if functional value equals to maximum, fix starting index:
                kFrom --;
            }
            double emx = 1 - x; // contributions from R are now {1-x, x}
            // TODO: hR[kNow] += emx; hR[kNow + 1] += x;
            for (int j = 0; j < 4; j++) {
                int jNow = jFrom + j;
                h[kFrom][jNow] += emx * contributionsFromT[j];
                h[kFrom + 1][jNow] += x * contributionsFromT[j];
                double prod0 = emx * contributionsFromT[4 + j];
                double prod1 = x * contributionsFromT[4 + j];
                for (int i = 0; i < parameters; i++) {
                    dp[i][kFrom][jNow] += prod0 * gradfT[i];
                    dp[i][kFrom + 1][jNow] += prod1 * gradfT[i];
                }
            }
        } else if (parzenWindow == ParzenWindow.NO_WINDOW) { // the only remaining option
            if (kFrom == nKappa - 1) { // if functional value equals to maximum, fix starting index:
                kFrom --;
            }
            if (x >= 0.5) { // decides which bin to increment by 1, left to fR (x < 0.5) or right to fR (x >= 0.5)
                kFrom ++;
            }
            // TODO: hR[kFrom] += 1;
            for (int j = 0; j < 4; j++) {
                int jNow = jFrom + j;
                h[kFrom][jNow] += contributionsFromT[j];
                double prod = contributionsFromT[4 + j];
                for (int i = 0; i < parameters; i++) {
                    dp[i][kFrom][jNow] += prod * gradfT[i];
                }
            }
        }

        cardinality ++; // one more sample was included
        if (evaluated) {
            evaluated = false;
            // these two need to be calculated again in method evaluate() and are set to initial value here:
            Arrays.fill(hT, 0);
            Arrays.fill(gradient, 0);
        }
    }

    /**
     * Finalizes the computation of gradient of mutual information, after functional values of all sample points have
     * been included using include() method.
     */
    public void evaluate() {
        if (!evaluated) {
            // calculate marginal histogram hT
            for (int j = 0; j < nIota; j++) {
                for (int k = 0; k < nKappa; k++) {
                    hT[j] += h[k][j];
                }
            }
            // calculate gradient
            for (int k = 0; k < nKappa; k++) {
                for (int j = 0; j < nIota; j++) {
                    double factor = Math.log(h[k][j] / hT[j]); // FIXME: hT[j] can be zero!
                    for (int i = 0; i < parameters; i++) {
                        gradient[i] += dp[i][k][j] * factor;
                    }
                }
            }
            double c = 1 / (cardinality * epsT);
            for (int i = 0; i < parameters; i++) {
                gradient[i] = c * gradient[i];
            }

            evaluated = true;
        }
    }

    /**
     * Returns value of gradient of mutual information, if already calculated.
     * @return value of gradient of mutual information, if already calculated, otherwise empty Optional.
     */
    public Optional<double[]> get() {
        if (evaluated) {
            return Optional.of(Arrays.copyOf(gradient, gradient.length));
        } else {
            return Optional.empty();
        }
    }

}
