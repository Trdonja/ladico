import java.util.Arrays;
import java.util.Optional;

/**
 * Gradient of mutual information with respect to transformation parameter vector.
 */
public class MutualInformationGradient {

    // probability arrays
    private double[] hT; // histogram for target image, smoothed with cubic Parzen window
    private double[][] h; // joint histogram for reference and target image
    /* Quasi gradient of joint probability with respect to transformation parameter vector,
       ∂h/∂μ. Needs to be multiplied with proper constant to be true probability density function. */
    private double[][][] dp;
    // number of bins and width of bins used in histograms
    private int nKappa; // number of bins in histogram for reference image
    private int nIota; // number of bins in histogram for target image
    private double epsR; // bin interval length in histogram for reference image
    private double epsT; // bin interval length in histogram for target image
    // number of parameters of geometric transformation (dimension of transformation model)
    private int parameters;
    // number of sample points used for evaluation of MI gradient
    private int cardinality;
    // values of MI gradient
    private boolean evaluated; // indicator, if MI gradient is evaluated or not
    private double[] gradient; // gradient of MI gradient, if it is evaluated


    /**
     * Initializes object for mutual information gradient calculation.
     * @param binsR number of bins used in histogram for reference image.
     * @param binsT number of bins used in histogram for target image.
     * @param maxfR maximum intensity gradient of reference image (minimum is assumed to be 0). I trust you!
     * @param maxfT maximum intensity gradient of target image (minimum is assumed to be 0). I trust you!
     * @param parameters number of geometric transformation parameters (degree of freedom of transformation model)
     */
    public MutualInformationGradient(int binsR, int binsT, double maxfR, double maxfT, int parameters) {
        nKappa = binsR + 2; // we add two more bins, one at the beginning and one at the end
        nIota = binsT + 2; // same here (because of cubic Parzen window)
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
        double r = fR/epsR;
        int kFrom = (int) r; // to know which bins should be updated
        double x = r - Math.floor(r);
        double[] kappa = Interpolation.coxDeBoor3(x);

        r = fT/epsT;
        int jFrom = (int) r; // to know which bins should be updated
        x = r - Math.floor(r);
        double[] iota = Interpolation.coxDeBoor3WithDerivatives(x);
        for (int i = 4; i < 8; i++) { // to get actual derivatives, raw derivative values must be divided by bin width
            iota[i] /= epsT;
        }

        // in case when one of functional values is equal to maximum, fix starting indices:
        if (kFrom == nKappa - 3) {
            kFrom --;
        }
        if (jFrom == nIota - 3) {
            jFrom --;
        }

        // updating histogram of h and histogram of of ∂h/∂μ
        for (int k = 0; k < 4; k++) {
            int kNow = kFrom + k;
            for (int j = 0; j < 4; j++) {
                int jNow = jFrom + j;
                // update h
                h[kNow][jNow] += kappa[k]*iota[j];
                // update gradient of h
                double prod = kappa[k]*iota[4 + j];
                for (int i = 0; i < parameters; i++) {
                    dp[i][kNow][jNow] += prod*gradfT[i];
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
                    double factor = Math.log(h[k][j] / hT[j]);
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
