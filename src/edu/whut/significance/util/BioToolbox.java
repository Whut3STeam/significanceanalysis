package edu.whut.significance.util;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.util.MathArrays;

import java.util.Arrays;

/**
 * Created by Justin on 2017/6/2.
 */
public class BioToolbox {
    /**
     * Gaussian blur.
     *
     * @param x     the x
     * @param width the width
     * @param sigma the sigma
     * @return the double [ ]
     */
    public static double[] GaussianBlur(double[] x, int width, double sigma) {

        int len = width | 0x01;
        int halflen = len >> 1;

        double[] h = new double[len];

        Gaussian g = new Gaussian(0, sigma);
        double step = 3.0 / halflen;
        double pos = step;

        h[halflen] = g.value(0);
        for (int i = 1; i <= halflen; i++, pos += step) {
            double v = g.value(pos);
            h[halflen + i] = v;
            h[halflen - i] = v;
        }

        h = MathArrays.normalizeArray(h, 1);

        double[] result = MathArrays.convolve(x, h);

        return Arrays.copyOfRange(result, halflen, x.length + halflen);

    }

    /**
     * Blur double [ ].
     *
     * @param x     the x
     * @param width the width of window
     * @return the double [ ]
     */
    public static double[] blur(double[] x, int width) {

        int len = width | 0x01; //make sure the width of window is odd.
        int halfLen = len >> 1;

        double[] h = new double[len];
        double f = 1.0 / len;
        for (int i = 0; i < len; i++) {
            h[i] = f;
        }

        double[] result = MathArrays.convolve(x, h);

        return Arrays.copyOfRange(result, halfLen, x.length + halfLen);

    }
}
