/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.util.math.MvFunction;

/**
 * This class provides the spectrum model function to be minimised by Powell.
 * (see ATBD (4), (6))
 *
 *
 */
public class emodSpec implements MvFunction {

    private final float[] specSoil;
    private final float[] specVeg;
    private final float[] surfReflec;
    private final float[] specWeights;
    private final int nSpecChannels;

    public emodSpec(float[] specSoil, float[] specVeg, float[] surfReflec, float[] specWeights) {
        this.specSoil = specSoil;
        this.specVeg = specVeg;
        this.surfReflec = surfReflec;
        this.specWeights = specWeights;
        this.nSpecChannels = surfReflec.length;
    }

    @Override
    public double f(double[] p) {
        double resid = 0.0;
        double[] mval = new double[nSpecChannels];
        float[] weight = normalize(specWeights);
        for (int iwvl = 0; iwvl < nSpecChannels; iwvl++) {
            // mval: rho_spec_mod in ATBD (p. 22) (model function)
            mval[iwvl] = p[0] * specVeg[iwvl] + p[1] * specSoil[iwvl] + p[2];
            // difference to measurement:
            double k = surfReflec[iwvl] - mval[iwvl];
            // residual:
            resid += weight[iwvl] * k * k;
        }

        // constraints for fit parameter p
        // specSoil and specVeg should not be scaled negative
        double limit;
        if (p[0] < 0.0) resid = resid + p[0] * p[0] * 1000;
        if (p[1] < 0.0) resid = resid + p[1] * p[1] * 1000;
        limit = -0.07;
        if (p[2] < limit) resid = resid + (p[2]-limit) * (p[2]-limit) * 1000;
        limit = -0.02;
        if (p[2] > limit) resid = resid + (p[2]-limit) * (p[2]-limit) * 1000;

        return(resid);
    }

    @Override
    public void g(double[] x, double[] g) {
    }

    private float[] normalize(float[] fa) {
    float sum = 0;
    for (float f : fa) sum += f;
    for (int i=0; i<fa.length; i++) fa[i] /= sum;
    return fa;
}

}

