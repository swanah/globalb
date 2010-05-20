/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.util.math.Brent;
import org.esa.beam.globalbedo.sdr.util.math.Function;
import org.esa.beam.globalbedo.sdr.util.math.MvFunction;
import org.esa.beam.globalbedo.sdr.util.math.Powell;
import org.esa.beam.util.math.LookupTable;

/**
 *
 * @author akheckel
 */
public class PointRetrieval {
    private final LookupTable lut;
    private final int nSpecChannels;
    private final float[] specSoil;     // [nSpecChannels]
    private final float[] specVeg;      // [nSpecChannels]
    private float[] surfReflec;
    private float[] toaReflec;
    private final float[] specWvl;
    //private final float[] geometry; //{SZA, RAZI, VZA, PRES}
    private float sza;
    private float razi;
    private float vza;
    private float surfPres;

    private final double[] specWeight;
    private final double ndviRed;
    private final double ndviNir;
    private double ndvi;
    private double aotStart;
    private double[] pSpec;
    private final Powell powell;

    public PointRetrieval(LookupTable lut, float[] soilSpec, float[] vegSpec, float[] specWvl) {
        this.lut = lut;
        this.nSpecChannels = specWvl.length;
        this.specSoil = soilSpec;
        this.specVeg = vegSpec;
        this.specWvl = specWvl;

        this.specWeight = new double[]{1.0, 1.0, 1.0, 1.0};
        this.ndviRed = 680.0;
        this.ndviNir = 880.0;
        this.aotStart = 0.1;
        this.pSpec = new double[]{0.7, 0.3};
        this.powell = new Powell();
    }


// public methods

    public RetrievalResults runRetrieval(float surfPres, float[] geometry, float[] toaReflec) {
        this.sza = geometry[0];
        this.razi = geometry[1];
        this.vza = geometry[2];
        this.surfPres = surfPres;
        this.toaReflec = toaReflec;
        this.pSpec = new double[]{0.8, 0.2};
        this.aotStart = 0.1;

        return runRetrieval();
    }

    public RetrievalResults runRetrieval() {
        boolean failed = false;

        final Brent b = new Brent(0.0, aotStart, 2.0, new emodSpecTau(), 5e-4);
        float optAOT = (float) b.getXmin();
        float optErr = (float) b.getFx();
        float retrievalErr = calcRetrievalErr();
        aotStart = optAOT;
        RetrievalResults results = new RetrievalResults(failed, optAOT, optErr, retrievalErr);
        results.scaleSoil = (float) pSpec[1];
        results.scaleVeg  = (float) pSpec[0];
        results.sdr = surfReflec;
        return results;
    }

// private methods

    private double calcNdvi() {
        int iRed = 0;
        int iNir = 0;
        for(int i=0; i<nSpecChannels; i++){
            if (Math.abs(specWvl[i]-ndviRed) < Math.abs(specWvl[iRed]-ndviRed)) iRed = i;
            if (Math.abs(specWvl[i]-ndviNir) < Math.abs(specWvl[iNir]-ndviNir)) iRed = i;
        }
        return (toaReflec[iNir] - toaReflec[iRed]) / (toaReflec[iNir] + toaReflec[iRed]);
    }

    private float calcRetrievalErr() {
        //TODO calcRetrievalError
        float error = 0;
        return error;
    }

    private float[] invertToaRefl(double tau) {
        //TODO invertToaRefl
        float[] sdr = new float[nSpecChannels];
        for (int i=0; i<nSpecChannels; i++){
            double[] x = {surfPres, vza, sza, razi, specWvl[i], tau, 0};
            double rhoPath = lut.getValue(x);
            x[6] = 1;
            double tupTdown = lut.getValue(x);
            x[6] = 2;
            double spherAlb = lut.getValue(x);
            double a = (toaReflec[i] - rhoPath) / tupTdown;
            sdr[i] = (float) (a / (1 + spherAlb * a));
        }
        return sdr;
    }

    
// private embedded classes

    /**
     * 
     */
    private class emodSpecTau implements Function {

        @Override
        public double f(double tau) {
            float fmin = 0;

            surfReflec = invertToaRefl(tau);

            // inversion can lead to overcorrection of atmosphere
            // and thus to too small surface reflectances
            // defining a steep but smooth function guides the optimization 

            final float llimit = 5e-6f;
            final float penalty = 1000f;
            for (int iwvl = 0; iwvl < nSpecChannels; iwvl++) {
                if (surfReflec[iwvl] < llimit) {
                    fmin += (surfReflec[iwvl]-llimit) * (surfReflec[iwvl]-llimit) * penalty;
                }
            }

            if (fmin <= 0.0f) {
                // initial vector p to start Powell optimization
                //pSpec[0] = ndvi; pSpec[1] = 1.0-ndvi;
                //if (pSpec.length > 2) pSpec[2] = 0.025;

                // defining unit matrix as base of the parameter space
                // needed for Powell
                double xi[][] = new double[pSpec.length][pSpec.length];
                for (int i = 0; i < pSpec.length; i++) xi[i][i] = 1.0;

                double ftol = 0.5e-2;   // change of fmin below ftol defines the end of optimization

                powell.powell(pSpec, xi, ftol, new emodSpec());
                fmin = (float) powell.fret;
            }
            else {
                //fmin += 1e-5;
                fmin += 1e-8;
            }
            return fmin;
        }
    }

    /**
     * This class provides the spectrum model function to be minimised by Powell.
     * (see ATBD (4), (6))
     *
     *
     */
    private class emodSpec implements MvFunction {

        public double f(double[] p) {
            double resid = 0.0;
            double[] mval = new double[nSpecChannels];
            double[] weight = normalize(specWeight);
            double k;
            for (int iwvl = 0; iwvl < nSpecChannels; iwvl++) {
                // mval: rho_spec_mod in ATBD (p. 22) (model function)
                if (p.length == 2) {
                    mval[iwvl] = p[0] * specVeg[iwvl] + p[1] * specSoil[iwvl];
                } else {
                    mval[iwvl] = p[0] * specVeg[iwvl] + p[1] * specSoil[iwvl] + p[2];
                }
                // difference to measurement:
                k = surfReflec[iwvl] - mval[iwvl];
                // residual:
                resid = resid + weight[iwvl] * k * k;
            }

            if (p[0] < 0.0) resid = resid + p[0] * p[0] * 1000;
            if (p[1] < 0.0) resid = resid + p[1] * p[1] * 1000;
            if (p.length > 2)
                if (p[2] < 0.0 || p[2] > 0.05) resid = resid + p[2] * p[2] * 1000;

            return(resid);
        }

        public void g(double[] x, double[] g) throws UnsupportedOperationException {
            throw new UnsupportedOperationException("Not supported yet.");
        }

        private double[] normalize(double[] fa) {
        float sum = 0;
        for (double f : fa) sum += f;
        for (int i=0; i<fa.length; i++) fa[i] /= sum;
        return fa;
    }

    }

}
