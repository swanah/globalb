/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.util.math.Powell;
import org.esa.beam.util.math.LookupTable;

/**
 *
 */
public class emodSpecTau implements UnivRetrievalFunction {

    private final LookupTable lut;
    private final InputPixelData inPix;
    private final float[] soilSpec;
    private final float[] vegSpec;
    private final float[] specWeights;
    private final float llimit;
    private final float penalty;
    private final float ndvi;
    private double[] pAtMin;

    public emodSpecTau(InputPixelData inPix, LookupTable lut,
                       float[] soilSpec, float[] vegSpec, float[] specWeights) {
        this.lut = lut;
        this.inPix = inPix;
        this.soilSpec = soilSpec;
        this.vegSpec  = vegSpec;
        this.specWeights = specWeights;

        // calculate NDVI
        // to make start point of optimization
        // dependent of NDVI
        this.ndvi = calcNdvi(685.0f, 865.0f);

        // constants used in f()
        // determining constraints for the optimisation
        this.llimit = 5e-6f;
        this.penalty = 1000f;
        this.pAtMin = null;
    }


    @Override
    public double f(double tau) {
        float fmin = 0;

        float[] surfReflec = invertToaRefl(tau);

        // inversion can lead to overcorrection of atmosphere
        // and thus to too small surface reflectances
        // defining a steep but smooth function guides the optimization
        for (int iwvl = 0; iwvl < inPix.nSpecWvl; iwvl++) {
            if (surfReflec[iwvl] < llimit) {
                fmin += (surfReflec[iwvl]-llimit) * (surfReflec[iwvl]-llimit) * penalty;
            }
        }

        if (fmin <= 0.0f) {
            // initial vector p to start Powell optimization
            
            double[] pSpec = {ndvi, 1.0-ndvi};

            // defining unit matrix as base of the parameter space
            // needed for Powell
            double xi[][] = new double[pSpec.length][pSpec.length];
            for (int i = 0; i < pSpec.length; i++) xi[i][i] = 1.0;

            double ftol = 0.5e-2;   // change of fmin below ftol defines the end of optimization
            // change of fmin below ftol defines the end of optimization

            emodSpec powellFitFct = new emodSpec(soilSpec, vegSpec, surfReflec, specWeights);

            Powell powell = new Powell(pSpec, xi, ftol, powellFitFct);
            fmin = (float) powell.getFmin();
            pAtMin = powell.getP();
        }
        else {
            //fmin += 1e-5;
            fmin += 1e-8;
        }
        return fmin;
    }

    @Override
    public float[] getModelReflec() {
        float[] model = new float[soilSpec.length];
        for (int i=0; i<model.length; i++) {
            model[i] = (float) (pAtMin[0] * vegSpec[i] + pAtMin[1] * soilSpec[i]);
        }
        return model;
    }

    @Override
    public double[] getpAtMin() {
        return pAtMin;
    }

    @Override
    public float[] getSurfReflec(float aot) {
        return invertToaRefl(aot);
    }


    private float calcNdvi(float ndviRed, float ndviNir) {
        int iRed = 0;
        int iNir = 0;
        for(int i=0; i<inPix.nSpecWvl; i++){
            if (Math.abs(inPix.specWvl[i]-ndviRed) < Math.abs(inPix.specWvl[iRed]-ndviRed)) iRed = i;
            if (Math.abs(inPix.specWvl[i]-ndviNir) < Math.abs(inPix.specWvl[iNir]-ndviNir)) iRed = i;
        }
        return (inPix.toaReflec[iNir] - inPix.toaReflec[iRed])
                / (inPix.toaReflec[iNir] + inPix.toaReflec[iRed]);
    }

    private float[] invertToaRefl(double tau) {
        float[] sdr = new float[inPix.nSpecWvl];
        for (int i=0; i<inPix.nSpecWvl; i++){
            double[] x = {inPix.surfPressure,
                          inPix.geom.vza,
                          inPix.geom.sza,
                          inPix.geom.razi,
                          inPix.specWvl[i],
                          tau, 0};
            double rhoPath = lut.getValue(x);
            x[6] = 1;
            double tupTdown = lut.getValue(x);
            x[6] = 2;
            double spherAlb = lut.getValue(x);
            double a = (inPix.toaReflec[i] - rhoPath) / tupTdown;
            sdr[i] = (float) (a / (1 + spherAlb * a));
        }
        return sdr;
    }

}