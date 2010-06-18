/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.lutUtils.MomoLut;
import org.esa.beam.globalbedo.sdr.util.math.Powell;
import org.esa.beam.util.math.LookupTable;

/**
 *
 */
public class emodSpecTau implements UnivRetrievalFunction {

    private final MomoLut lut;
    private final LookupTable[] synLut;
    private final float[] lutAlbedo;
    //private final InputPixelData inPix;
    private final SourcePixelField spf;
    private final float[] soilSpec;
    private final float[] vegSpec;
    private final float[] specWeights;
    private final float llimit;
    private final float penalty;
    //private final float ndvi;
    private double[] pAtMin;
    private final boolean isSyn;
    //private final double rad2rfl;

    public emodSpecTau(SourcePixelField spf, MomoLut lut,
                       float[] soilSpec, float[] vegSpec, float[] specWeights) {
        this.lut = lut;
        this.synLut = null;
        this.lutAlbedo = null;
        this.isSyn = false;
        this.spf = spf;
        this.soilSpec = soilSpec;
        this.vegSpec  = vegSpec;
        this.specWeights = specWeights;

        // calculate NDVI
        // to make start point of optimization
        // dependent of NDVI
        //this.ndvi = calcNdvi(685.0f, 865.0f);

        // constants used in f()
        // determining constraints for the optimisation
        this.llimit = 5e-6f;
        this.penalty = 1000f;
        this.pAtMin = new double[]{-1, -1, -1};
        //this.rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));

    }

    public emodSpecTau(SourcePixelField spf, LookupTable[] lut, float[] lutAlbedo,
                       float[] soilSpec, float[] vegSpec, float[] specWeights) {
        this.lut = null;
        this.synLut = lut;
        this.lutAlbedo = lutAlbedo;
        this.isSyn = true;
        this.spf = spf;
        this.soilSpec = soilSpec;
        this.vegSpec  = vegSpec;
        this.specWeights = specWeights;

        // calculate NDVI
        // to make start point of optimization
        // dependent of NDVI
        //this.ndvi = calcNdvi(685.0f, 865.0f);

        // constants used in f()
        // determining constraints for the optimisation
        this.llimit = 5e-6f;
        this.penalty = 1000f;
        this.pAtMin = new double[]{-1, -1, -1};
        //this.rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));
    }


    @Override
    public double f(double tau) {
        float fmin = 0;
        for (int i=0; i<spf.inPixArr.length; i++){
            if (spf.valid[i]) fmin += fPix(tau, spf.inPixArr[i]);
        }
        return fmin;
    }

    @Override
    public float[] getModelReflec() {
        float[] model = new float[soilSpec.length];
        for (int i=0; i<model.length; i++) {
            model[i] = (float) (pAtMin[0] * vegSpec[i] + pAtMin[1] * soilSpec[i] + pAtMin[2]);
        }
        return model;
    }

    @Override
    public double[] getpAtMin() {
        return pAtMin;
    }

    @Override
    public float[] getSurfReflec(float aot) {
        return (isSyn) ? invertSynLut(spf.inPixArr[spf.pixelIndex], aot) : invertToaRefl(spf.inPixArr[spf.pixelIndex], aot);
    }

    private float[] getSurfReflec(InputPixelData inPix, float aot) {
        return (isSyn) ? invertSynLut(inPix, aot) : invertToaRefl(inPix, aot);
    }

    private double fPix (double tau, InputPixelData inPix){
        float fmin = 0;

        float[] surfReflec = getSurfReflec(inPix, (float) tau);

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
            float ndvi = calcNdvi(inPix, 685.0f, 865.0f);
            double[] pSpec = {ndvi, 1.0-ndvi, 0};

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

    private float calcNdvi(InputPixelData inPix, float ndviRed, float ndviNir) {
        int iRed = 0;
        int iNir = 0;
        for(int i=0; i<inPix.nSpecWvl; i++){
            if (Math.abs(inPix.specWvl[i]-ndviRed) < Math.abs(inPix.specWvl[iRed]-ndviRed)) iRed = i;
            if (Math.abs(inPix.specWvl[i]-ndviNir) < Math.abs(inPix.specWvl[iNir]-ndviNir)) iNir = i;
        }
        return (inPix.toaReflec[iNir] - inPix.toaReflec[iRed])
                / (inPix.toaReflec[iNir] + inPix.toaReflec[iRed]);
    }

    private float[] invertToaRefl(InputPixelData inPix, double tau) {
        float[] sdr = new float[inPix.nSpecWvl];
        double rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));
        float[] o3corr = lut.getO3corr();
        final double geomAMF = (1/Math.cos(Math.toRadians(inPix.geom.sza))+1/Math.cos(Math.toRadians(inPix.geom.vza)))/2;
        for (int i=0; i<inPix.nSpecWvl; i++){
            double[] x = {inPix.surfPressure,
                          inPix.geom.vza,
                          inPix.geom.sza,
                          inPix.geom.razi,
                          inPix.specWvl[i],
                          tau, 0};
            double rhoPath = lut.getLookupTable().getValue(x);
            x[6] = 1;
            double tupTdown = lut.getLookupTable().getValue(x);
            x[6] = 2;
            double spherAlb = lut.getLookupTable().getValue(x);
            double toaCorr = inPix.toaReflec[i] / (Math.exp(inPix.o3du * o3corr[i] * geomAMF));
            double a = (toaCorr - (rhoPath*rad2rfl)) / (tupTdown/Math.PI*rad2rfl);
            sdr[i] = (float) (a / (1 + spherAlb * a));
        }
        return sdr;
    }

    private float[] invertSynLut(InputPixelData inPix, double tau) {
        float[] sdr = new float[inPix.nSpecWvl];
        double rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));
        for (int i=0; i<inPix.nSpecWvl; i++){
            int iAlb = 0;
            double[] x1 = {Math.log(inPix.surfPressure),
                           inPix.geom.vza,
                           inPix.geom.razi,
                           inPix.geom.sza,
                           tau,
                           lutAlbedo[iAlb]};
            double toa1 = rad2rfl * synLut[i].getValue(x1);
            while (iAlb < lutAlbedo.length-1 && toa1 < inPix.toaReflec[i]){
                iAlb++;
                x1[5] = lutAlbedo[iAlb];
                toa1 = rad2rfl * synLut[i].getValue(x1);
            }
            if (iAlb <= 0) {
                iAlb = 1;
                x1[5] = lutAlbedo[iAlb];
                toa1 = rad2rfl * synLut[i].getValue(x1);
            }
            if (iAlb >= lutAlbedo.length) {
                iAlb = lutAlbedo.length-1;
                x1[5] = lutAlbedo[iAlb];
                toa1 = rad2rfl * synLut[i].getValue(x1);
            }
            double[] x0 = {Math.log(inPix.surfPressure),
                           inPix.geom.vza,
                           inPix.geom.razi,
                           inPix.geom.sza,
                           tau,
                           lutAlbedo[iAlb-1]};
            double toa0 = rad2rfl * synLut[i].getValue(x0);
            double alb0 = lutAlbedo[iAlb-1];
            double alb1 = lutAlbedo[iAlb];
            
            sdr[i] = (float) (alb0 + (inPix.toaReflec[i]-toa0)*(alb1-alb0)/(toa1-toa0));
        }
        return sdr;
    }

}
