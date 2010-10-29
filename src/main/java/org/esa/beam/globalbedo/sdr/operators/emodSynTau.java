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
 * @author akheckel
 */
class emodSynTau implements UnivRetrievalFunction{

    private final MomoLut lut;
    private final LookupTable[] synLut;
    private final float[] lutAlbedo;
    //private final InputPixelData inPix;
    private final SourcePixelField spf;
    private final float[] soilSpec;
    private final float[] vegSpec;
    private final double[] specWeights;
    private final float llimit;
    private final float penalty;
    //private final float ndvi;
    private double[] pAngAtMin;
    private double[] pSpecAtMin;
    private final boolean isSyn;
    //private final double rad2rfl;

    public emodSynTau(SourcePixelField spf, MomoLut lut,
                      float[] soilSpec, float[] vegSpec, double[] specWeights) {
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
        this.pAngAtMin = new double[]{-1, -1, -1};
        //this.rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));

    }

    public emodSynTau(SourcePixelField spf, LookupTable[] lut, float[] lutAlbedo,
                      float[] soilSpec, float[] vegSpec, double[] specWeights) {
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
        this.pAngAtMin = new double[]{-1, -1, -1};
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
    public float[] getModelReflec(float aot) {
        double[][] diffuseFraction = getDiffuseFraction(spf.inPixArr[spf.pixelIndex], aot);
        int nSpecChannels = spf.inPixArr[0].nSpecWvl;
        float[][] model = new float[2][nSpecChannels];
        double DF = 1.0f;
        double gamma = 0.35f;
        double dir, g, dif, k;

        for (int iwvl = 0; iwvl < nSpecChannels; iwvl++){
            for (int iview = 0; iview < 2; iview++) {
                dir = (1.0 - DF * diffuseFraction[iview][iwvl]) * pAngAtMin[nSpecChannels+iview] * pAngAtMin[iwvl];
                g   = (1.0 - gamma) * pAngAtMin[iwvl];
                dif = (DF * diffuseFraction[iview][iwvl]
                        + g * (1.0 - DF * diffuseFraction[iview][iwvl])) * gamma * pAngAtMin[iwvl] / (1.0 - g);
                // mval: rho_spec_ang in ATBD (p. 23) (model function)
                model[iview][iwvl] = (float) (dir + dif);
            }
        }
        return model[0];
    }

    @Override
    public double[] getpAtMin() {
        return pAngAtMin;
    }

    @Override
    public double[] getSurfReflec(float aot) {
        double[] sdr = new double[spf.inPixArr[spf.pixelIndex].toaReflec.length];
        if (isSyn) {
            return invertSynLut(spf.inPixArr[spf.pixelIndex], aot);
        }else{
            double[][] sdrNadFwd = invertToaRefl(spf.inPixArr[spf.pixelIndex], aot);
            sdr = sdrNadFwd[0];
        }
        return sdr;
    }

    private double[][] getSurfReflec(InputPixelData inPix, float aot) {
        //return (isSyn) ? invertSynLut(inPix, aot) : invertToaRefl(inPix, aot);
        return invertToaRefl(inPix, aot);
    }

    private double fPix (double tau, InputPixelData inPix){
        float fmin = 0;

        double[][] diffFrac = getDiffuseFraction(inPix, tau);
        double[][] sdr = getSurfReflec(inPix, (float) tau);

        // inversion can lead to overcorrection of atmosphere
        // and thus to too small surface reflectances
        // defining a steep but smooth function guides the optimization
        for (int iwvl = 0; iwvl < inPix.nSpecWvl; iwvl++) {
            if (sdr[0][iwvl] < llimit) {
                fmin += (sdr[0][iwvl]-llimit) * (sdr[0][iwvl]-llimit) * penalty;
            }
            if (sdr[1][iwvl] < llimit) {
                fmin += (sdr[1][iwvl]-llimit) * (sdr[1][iwvl]-llimit) * penalty;
            }
        }

        if (fmin <= 0.0f) {
            // initial vector p to start Powell optimization
            double[] pAng = {0.1f, 0.1f, 0.1f, 0.1f, 0.5f, 0.3f};

            // defining unit matrix as base of the parameter space
            // needed for Powell
            double xi[][] = new double[pAng.length][pAng.length];
            for (int i = 0; i < pAng.length; i++) xi[i][i] = 1.0;

            double ftol = 0.5e-2;   // change of fmin below ftol defines the end of optimization
            // change of fmin below ftol defines the end of optimization

            emodAng powellAngFct = new emodAng(diffFrac, sdr, specWeights);

            Powell powell = new Powell(pAng, xi, ftol, powellAngFct);
            float fminAng = (float) powell.getFmin();
            pAngAtMin = powell.getP();

            double[] pSpec = {0.05f, 0.95f, 0f};
            xi = new double[pSpec.length][pSpec.length];
            for (int i = 0; i < pSpec.length; i++) xi[i][i] = 1.0;
            emodSpec powellSpecFct = new emodSpec(soilSpec, vegSpec, sdr[1], specWeights);

            powell = new Powell(pSpec, xi, ftol, powellSpecFct);
            float fminSpec = (float) powell.getFmin();
            pSpecAtMin = powell.getP();

            float ndvi = calcNdvi(inPix, 685.0f, 865.0f);
            float weight = (ndvi > 0 && ndvi < 0.7) ? ndvi/0.7f*0.5f : 0.5f;
            fmin = (1.0f-weight)*fminAng + weight*fminSpec;
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
        return (float)((inPix.toaReflec[iNir] - inPix.toaReflec[iRed])
                / (inPix.toaReflec[iNir] + inPix.toaReflec[iRed]));
    }

    private double[][] invertToaRefl(InputPixelData inPix, double tau) {
        lut.getSdrAndDiffuseFrac(inPix, tau);
        return inPix.surfReflec;
    }

    private double[] invertSynLut(InputPixelData inPix, double tau) {
        double[] sdr = new double[inPix.nSpecWvl];
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

            sdr[i] = alb0 + (inPix.toaReflec[i]-toa0)*(alb1-alb0)/(toa1-toa0);
        }
        return sdr;
    }

    private double[][] getDiffuseFraction(InputPixelData inPix, double tau) {
        lut.getSdrAndDiffuseFrac(inPix, tau);
        return inPix.diffuseFrac;
    }
}
