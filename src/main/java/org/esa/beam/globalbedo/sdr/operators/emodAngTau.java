/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.lutUtils.Aardvarc4DLut;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoLut;
import org.esa.beam.globalbedo.sdr.util.math.Powell;
import org.esa.beam.util.math.FracIndex;
import org.esa.beam.util.math.LookupTable;

/**
 *
 */
public class emodAngTau implements UnivRetrievalFunction {

    private final MomoLut lut;
    private final LookupTable[] synLut;
    private final Aardvarc4DLut aardvarcLut;
    private final float[] lutAlbedo;
    //private final InputPixelData inPix;
    private final SourcePixelField spf;
    private final double[] specWeights;
    private final float llimit;
    private final float penalty;
    //private final float ndvi;
    private double[] pAtMin;
    private final boolean isSyn;
    //private final double rad2rfl;

    public emodAngTau(SourcePixelField spf, MomoLut lut,
                      double[] specWeights) {
        this.lut = lut;
        this.synLut = null;
        this.aardvarcLut = null;
        this.lutAlbedo = null;
        this.isSyn = false;
        this.spf = spf;
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

    public emodAngTau(SourcePixelField spf, Aardvarc4DLut alut,
                      double[] specWeights) {
        this.lut = null;
        this.synLut = null;
        this.aardvarcLut = alut;
        this.lutAlbedo = null;
        this.isSyn = false;
        this.spf = spf;
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

    public emodAngTau(SourcePixelField spf, LookupTable[] lut, float[] lutAlbedo,
                      double[] specWeights) {
        this.lut = null;
        this.aardvarcLut = null;
        this.synLut = lut;
        this.lutAlbedo = lutAlbedo;
        this.isSyn = true;
        this.spf = spf;
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
    public float[] getModelReflec(float aot) {
        //float[][] diffuseFraction = getDiffuseFraction(spf.inPixArr[spf.pixelIndex], aot);
        InputPixelData ipd = spf.inPixArr[spf.pixelIndex];
        lut.getSdrAndDiffuseFrac(ipd, aot);
        int nSpecChannels = ipd.nSpecWvl;
        float[][] model = new float[2][nSpecChannels];
        double DF = 1.0f;
        double gamma = 0.35f;
        double dir, g, dif, k;

        for (int iwvl = 0; iwvl < nSpecChannels; iwvl++){
            for (int iview = 0; iview < 2; iview++) {
                dir = (1.0 - DF * ipd.diffuseFrac[iview][iwvl]) * pAtMin[nSpecChannels+iview] * pAtMin[iwvl];
                g   = (1.0 - gamma) * pAtMin[iwvl];
                dif = (DF * ipd.diffuseFrac[iview][iwvl]
                        + g * (1.0 - DF * ipd.diffuseFrac[iview][iwvl])) * gamma * pAtMin[iwvl] / (1.0 - g);
                // mval: rho_spec_ang in ATBD (p. 23) (model function)
                model[iview][iwvl] = (float) (dir + dif);
            }
        }
        return model[0];
    }

    @Override
    public double[] getpAtMin() {
        return pAtMin;
    }

    /**
     *
     * @param aot
     * @return nadir sdr only
     */
    @Override
    public double[] getSurfReflec(float aot) {
        double[] sdr = new double[spf.inPixArr[spf.pixelIndex].toaReflec.length];
        if (isSyn) {
            return invertSynLut(spf.inPixArr[spf.pixelIndex], aot);
        }else{
            InputPixelData ipd = spf.inPixArr[spf.pixelIndex];
            lut.getSdrAndDiffuseFrac(ipd, aot);
            sdr = ipd.surfReflec[0];
            //float[][] sdrDiff = invertToaRefl(spf.inPixArr[spf.pixelIndex], aot);
            //sdr = sdrDiff[0];
        }
        return sdr;
    }

    private double[][] getSurfReflec(InputPixelData inPix, float aot) {
        //return (isSyn) ? invertSynLut(inPix, aot) : invertToaRefl(inPix, aot);
        //return invertToaRefl(inPix, aot);
        lut.getSdrAndDiffuseFrac(inPix, aot);
        return inPix.surfReflec;
    }

    private double fPix (double tau, InputPixelData inPix){
        float fmin = 0;

        if (aardvarcLut != null) {
            aardvarcLut.getSdrAndDiffuseFrac(inPix, tau);
        }
        else {
            lut.getSdrAndDiffuseFrac(inPix, tau);
        }
        double[][] sdr = inPix.surfReflec;
        double[][] diffFrac = inPix.diffuseFrac;


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

            emodAng powellFitFct = new emodAng(diffFrac, sdr, specWeights);

            Powell powell = new Powell(pAng, xi, ftol, powellFitFct);
            fmin = (float) powell.getFmin();
            pAtMin = powell.getP();
        }
        else {
            //fmin += 1e-5;
            fmin += 1e-8;
        }
        return fmin;
    }

/*
    private float[][] invertToaRefl(InputPixelData inPix, double tau) {
        float[][] sdr = new float[2][inPix.nSpecWvl];
        double rad2rfl = Math.PI / Math.cos(Math.toRadians(inPix.geom.sza));

        //nadir
        float geomAMF = (float) ((1 / Math.cos(Math.toRadians(inPix.geom.sza))
                                        + 1 / Math.cos(Math.toRadians(inPix.geom.vza))));
        //float[] o3corr = lut.getO3corr();
        float[] gasT = lut.getGasTransmission(geomAMF, inPix.wvCol, inPix.o3du/1000);
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
            //double tgO3 = Math.exp(inPix.o3du * o3corr[i] * geomAMF/2); //my o3 correction scheme uses AMF=SC/VC not AMF=SC
            double toaCorr = inPix.toaReflec[i] / gasT[i];
            double a = (toaCorr - (rhoPath*rad2rfl)) / (tupTdown/Math.PI*rad2rfl);
            sdr[0][i] = (float) (a / (1 + spherAlb * a));
        }

        //fward
        if (inPix.geomFward != null){
            geomAMF = (float) ((1 / Math.cos(Math.toRadians(inPix.geomFward.sza))
                                            + 1 / Math.cos(Math.toRadians(inPix.geomFward.vza))));
            //float[] o3corr = lut.getO3corr();
            gasT = lut.getGasTransmission(geomAMF, inPix.wvCol, inPix.o3du/1000);
            for (int i=0; i<inPix.nSpecWvl; i++){
                double[] x = {inPix.surfPressure,
                              inPix.geomFward.vza,
                              inPix.geomFward.sza,
                              inPix.geomFward.razi,
                              inPix.specWvl[i],
                              tau, 0};
                double rhoPath = lut.getLookupTable().getValue(x);
                x[6] = 1;
                double tupTdown = lut.getLookupTable().getValue(x);
                x[6] = 2;
                double spherAlb = lut.getLookupTable().getValue(x);
                //double tgO3 = Math.exp(inPix.o3du * o3corr[i] * geomAMF/2); //my o3 correction scheme uses AMF=SC/VC not AMF=SC
                double toaCorr = inPix.toaReflecFward[i] / gasT[i];
                double a = (toaCorr - (rhoPath*rad2rfl)) / (tupTdown/Math.PI*rad2rfl);
                sdr[1][i] = (float) (a / (1 + spherAlb * a));
            }
        }
        return sdr;
    }
*/
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

/*
 private float[][] getDiffuseFraction(InputPixelData inPix, double tau) {
        float[][] result = new float[2][inPix.nSpecWvl];
        for (int i=0; i<inPix.nSpecWvl; i++){
            double[] x = {inPix.surfPressure,
                          inPix.geom.vza,
                          inPix.geom.sza,
                          inPix.geom.razi,
                          inPix.specWvl[i],
                          tau, 3};
            double diff2total = 1.0 - lut.getLookupTable().getValue(x);
            result[0][i] = (float) diff2total;
            x = new double[] {inPix.surfPressure,
                          inPix.geomFward.vza,
                          inPix.geomFward.sza,
                          inPix.geomFward.razi,
                          inPix.specWvl[i],
                          tau, 3};
            diff2total = 1.0 - lut.getLookupTable().getValue(x);
            result[1][i] = (float) diff2total;
        }
        return result;
    }
 *
 */

}
