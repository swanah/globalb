/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import Jama.Matrix;
import org.esa.beam.globalbedo.sdr.util.math.Brent;
import org.esa.beam.globalbedo.sdr.util.math.Function;

/**
 * Provides aerosol retrieval class
 * currently not thread safe !!!
 *     subsequent calls to Brent and Powell.
 *     Somewhere along that line thread safety is broken
 * --> instantiate locally in computeTileStack()
 * 
 * @author akheckel
 */
public class PointRetrieval {

    private final UnivRetrievalFunction brentFitFct;

    public PointRetrieval(UnivRetrievalFunction brentFitFct) {
        this.brentFitFct = brentFitFct;
    }


// public methods

    public RetrievalResults runRetrieval(String instrument) {
        final Brent b = new Brent(0.005, 0.1, 2.0, brentFitFct, 5e-4);
        float optAOT = (float) b.getXmin();
        float optErr = (float) b.getFx();
        float retrievalErr = calcRetrievalErr(optAOT, optErr);
        boolean failed = optAOT < 0.005;
        failed = failed || (retrievalErr/optAOT > 3);
        RetrievalResults results = new RetrievalResults(failed, optAOT, optErr, retrievalErr);
/*
        if (instrument.equals("VGT")){
            results.sdr = brentFitFct.getSurfReflec(optAOT);
            results.modelSpec = brentFitFct.getModelReflec(optAOT);
            results.pAtMin = brentFitFct.getpAtMin();
        }
 */
        return results;
    }

    public RetrievalResults retrieveSDR(double aot) {
        boolean failed = false;
        float optAOT = (float) aot;
        float optErr = (float) brentFitFct.f(aot);
        float retrievalErr = calcRetrievalErr(optAOT, optErr);
        RetrievalResults results = new RetrievalResults(failed, optAOT, optErr, retrievalErr);
/*
        results.sdr = brentFitFct.getSurfReflec(optAOT);
        results.modelSpec = brentUnivFct.getModelReflec(optAOT);
        results.pAtMin = brentUnivFct.getpAtMin();
 * 
 */
        return results;
    }

// private methods

    private float calcRetrievalErr(float optAOT, float optErr) {

        //backup solution for unphysical zero aot retrievals
        double p1 = (optAOT > 0.05) ? 0.8 : 1.0/0.8;
        double p2 = (optAOT > 0.05) ? 0.6 : 1.0/0.6;

        final double[] x0 = {1, 1, 1};
        final double[] x1 = {p1*optAOT, optAOT, p2*optAOT};
        final double[] x2 = {x1[0]*x1[0], x1[1]*x1[1], x1[2]*x1[2]};

        double optErrLow = brentFitFct.f(x1[0]);
        double optErrHigh = brentFitFct.f(x1[2]);
        final double[][] y = {{optErrLow}, {optErr}, {optErrHigh}};
        final double[][] xArr = {{x2[0], x1[0], x0[0]},{x2[1], x1[1], x0[1]},{x2[2], x1[2], x0[2]}};
        Matrix A = new Matrix(xArr);
        Matrix c = new Matrix(y);
        Matrix result = A.solve(c);
        final double[][] resultArr = result.getArray();
        final double a = resultArr[0][0]; // curvature term of parabola

        double retrievalError;
        if (a < 0) {
            retrievalError = Math.sqrt(optErr / 0.8 * 2 / 1e-4) + 0.03;
            //failed = true;
        }
        else {
            retrievalError = Math.sqrt(optErr / 0.8 * 2 / a) + 0.03;
        }
        return (float) retrievalError;
    }

}
