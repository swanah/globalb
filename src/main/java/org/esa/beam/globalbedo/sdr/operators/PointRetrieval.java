/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.util.math.Brent;
import ucar.nc2.iosp.netcdf3.SPFactory;

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

// public methods

    public RetrievalResults runRetrieval(UnivRetrievalFunction brentFitFct) {
        boolean failed = false;

        final Brent b = new Brent(0.0, 0.1, 2.0, brentFitFct, 5e-4);
        float optAOT = (float) b.getXmin();
        float optErr = (float) b.getFx();
        float retrievalErr = calcRetrievalErr();
        RetrievalResults results = new RetrievalResults(failed, optAOT, optErr, retrievalErr);
        results.sdr = brentFitFct.getSurfReflec(optAOT);
        results.modelSpec = brentFitFct.getModelReflec();
        results.pAtMin = brentFitFct.getpAtMin();
        return results;
    }

    public RetrievalResults retrieveSDR(UnivRetrievalFunction brentFitFct, double aot) {
        boolean failed = false;
        float optAOT = (float) aot;
        float optErr = (float) brentFitFct.f(aot);
        float retrievalErr = calcRetrievalErr();
        RetrievalResults results = new RetrievalResults(failed, optAOT, optErr, retrievalErr);
        results.sdr = brentFitFct.getSurfReflec(optAOT);
        results.modelSpec = brentFitFct.getModelReflec();
        results.pAtMin = brentFitFct.getpAtMin();
        return results;
    }

// private methods

    private float calcRetrievalErr() {
        //TODO calcRetrievalError
        float error = 0;
        return error;
    }

}
