/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

/**
 *
 * @author akheckel
 */
public class RetrievalResults {
    final boolean retrievalFailed;
    final float optAOT;
    final float optErr;
    final float retrievalErr;
    final float curvature;
    double[] pAtMin;
    float[] modelSpec;
    double[] sdr;

    public RetrievalResults() {
        this(false, -1.0f, -1.0f, -1.0f, -1.0f);
    }

    public RetrievalResults(boolean retrievalFailed, float optAOT, float optErr, float retrievalErr) {
        this(retrievalFailed, optAOT, optErr, retrievalErr, -1.0f);
    }
    
    public RetrievalResults(boolean retrievalFailed, float optAOT, float optErr, float retrievalErr, float curv) {
        this.retrievalFailed = retrievalFailed;
        this.optAOT = optAOT;
        this.optErr = optErr;
        this.retrievalErr = retrievalErr;
        this.curvature = curv;
    }
}
