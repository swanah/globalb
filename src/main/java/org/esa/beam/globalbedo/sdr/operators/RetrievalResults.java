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
    public final boolean retrievalFailed;
    public final float optAOT;
    public final float optErr;
    public final float retrievalErr;
    public float scaleVeg;
    public float scaleSoil;
    public float[] sdr;

    public RetrievalResults() {
        this(false, -1.0f, -1.0f, -1.0f);
    }

    public RetrievalResults(boolean retrievalFailed, float optAOT, float optErr, float retrievalErr) {
        this.retrievalFailed = retrievalFailed;
        this.optAOT = optAOT;
        this.optErr = optErr;
        this.retrievalErr = retrievalErr;
    }
}
