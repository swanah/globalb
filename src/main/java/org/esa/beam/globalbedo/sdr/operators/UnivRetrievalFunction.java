/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import org.esa.beam.globalbedo.sdr.util.math.Function;

/**
 *
 * @author akheckel
 */
public interface UnivRetrievalFunction extends Function {

    float[] getModelReflec();

    double[] getpAtMin();

    float[] getSurfReflec(float aot);

}
