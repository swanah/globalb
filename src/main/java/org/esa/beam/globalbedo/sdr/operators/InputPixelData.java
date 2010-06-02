/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

/**
 *
 * @author akheckel
 */
public class InputPixelData {
    final PixelGeometry geom;
    final float surfPressure;
    final int nSpecWvl;
    final float[] specWvl;
    final float[] toaReflec;

    public InputPixelData(PixelGeometry geom, float surfPressure, float[] specWvl, float[] toaReflec) {
        this.geom = geom;
        this.surfPressure = surfPressure;
        this.specWvl = specWvl;
        this.nSpecWvl = specWvl.length;
        this.toaReflec = toaReflec;
    }

}
