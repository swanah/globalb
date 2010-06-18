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
    final float o3du;
    final int nSpecWvl;
    final float[] specWvl;
    final float[] toaReflec;

    public InputPixelData(PixelGeometry geom, float surfPressure, float o3du, float[] specWvl, float[] toaReflec) {
        this.geom = geom;
        this.surfPressure = surfPressure;
        this.o3du = o3du;
        this.specWvl = specWvl;
        this.nSpecWvl = specWvl.length;
        this.toaReflec = toaReflec;
    }

}
