/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

/**
 *
 * @author akheckel
 */
public class SourcePixelField {
    public final InputPixelData[] inPixArr;
    public final boolean[] valid;
    //TODO: the concept of pixelIndex only works if the pixel to be retrieved is inside the pixelWindow!!!
    public final int pixelIndex;

    public SourcePixelField(InputPixelData[] inPixArr, boolean[] valid, int pixIndex) {
        this.inPixArr = inPixArr;
        this.valid = valid;
        this.pixelIndex = pixIndex;
    }

}
