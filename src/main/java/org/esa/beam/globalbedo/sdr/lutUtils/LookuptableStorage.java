/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.lutUtils;


import java.io.Serializable;
import org.esa.beam.util.math.LookupTable;

/**
 * Serializable Class used to store Lookuptables
 * @author akheckel
 */
public class LookuptableStorage implements Serializable {
    private final static long serialVersionUID = 25041973L;
    private final float[] values;
    private final float[][] dimensions;

    /**
     * empty Standart consturctor
     */
    public LookuptableStorage() {
        this.values = null;
        this.dimensions = null;
    }
    
    /**
     * Constructor for the lookup table
     * @param values - Values of the lookuptable
     * @param dimensions - sequence of 1d arrays providing the dimensions of the lookup table
     */
    public LookuptableStorage(float[] values, float[]... dimensions) {
        this.values = values;
        this.dimensions = dimensions;
    }

    public float[] getValues() {
        return values;
    }

    public float[][] getDimensions() {
        return dimensions;
    }

    public LookupTable getLookupTable() {
        LookupTable lut = new LookupTable(values, dimensions);
        return lut;
    }
}
