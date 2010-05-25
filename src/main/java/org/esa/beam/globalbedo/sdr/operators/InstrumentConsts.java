/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.dataio.envisat.EnvisatConstants;
import org.esa.beam.meris.brr.Rad2ReflOp;

/**
 * Instrument specific constants
 * @author akheckel
 */
public class InstrumentConsts {

    private static InstrumentConsts instance;

    private final String[] supportedInstruments;

    // band names

    private final String[] merisReflectanceNames = {
        "reflectance_1",
        "reflectance_2",
        "reflectance_3",
        "reflectance_4",
        "reflectance_5",
        "reflectance_6",
        "reflectance_7",
        "reflectance_8",
        "reflectance_9",
        "reflectance_10",
        "reflectance_12",
        "reflectance_13",
        "reflectance_14",
    };
    private final String[] merisGeomNames = {
        EnvisatConstants.MERIS_SUN_ZENITH_DS_NAME,
        EnvisatConstants.MERIS_SUN_AZIMUTH_DS_NAME,
        EnvisatConstants.MERIS_VIEW_ZENITH_DS_NAME,
        EnvisatConstants.MERIS_VIEW_AZIMUTH_DS_NAME
    };
    private final float[] merisFitWeights = {1.0f, 1.0f, 1.0f, 1.0f, 0.2f, 1.0f, 1.0f, 1.0f,
                                             0.05f, 0.05f, 0.05f, 0.05f, 0.05f};
    private final String merisValidExpr = "(!l1_flags.INVALID)";
    private final int merisNLutBands = 15;

    private final String[] vgtReflectanceNames = {"B0", "B2", "B3", "MIR"};
    private final String[] vgtGeomNames = {"SZA", "SAA", "VZA", "VAA"};
    private final float[] vgtFitWeights = {1.0f, 1.0f, 1.0f, 1.0f};
    private final String  vgtValidExpr = "(SM.B0_OK && SM.B2_OK && SM.B3_OK && SM.MIR_OK && SM.LAND && (MIR > 0.045))";
    private final int vgtNLutBands = 4;

    private final Map<String, String[]> reflecNames;
    private final Map<String, String[]> geomNames;
    private final Map<String, float[]> fitWeights;
    private final Map<String, String> validExpr;
    private final Map<String, Integer> nLutBands;

    private final String lutPath = "e:/model_data/momo/LUTs_Swansea/%INSTRUMENT%/%INSTRUMENT%_LUT_MOMO_ContinentalI_80_SU.bin";

    private InstrumentConsts() {
        this.supportedInstruments = new String[]{"MERIS", "VGT"};

        this.reflecNames = new HashMap<String, String[]>();
        reflecNames.put(supportedInstruments[0], merisReflectanceNames);
        reflecNames.put(supportedInstruments[1], vgtReflectanceNames);

        this.geomNames = new HashMap<String, String[]>();
        geomNames.put(supportedInstruments[0], merisGeomNames);
        geomNames.put(supportedInstruments[1], vgtGeomNames);

        this.fitWeights = new HashMap<String, float[]>();
        fitWeights.put(supportedInstruments[0], merisFitWeights);
        fitWeights.put(supportedInstruments[1], vgtFitWeights);

        this.validExpr = new HashMap<String, String>();
        validExpr.put(supportedInstruments[0], merisValidExpr);
        validExpr.put(supportedInstruments[1], vgtValidExpr);

        this.nLutBands = new HashMap<String, Integer>();
        nLutBands.put(supportedInstruments[0], merisNLutBands);
        nLutBands.put(supportedInstruments[1], vgtNLutBands);
    }



    public static InstrumentConsts getInstance() {
        if (instance == null) {
            instance = new InstrumentConsts();
        }
        return instance;
    }

    public float[] getFitWeights(String instrument) {
        return fitWeights.get(instrument);
    }

    public String[] getGeomBandNames(String instrument) {
        return geomNames.get(instrument);
    }

    public String getLutName(String instrument) {
        return lutPath.replace("%INSTRUMENT%", instrument);
    }

    public String[] getSpecBandNames(String instrument) {
        return reflecNames.get(instrument);
    }

    public String[] getSupportedInstruments() {
        return supportedInstruments;
    }

    public String getValidExpression(String instrument) {
        return validExpr.get(instrument);
    }

    public int getnLutBands(String instrument) {
        return nLutBands.get(instrument);
    }

}
