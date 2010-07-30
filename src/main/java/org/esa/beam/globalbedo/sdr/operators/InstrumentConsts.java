/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.dataio.envisat.EnvisatConstants;
import org.esa.beam.framework.gpf.OperatorException;

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
//    private final float[] merisFitWeights = {1.0f, 1.0f, 1.0f, 1.0f, 0.2f, 1.0f, 1.0f, 1.0f,
//                                             0.05f, 0.05f, 0.05f, 0.05f, 0.05f};
    private final float[] merisFitWeights = {1.0f, 1.0f, 1.0f, 1.0f, 0.2f, 1.0f, 1.0f, 1.0f,
                                             0.5f, 0.5f, 0.5f, 0.5f, 0.5f};
    private final String merisValidExpr = "(!l1_flags.INVALID)";
    private final int merisNLutBands = 15;
    private final String merisSurfPressureName = "atm_press";
    private final String merisOzoneName = "ozone";

    private final String[] vgtReflectanceNames = {"B0", "B2", "B3", "MIR"};
    private final String[] vgtGeomNames = {"SZA", "SAA", "VZA", "VAA"};
    private final float[] vgtFitWeights = {1.0f, 1.0f, 1.0f, 1.0f};
    private final String  vgtCloudExpr = "( MIR<180*0.0005 )";
    //private final String  vgtValidExpr = "(SM.B0_GOOD && SM.B2_GOOD && SM.B3_GOOD && SM.MIR_GOOD && SM.LAND && " + vgtCloudExpr + ")";
    private final String  vgtValidExpr = "(SM.B0_GOOD && SM.B2_GOOD && SM.B3_GOOD && SM.MIR_GOOD && " + vgtCloudExpr + ")";
    private final int vgtNLutBands = 4;
    private final String vgtSurfPressureName = "surfPressEstimate";
    private final String vgtOzoneName = "OG";

    private final Map<String, String[]> reflecNames;
    private final Map<String, String[]> geomNames;
    private final Map<String, float[]> fitWeights;
    private final Map<String, String> validExpr;
    private final Map<String, Integer> nLutBands;
    private final Map<String, String> surfPressureName;
    private final Map<String, String> ozoneName;

    private final String lutLocaFile = System.getProperty("user.home")
                    + File.separator + ".beam"
                    + File.separator + "ga-aerosol"
                    + File.separator + "lut.location";
    private final String lutPattern = "%INSTRUMENT%/%INSTRUMENT%_LUT_MOMO_ContinentalI_80_SDR_noG.bin";
    //private final String lutPattern = "%INSTRUMENT%/%INSTRUMENT%_LUT_MOMO_ContinentalI_80_SU_noG_Kx-AOD.bin";
    private final String lutPath;
    //private final String lutPath = "e:/model_data/momo/LUTs_Swansea/%INSTRUMENT%/%INSTRUMENT%_LUT_MOMO_ContinentalI_80_SU_noG.bin";

    private InstrumentConsts() {
        lutPath = getLutPath() + "/" + lutPattern;

        this.supportedInstruments = new String[]{"MERIS", "VGT"};

        this.reflecNames = new HashMap<String, String[]>(supportedInstruments.length);
        reflecNames.put(supportedInstruments[0], merisReflectanceNames);
        reflecNames.put(supportedInstruments[1], vgtReflectanceNames);

        this.geomNames = new HashMap<String, String[]>(supportedInstruments.length);
        geomNames.put(supportedInstruments[0], merisGeomNames);
        geomNames.put(supportedInstruments[1], vgtGeomNames);

        this.fitWeights = new HashMap<String, float[]>(supportedInstruments.length);
        fitWeights.put(supportedInstruments[0], merisFitWeights);
        fitWeights.put(supportedInstruments[1], vgtFitWeights);

        this.validExpr = new HashMap<String, String>(supportedInstruments.length);
        validExpr.put(supportedInstruments[0], merisValidExpr);
        validExpr.put(supportedInstruments[1], vgtValidExpr);

        this.nLutBands = new HashMap<String, Integer>(supportedInstruments.length);
        nLutBands.put(supportedInstruments[0], merisNLutBands);
        nLutBands.put(supportedInstruments[1], vgtNLutBands);

        this.surfPressureName = new HashMap<String, String>(supportedInstruments.length);
        surfPressureName.put(supportedInstruments[0], merisSurfPressureName);
        surfPressureName.put(supportedInstruments[1], vgtSurfPressureName);

        this.ozoneName = new HashMap<String, String>(supportedInstruments.length);
        ozoneName.put(supportedInstruments[0], merisOzoneName);
        ozoneName.put(supportedInstruments[1], vgtOzoneName);
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

    public String getSurfPressureName(String instrument) {
        return surfPressureName.get(instrument);
    }

    public String getOzoneName(String instrument) {
        return ozoneName.get(instrument);
    }

    private String getLutPath() {
        String lutP = null;
        BufferedReader reader = null;
        try{
            reader = new BufferedReader(new FileReader(lutLocaFile));
            String line;
            while ((line = reader.readLine())!= null) {
                if (line.startsWith("ga.lutInstallDir")) {
                    String[] split = line.split("=");
                    if(split.length > 1) {
                        lutP = split[1].trim();
                        break;
                    }
                }
            }
        } catch (IOException ex) {
            throw new OperatorException(ex);
        }
        if (lutP == null) throw new OperatorException("Lut install dir  not found");
        if (lutP.endsWith(File.separator) || lutP.endsWith("/")) lutP = lutP.substring(0, lutP.length()-1);
        return lutP;
    }

}
