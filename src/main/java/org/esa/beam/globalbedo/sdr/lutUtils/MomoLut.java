/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.lutUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.util.math.LookupTable;

/**
 * class to read the lookuptables for the atmospheric correction
 * for GlobAlbedo provided by L.Guanter, FUB, Berlin
 * The LUTs contain values of the following parameters:
 *  1 - atmospheric path radiance
 *  2 - Tdown * Tup = product of downward and upward atmospheric transmission
 *  3 - spherical albedo of the atmosphere at ground level
 *  4 - fraction of direct to diffuse radiation at ground level
 *
 * @author akheckel
 */
public class MomoLut {
    private final int nWvl;
    private final int nParameter;
    private final int nVza;
    private final int nSza;
    private final int nAzi;
    private final int nHsf;
    private final int nAot;
    private final float[] vza;
    private final float[] sza;
    private final float[] azi;
    private final float[] hsf;
    private final float[] aot;
    private final float[] wvl;
    private final float[] values;
    private final float solIrrad;
    private final LookupTable lut;
    private final float[] o3corr;


    /**
     * standart constructor reading the binary LUT file from "lutName"
     * the number of channels or wavelength for which the LUTs is given
     * is not contained in the file
     * @param instrument - instrument name
     * @param lutName - file name of the binary LUTs (original format from FUB)
     * @param nWvl - number of spectral channels
     */
    public MomoLut(String lutName, int nWvl) {
        this.nParameter = 4;  // the 4 parameter i the LUT as described above
        this.nWvl = nWvl;

        ByteBuffer bb = readFileToByteBuffer(lutName);

        // read LUT dimensions and values
        this.vza = readDimension(bb); this.nVza = vza.length;
        this.sza = readDimension(bb); this.nSza = sza.length;
        this.azi = readDimension(bb); this.nAzi = azi.length;
        this.hsf = readDimension(bb); this.nHsf = hsf.length;
        this.aot = readDimension(bb); this.nAot = aot.length;

        // attention readValues depends on the n??? fields !!!
        this.values = readValues(bb);

        this.wvl = readDimension(bb, nWvl);
        this.solIrrad = bb.getFloat();

        // invert order to store stirctly increasing dimension in lut
        // swapping the actual values in the lut is taken care off in readValues()
        for (int i=0; i<nHsf/2; i++){
            float swap = hsf[nHsf-1-i];
            hsf[nHsf-1-i] = hsf[i];
            hsf[i] = swap;
        }
        this.lut = new LookupTable(values, getDimensions());
        this.o3corr = readO3corr();
    }

    // public methods

    public float[][] getDimensions() {
        return new float[][]{hsf, vza, sza, azi, wvl, aot, new float[]{0,1,2,3}};
    }

    public LookupTable getLookupTable() {
        return this.lut;
    }

    public float[] getValues() {
        return values;
    }

    public float getSolIrrad() {
        return solIrrad;
    }

    public float[] getAot() {
        return aot;
    }

    public float[] getAzi() {
        return azi;
    }

    public float[] getHsf() {
        return hsf;
    }

    public float[] getO3corr() {
        return o3corr;
    }

    public float[] getSza() {
        return sza;
    }

    public float[] getVza() {
        return vza;
    }

    public float[] getWvl() {
        return wvl;
    }


    // private methods

    private int calcPosition(int[] indices, int[] sizes) {
        //(((((iHsf*nVza + iVza)* nSza + iSza)* nAzi + iAzi)* nWvl + iWvl)* nAot + iAot)* nParameter + iPar
        int pos = 0;
        for (int i=0; i<sizes.length; i++) pos = (pos * sizes[i] + indices[i]);
        return pos;
    }

    private float[] readDimension(ByteBuffer bb) {
        int len = bb.getInt();
        return readDimension(bb, len);
    }

    private float[] readDimension(ByteBuffer bb, int len) {
        float[] dim = new float[len];
        for (int i = 0; i < len; i++) {
            dim[i] = bb.getFloat();
        }
        return dim;
    }

    private ByteBuffer readFileToByteBuffer(String lutName) {
        ByteBuffer bb = null;
        File momoFile = new File(lutName);
        FileInputStream fis = null;
        try {
            fis = new FileInputStream(momoFile);
            byte[] buffer = new byte[(int) momoFile.length()];
            fis.read(buffer, 0, (int) momoFile.length());
            bb = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN);
        } catch (Exception ex) {
            System.err.println(ex.getMessage());
        } finally {
            if (fis != null) {
                try {
                    fis.close();
                } catch (IOException ex) {
                    System.err.println(ex.getMessage());
                }
            }
        }
        return bb;
    }

    private float[] readValues(ByteBuffer bb) {
        int len = nWvl * nAot * nHsf * nAzi * nSza * nVza * nParameter;
        float[] val = new float[len];
        for (int iWvl=0; iWvl<nWvl; iWvl++){
            for (int iAot = 0; iAot < nAot; iAot++) {
                for (int iHsf = nHsf-1; iHsf >= 0; iHsf--) {
                    for (int iAzi = nAzi-1; iAzi >= 0; iAzi--) {
                        for (int iSza = 0; iSza < nSza; iSza++) {
                            for (int iVza = 0; iVza < nVza; iVza++) {
                                for (int iPar = 0; iPar < nParameter; iPar++) {
                                    int pos = calcPosition(new int[]{iHsf, iVza, iSza, iAzi, iWvl, iAot, iPar},
                                                           new int[]{nHsf, nVza, nSza, nAzi, nWvl, nAot, nParameter});
                                    val[pos] = bb.getFloat();
                                }
                            }
                        }
                    }
                }
            }
        }
        return val;
    }

    private float[] readO3corr() {
        final InputStream inputStream = MomoLut.class.getResourceAsStream("o3Correction.asc");
        BufferedReader reader = null;
        reader = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        float[] o3cwvl = new float[15+4+4];
        float[] o3c = new float[15+4+4];
        try {
            int i = 0;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!(line.isEmpty() || line.startsWith("#") || line.startsWith("*"))) {
                    String[] stmp = line.split("[ \t]+");
                    o3cwvl[i] = Float.valueOf(stmp[1]);
                    o3c[i] = Float.valueOf(stmp[2]);
                    i++;
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(MomoLut.class.getName()).log(Level.SEVERE, null, ex);
            throw new OperatorException(ex.getMessage(), ex.getCause());
        }
        
        return getO3CorrOfBands(o3c, o3cwvl, wvl);
    }

    private float[] getO3CorrOfBands(float[] o3c, float[] o3cWvl, float[] lutWvl) {
        int jmin = 0;
        float[] lutO3c = new float[lutWvl.length];
        for (int i=0; i<lutWvl.length; i++){
            for (int j=0; j<o3cWvl.length; j++) {
                if (Math.abs(o3cWvl[j]-lutWvl[i]) < Math.abs(o3cWvl[jmin]-lutWvl[i])) jmin = j;
            }
            lutO3c[i] = o3c[jmin];
        }
        return lutO3c;
    }
}
