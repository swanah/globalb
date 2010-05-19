/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.lutUtils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import org.esa.beam.util.math.LookupTable;

/**
 *
 * @author akheckel
 */
public class MomoLut {
    private final int nWvl = 15;
    private final int nParameter = 4;
    private int nVza;
    private int nSza;
    private int nAzi;
    private int nHsf;
    private int nAot;
    private float[] vza;
    private float[] sza;
    private float[] azi;
    private float[] hsf;
    private float[] aot;
    private float[] wvl;
    private float[] values;
    private float solIrrad;


    public MomoLut(String lutName) {
        File momoFile = new File(lutName);
        FileInputStream fis = null;
        try {
            fis = new FileInputStream(momoFile);
            byte[] buffer = new byte[(int) momoFile.length()];
            fis.read(buffer, 0, (int) momoFile.length());

            ByteBuffer bb = ByteBuffer.wrap(buffer).order(ByteOrder.LITTLE_ENDIAN);

            // read LUT
            vza = readDimension(bb); nVza = vza.length;
            sza = readDimension(bb); nSza = sza.length;
            azi = readDimension(bb); nAzi = azi.length;
            hsf = readDimension(bb); nHsf = hsf.length;
            aot = readDimension(bb); nAot = aot.length;
            values = readValues(bb);
            wvl = readDimension(bb, nWvl);
            solIrrad = bb.getFloat();
            
            // invert order to store stirctly increasing dimension in lut
            // swapping the actual values in the lut is taken care off in readValues()
            for (int i=0; i<nHsf/2; i++){
                float swap = hsf[nHsf-1-i];
                hsf[nHsf-1-i] = hsf[i];
                hsf[i] = swap;
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            if (fis != null) {
                try {
                    fis.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }
        }

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

    private float[] readValues(ByteBuffer bb) {
        int len = nWvl * aot.length * hsf.length * azi.length * sza.length * vza.length * nParameter;
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


    public int calcPosition(int[] indices, int[] sizes) {
        //(((((iHsf*nVza + iVza)* nSza + iSza)* nAzi + iAzi)* nWvl + iWvl)* nAot + iAot)* nParameter + iPar
        int pos = 0;
        for (int i=0; i<sizes.length; i++) pos = (pos * sizes[i] + indices[i]);
        return pos;
    }

    public float[] getAot() {
        return aot;
    }

    public float[] getAzi() {
        return azi;
    }

    public float[][] getDimensions() {
        return new float[][]{hsf, vza, sza, azi, wvl, aot, new float[]{0,1,2,3}};
    }

    public float[] getHsf() {
        return hsf;
    }

    public LookupTable getLookupTable() {
        return new LookupTable(values, getDimensions());
    }

    public float getSolIrrad() {
        return solIrrad;
    }

    public float[] getSza() {
        return sza;
    }

    public float[] getValues() {
        return values;
    }

    public float[] getVza() {
        return vza;
    }

    public float[] getWvl() {
        return wvl;
    }

}
