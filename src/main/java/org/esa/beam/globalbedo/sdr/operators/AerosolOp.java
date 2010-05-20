/*
 * Copyright (C) 2002-2007 by ?
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.globalbedo.sdr.operators;

import com.bc.ceres.core.ProgressMonitor;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.ProductUtils;

import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoLut;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.gpf.operators.standard.BandMathsOp;
import org.esa.beam.util.math.LookupTable;

// TODO - Rename this operator
// TODO - Adapt OperatorMetadata
// TODO - Adapt javadoc

/**
 * The sample operator implementation for an algorithm that outputs
 * all bands of the target product at once.
 */
@OperatorMetadata(alias = "AerosolOp",
                  description = "Computes <X> and <Y> from <Z>",
                  authors = "Me, myself and I",
                  version = "1.0",
                  copyright = "(C) 2010 by Brockmann Consult GmbH (beam@brockmann-consult.de)")
public class AerosolOp extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;

    //@Parameter
    private String SurfaceSpecName = "surface_reflectance_spec.asc";

    private String productName = "pname";
    private String productType = "ptype";

    private final String[] specBandNames = {"B0", "B2", "B3", "MIR"};
    private final String[] geomBandNames = {"SZA", "SAA", "VZA", "VAA"};
    private final String   surfPresName = "";

    private int rasterWidth;
    private int rasterHeight;
    private ArrayList<Band> specBandList;
    private ArrayList<Band> geometryBandList;
    private float[] specWvl;
    private int nSpecBands;
    private float[] soilSurfSpec;
    private float[] vegSurfSpec;
    private LookupTable lut;
    //private PointRetrieval retrieval;
    private Band validBand;

    /**
     * Default constructor. The graph processing framework
     * requires that an operator has a default constructor.
     */
    public AerosolOp() {
    }

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link org.esa.beam.framework.datamodel.Product}
     * annotated with the {@link org.esa.beam.framework.gpf.annotations.TargetProduct TargetProduct} annotation or
     * by calling {@link #setTargetProduct} method.</p>
     * <p>The framework calls this method after it has created this operator.
     * Any client code that must be performed before computation of tile data
     * should be placed here.</p>
     *
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          If an error occurs during operator initialisation.
     * @see #getTargetProduct()
     */
    @Override
    public void initialize() throws OperatorException {

        specBandList = getBandList(specBandNames);
        nSpecBands = specBandList.size();
        specWvl = getSpectralWvl(specBandList);

        geometryBandList = getBandList(geomBandNames);

        rasterHeight = sourceProduct.getSceneRasterHeight();
        rasterWidth  = sourceProduct.getSceneRasterWidth();

        readSurfaceSpectra(SurfaceSpecName);

        String lutName = "e:/model_data/momo/LUTs_Swansea/MERIS/MERIS_LUT_MOMO_ContinentalI_80_SU.bin";
        lut = new MomoLut(lutName).getLookupTable();
        

        final BandMathsOp validBandOp = BandMathsOp.createBooleanExpressionBand("(SM.B0_OK && SM.B2_OK && SM.B3_OK && SM.MIR_OK && SM.LAND && (MIR > 0.045))", sourceProduct);
        validBand = validBandOp.getTargetProduct().getBandAt(0);

        createTargetProduct();
    }

    /**
     * Called by the framework in order to compute the stack of tiles for the given target bands.
     * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
     *
     * @param targetTiles     The current tiles to be computed for each target band.
     * @param targetRectangle The area in pixel coordinates to be computed (same for all rasters in
     *                        <code>targetRasters</code>).
     * @param pm              A progress monitor which should be used to determine computation cancellation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          if an error occurs during computation of the target rasters.
     */
    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {
        //TODO create surface elevation map
        System.out.println("Tile:"+targetRectangle.x+"/"+rasterWidth+", "+targetRectangle.y+"/"+rasterHeight+"");
        Map<String, Tile> geomTiles = getTiles(geometryBandList, targetRectangle);
        Map<String, Tile> specTiles = getTiles(specBandList, targetRectangle);
        //Tile surfPressure = getSourceTile(sourceProduct.getBand(surfPresName), targetRectangle, ProgressMonitor.NULL);
        Tile validTile = getSourceTile(validBand, targetRectangle, ProgressMonitor.NULL);

        Tile aotTile = targetTiles.get(targetProduct.getBand("aot"));
        Tile errTile = targetTiles.get(targetProduct.getBand("aot_err"));
        Tile scaleVegTile = targetTiles.get(targetProduct.getBand("scaleVeg"));
        Tile scaleSoilTile = targetTiles.get(targetProduct.getBand("scaleSoil"));
        Tile[] sdrTiles = new Tile[nSpecBands];
        Tile[] modelTiles = new Tile[nSpecBands];
        for (int i=0; i<nSpecBands; i++) {
            sdrTiles[i] = targetTiles.get(targetProduct.getBand(String.format("sdr_%1d",i)));
            modelTiles[i] = targetTiles.get(targetProduct.getBand(String.format("model_%1d",i)));
        }
        sourceProduct.getBandAt(0).getValidMaskImage();

        PointRetrieval retrieval = new PointRetrieval(lut, soilSurfSpec, vegSurfSpec, specWvl);

        int x0 = (int) targetRectangle.getX();
        int y0 = (int) targetRectangle.getY();
        int width = (int) targetRectangle.getWidth() + x0 - 1;
        int height = (int) targetRectangle.getHeight() + y0 - 1;
        for (int iY=y0; iY <= height; iY++) {
            for (int iX=x0; iX <= width; iX++) {
                float aot = (float) aotTile.getRasterDataNode().getGeophysicalNoDataValue();
                float err = (float) errTile.getRasterDataNode().getGeophysicalNoDataValue();
                float scaleVeg = (float) scaleVegTile.getRasterDataNode().getGeophysicalNoDataValue();
                float scaleSoil = (float) scaleSoilTile.getRasterDataNode().getGeophysicalNoDataValue();
                float[] sdr = new float[nSpecBands];
                float[] model = new float[nSpecBands];
                for (int i=0; i<nSpecBands; i++) {
                    sdr[i] = (float) sdrTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                    model[i] = (float) modelTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                }

                if (validTile.getSampleBoolean(iX, iY)) {
                    float sza = (float) geomTiles.get(geomBandNames[0]).getSampleDouble(iX, iY);
                    float saa = (float) geomTiles.get(geomBandNames[1]).getSampleDouble(iX, iY);
                    float vza = (float) geomTiles.get(geomBandNames[2]).getSampleDouble(iX, iY);
                    float vaa = (float) geomTiles.get(geomBandNames[3]).getSampleDouble(iX, iY);
                    float surfPressure = 1013.25f;
                    float[] geometry = {sza, getRelativeAzi(saa, vaa), vza};
                    float[] toaReflec = getSpectrum(specTiles, iX, iY);

                    RetrievalResults retrievalResult = retrieval.runRetrieval(surfPressure, geometry, toaReflec);
                    aot = retrievalResult.optAOT;
                    err = retrievalResult.optErr;
                    scaleVeg = retrievalResult.scaleVeg;
                    scaleSoil = retrievalResult.scaleSoil;
                    sdr = retrievalResult.sdr;
                    for (int i=0; i<nSpecBands; i++) {
                        model[i] = retrievalResult.scaleSoil * soilSurfSpec[i];
                        model[i] += retrievalResult.scaleVeg * vegSurfSpec[i];
                    }

                }
                aotTile.setSample(iX, iY, aot);
                errTile.setSample(iX, iY, err);
                scaleVegTile.setSample(iX, iY, scaleVeg);
                scaleSoilTile.setSample(iX, iY, scaleSoil);
                for (int i=0; i<nSpecBands; i++) {
                    sdrTiles[i].setSample(iX, iY, sdr[i]);
                    modelTiles[i].setSample(iX, iY, model[i]);
                }
            }
        }

    }

    private void createTargetProduct() {

        targetProduct = new Product(productName, productType, rasterWidth, rasterHeight);

        ProductUtils.copyMetadata(sourceProduct, targetProduct);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        ProductUtils.copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyFlagBands(sourceProduct, targetProduct);

        //AerosolHelpers.addAerosolFlagBand(targetProduct, downscaledRasterWidth, downscaledRasterHeight);

        createTargetProductBands();

        targetProduct.setPreferredTileSize(128,128);
        setTargetProduct(targetProduct);

    }

    private void createTargetProductBands() {

        for (Band srcBand : sourceProduct.getBands()) {
            if (!srcBand.isFlagBand()) {
                ProductUtils.copyBand(srcBand.getName(), sourceProduct, targetProduct);
            }
            Band tarBand = targetProduct.getBand(srcBand.getName());
            tarBand.setSourceImage(srcBand.getSourceImage());
        }

        Band targetBand = new Band("aot", ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
        targetBand.setDescription("best fitting aot Band");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setValidPixelExpression(targetBand.getName() + ">= 0 AND " + targetBand.getName() + "<= 2");
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("aot_err", ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("scaleVeg", ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("scaleSoil", ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        for (int i=0; i<nSpecBands; i++) {
            targetBand = new Band(String.format("sdr_%1d",i), ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setUnit("dl");
            targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
            targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
            targetProduct.addBand(targetBand);
        }
        for (int i=0; i<nSpecBands; i++) {
            targetBand = new Band(String.format("model_%1d",i), ProductData.TYPE_FLOAT32, rasterWidth, rasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setUnit("dl");
            targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
            targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
            targetProduct.addBand(targetBand);
        }
    }

    private ArrayList<Band> getBandList(String[] bandNames) {
        ArrayList<Band> list = new ArrayList<Band>(bandNames.length);
        for (int i=0; i<bandNames.length; i++) {
            list.add(sourceProduct.getBand(bandNames[i]));
        }
        if (list.get(0).getSpectralWavelength() > 0) {
            Comparator<Band> byWavelength = new WavelengthComparator();
            Collections.sort(list, byWavelength);
        }

        return list;
    }

    private float getRelativeAzi(float saa, float vaa) {
        float relAzi = Math.abs(saa - vaa);
        relAzi = (relAzi > 180.0f) ? 180 - (360 - relAzi) : 180 - relAzi;
        return relAzi;
    }

    private float[] getSpectralWvl(ArrayList<Band> specBandList) {
        float[] wvl = new float[specBandList.size()];
        for (int i=0; i<specBandList.size(); i++) {
            wvl[i] = specBandList.get(i).getSpectralWavelength();
        }
        return wvl;
    }

    private float[] getSpectrum(Map<String, Tile> specTiles, int iX, int iY) {
        float[] spec = new float[specTiles.size()];
        for (int i=0; i<specTiles.size(); i++) {
            spec[i] = specTiles.get(specBandNames[i]).getSampleFloat(iX, iY);
        }
        return spec;
    }

    private Map<String, Tile> getTiles(ArrayList<Band> bandList, Rectangle rec) {
        Map<String, Tile> tileMap = new HashMap<String, Tile>(bandList.size());
        for (Band b : bandList) {
            tileMap.put(b.getName(), getSourceTile(b, rec, ProgressMonitor.NULL));
        }
        return tileMap;
    }

    private void readSurfaceSpectra(String fname) {
        final InputStream inputStream = AerosolOp.class.getResourceAsStream(fname);
        BufferedReader reader = null;
        reader = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        float[] fullWvl  = new float[10000];
        float[] fullSoil = new float[10000];
        float[] fullVeg  = new float[10000];
        int nWvl = 0;
        try {
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!(line.isEmpty() || line.startsWith("#") || line.startsWith("*"))) {
                    String[] stmp = line.split("[ \t]+");
                    fullWvl[nWvl] = Float.valueOf(stmp[0]);
                    if (fullWvl[nWvl] < 100) fullWvl[nWvl] *= 1000; // conversion from um to nm
                    fullSoil[nWvl] = Float.valueOf(stmp[1]);
                    fullVeg[nWvl] = Float.valueOf(stmp[2]);
                    nWvl++;
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(AerosolOp.class.getName()).log(Level.SEVERE, null, ex);
            throw new OperatorException(ex.getMessage(), ex.getCause());
        }

        soilSurfSpec = new float[nSpecBands];
        vegSurfSpec  = new float[nSpecBands];
        int j = 0;
        for (int i=0; i<nSpecBands; i++) {
            float wvl = specBandList.get(i).getSpectralWavelength();
            float width = specBandList.get(i).getSpectralBandwidth();
            int count = 0;
            while (j < nWvl && fullWvl[j] < wvl-width/2) j++;
            if (j == nWvl) throw new OperatorException("wavelength not found reading surface spectra");
            while (fullWvl[j] < wvl+width/2) {
                soilSurfSpec[i] += fullSoil[j];
                vegSurfSpec[i]  += fullVeg[j];
                count++; j++;
            }
            if (j == nWvl) throw new OperatorException("wavelength window exceeds surface spectra range");
            if (count > 0) {
                soilSurfSpec[i] /= count;
                vegSurfSpec[i]  /= count;
            }
        }
    }


    /**
     * The SPI is used to register this operator in the graph processing framework
     * via the SPI configuration file
     * {@code META-INF/services/org.esa.beam.framework.gpf.OperatorSpi}.
     * This class may also serve as a factory for new operator instances.
     *
     * @see OperatorSpi#createOperator()
     * @see OperatorSpi#createOperator(java.util.Map, java.util.Map)
     */
    public static class Spi extends OperatorSpi {
        public Spi() {
            super(AerosolOp.class);
        }
    }
}
