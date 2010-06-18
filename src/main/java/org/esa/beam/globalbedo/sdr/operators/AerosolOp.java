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
import javax.media.jai.BorderExtender;
import org.esa.beam.akh.createelevationband.CreateElevationBandOp;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoLut;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoSynLut;
import org.esa.beam.gpf.operators.standard.BandMathsOp;
import org.esa.beam.util.math.LookupTable;

// TODO - Adapt javadoc

/**
 * Operator producing AOT for GlobAlbedo
 */
@OperatorMetadata(alias = "AerosolOp",
                  description = "Computes aerosol optical thickness",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class AerosolOp extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;
    @Parameter(defaultValue="false")
    private boolean retrieveAOT = false;
    //@Parameter
    private String SurfaceSpecName = "surface_reflectance_spec.asc";

    private String productName = "pname";
    private String productType = "ptype";

    private String instrument;

    private String[] specBandNames;
    private String[] geomBandNames;
    private String surfPresName;
    private String ozoneName;

    private int srcRasterWidth;
    private int srcRasterHeight;
    private ArrayList<Band> specBandList;
    private ArrayList<RasterDataNode> geometryBandList;
    private float[] specWvl;
    private int nSpecBands;
    private float[] soilSurfSpec;
    private float[] vegSurfSpec;
    private float[] specWeights;
    private LookupTable lut;
    private MomoLut momo;
    //private PointRetrieval retrieval;
    private Band validBand;
    private VirtualBand surfPresBand;
    private LookupTable[] synLut;
    private float[] lutAlb;
    private boolean useSynLut = false;
    private BorderExtender borderExt;
    private Rectangle pixelWindow;

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

        srcRasterHeight = sourceProduct.getSceneRasterHeight();
        srcRasterWidth  = sourceProduct.getSceneRasterWidth();

        instrument = getInstrument();

        specBandList = getSpecBandList(instrument);
        nSpecBands = specBandList.size();
        specWvl = getSpectralWvl(specBandList);
        specBandNames = InstrumentConsts.getInstance().getSpecBandNames(instrument);

        geometryBandList = getGeomBandList(instrument);
        geomBandNames = InstrumentConsts.getInstance().getGeomBandNames(instrument);

        surfPresName = InstrumentConsts.getInstance().getSurfPressureName(instrument);
        ozoneName    = InstrumentConsts.getInstance().getOzoneName(instrument);
        if (instrument.equals("VGT") && !sourceProduct.containsBand(surfPresName)) {
            sourceProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(CreateElevationBandOp.class), GPF.NO_PARAMS, sourceProduct);
            String presExpr = "(1013.25 * exp(-elev/8400))";
            surfPresBand = new VirtualBand(surfPresName, ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight, presExpr);
            surfPresBand.setDescription("estimated sea level pressure (p0=1013.25hPa, hScale=8.4km)");
            surfPresBand.setNoDataValue(0);
            surfPresBand.setNoDataValueUsed(true);
            surfPresBand.setUnit("hPa");
            sourceProduct.addBand(surfPresBand);
        }

        String validExpression = InstrumentConsts.getInstance().getValidExpression(instrument);
        final BandMathsOp validBandOp = BandMathsOp.createBooleanExpressionBand(validExpression, sourceProduct);
        validBand = validBandOp.getTargetProduct().getBandAt(0);

        readSurfaceSpectra(SurfaceSpecName);
        specWeights = InstrumentConsts.getInstance().getFitWeights(instrument);
        if (useSynLut){
            lut = null;
            momo = null;
            MomoSynLut momoSynLut = new MomoSynLut("e:/model_data/momo/bin", 8, "MERIS", specWvl);
            lutAlb = momoSynLut.getAlbDim();
            synLut = momoSynLut.getLut();
        }
        else {
            String lutName = InstrumentConsts.getInstance().getLutName(instrument);
            int nLutBands = InstrumentConsts.getInstance().getnLutBands(instrument);
            momo = new MomoLut(lutName, nLutBands);
            lut = momo.getLookupTable();
            lutAlb = null;
            synLut = null;
        }

        createTargetProduct();
        borderExt = BorderExtender.createInstance(BorderExtender.BORDER_COPY);
        pixelWindow = new Rectangle(0, 0, 2, 2);
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
        pm.beginTask("begin Aerosol retrieval", srcRasterHeight);
        System.out.println("Tile:"+targetRectangle.x+"/"+srcRasterWidth+", "+targetRectangle.y+"/"+srcRasterHeight);

        Rectangle srcRec = addBorder(targetRectangle, pixelWindow);

        Map<String, Tile> inputTiles = getTiles(geometryBandList, srcRec);
        inputTiles.putAll(getTiles(specBandList, srcRec));
        Tile surfPressure = getSourceTile(sourceProduct.getBand(surfPresName), srcRec, borderExt, ProgressMonitor.NULL);
        inputTiles.put(surfPresName, surfPressure);
        Tile ozoneTile = getSourceTile(sourceProduct.getBand(ozoneName), srcRec, borderExt, ProgressMonitor.NULL);
        inputTiles.put(ozoneName, ozoneTile);
        Tile validTile = getSourceTile(validBand, srcRec, borderExt, ProgressMonitor.NULL);
        inputTiles.put(validBand.getName(), validTile);

        Tile srcAotTile = null;
        if (instrument.equals("VGT") && !retrieveAOT){
            srcAotTile =getSourceTile(sourceProduct.getBand("AG"), targetRectangle, ProgressMonitor.NULL);
        }

        Tile aotTile = targetTiles.get(targetProduct.getBand("aot"));
        Tile errTile = targetTiles.get(targetProduct.getBand("aot_err"));
        Tile scaleVegTile = targetTiles.get(targetProduct.getBand("scaleVeg"));
        Tile scaleSoilTile = targetTiles.get(targetProduct.getBand("scaleSoil"));
        Tile addConstTile = targetTiles.get(targetProduct.getBand("addConst"));
        Tile[] sdrTiles = new Tile[nSpecBands];
        Tile[] modelTiles = new Tile[nSpecBands];
        for (int i=0; i<nSpecBands; i++) {
            sdrTiles[i] = targetTiles.get(targetProduct.getBand(specBandNames[i]+"_sdr"));
            modelTiles[i] = targetTiles.get(targetProduct.getBand(specBandNames[i]+"_model"));
        }

        int x0 = (int) targetRectangle.getX();
        int y0 = (int) targetRectangle.getY();
        int width = (int) targetRectangle.getWidth() + x0 - 1;
        int height = (int) targetRectangle.getHeight() + y0 - 1;
        for (int iY=y0; iY <= height; iY++) {
            //System.out.printf(outS+" %f %% \r", (100.0*iY/height));
            for (int iX=x0; iX <= width; iX++) {
                float aot = (float) aotTile.getRasterDataNode().getGeophysicalNoDataValue();
                float err = (float) errTile.getRasterDataNode().getGeophysicalNoDataValue();
                float scaleVeg = (float) scaleVegTile.getRasterDataNode().getGeophysicalNoDataValue();
                float scaleSoil = (float) scaleSoilTile.getRasterDataNode().getGeophysicalNoDataValue();
                float addConst = (float) addConstTile.getRasterDataNode().getGeophysicalNoDataValue();
                float[] sdr = new float[nSpecBands];
                float[] model = new float[nSpecBands];
                for (int i=0; i<nSpecBands; i++) {
                    sdr[i] = (float) sdrTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                    model[i] = (float) modelTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                }

                if (validTile.getSampleBoolean(iX, iY)) {
                    //float[] saa = (float) geomTiles.get(geomBandNames[0]).getSampleDouble(iX, iY);
                    //float[] saa = (float) geomTiles.get(geomBandNames[1]).getSampleDouble(iX, iY);
                    //float[] vza = (float) geomTiles.get(geomBandNames[2]).getSampleDouble(iX, iY);
                    //float[] vaa = (float) geomTiles.get(geomBandNames[3]).getSampleDouble(iX, iY);
                    //float surfPressure = 1013.25f;
                    //float[] toaReflec = getSpectrum(specTiles, iX, iY);
                    //PixelGeometry geom = new PixelGeometry(sza, saa, vza, vaa);
                    //InputPixelData inPixData = new InputPixelData(geom, surfPressure, specWvl, toaReflec);
                    SourcePixelField inPixData = readPixelData(inputTiles, iX, iY, pixelWindow);
                    emodSpecTau brentFitFct;
                    if (useSynLut){
                        brentFitFct = new emodSpecTau(inPixData, synLut, lutAlb, soilSurfSpec, vegSurfSpec, specWeights);
                    }
                    else {
                        brentFitFct = new emodSpecTau(inPixData, momo, soilSurfSpec, vegSurfSpec, specWeights);
                    }

                    RetrievalResults retrievalResult;
                    if (retrieveAOT) retrievalResult = new PointRetrieval().runRetrieval(brentFitFct);
                    else if (instrument.equals("VGT")) {
                        retrievalResult = new PointRetrieval().retrieveSDR(brentFitFct, srcAotTile.getSampleFloat(iX,iY));
                    }
                    else {
                        retrievalResult = new PointRetrieval().retrieveSDR(brentFitFct, 0.1);
                    }

                    aot = retrievalResult.optAOT;
                    err = retrievalResult.optErr;
                    scaleVeg = (float) retrievalResult.pAtMin[0];
                    scaleSoil = (float) retrievalResult.pAtMin[1];
                    addConst = (float) retrievalResult.pAtMin[2];
                    sdr = retrievalResult.sdr;
                    model = retrievalResult.modelSpec;
                }
                aotTile.setSample(iX, iY, aot);
                errTile.setSample(iX, iY, err);
                scaleVegTile.setSample(iX, iY, scaleVeg);
                scaleSoilTile.setSample(iX, iY, scaleSoil);
                addConstTile.setSample(iX, iY, addConst);
                for (int i=0; i<nSpecBands; i++) {
                    sdrTiles[i].setSample(iX, iY, sdr[i]);
                    modelTiles[i].setSample(iX, iY, model[i]);
                }
                if(pm.isCanceled()) return;
            }
            pm.worked(1);

        }
        pm.done();

    }

    private void createTargetProduct() {

        targetProduct = new Product(productName, productType, srcRasterWidth, srcRasterHeight);

        ProductUtils.copyMetadata(sourceProduct, targetProduct);
        ProductUtils.copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
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

        Band targetBand = new Band("aot", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("best fitting aot Band");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setValidPixelExpression(targetBand.getName() + ">= 0 AND " + targetBand.getName() + "<= 2");
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("aot_err", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("scaleVeg", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("scaleSoil", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("addConst", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        for (int i=0; i<nSpecBands; i++) {
            targetBand = new Band(specBandNames[i]+"_sdr", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setUnit("dl");
            targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
            targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
            targetProduct.addBand(targetBand);
        }
        for (int i=0; i<nSpecBands; i++) {
            targetBand = new Band(specBandNames[i]+"_model", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setUnit("dl");
            targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
            targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
            targetProduct.addBand(targetBand);
        }
    }

    private ArrayList<RasterDataNode> getGeomBandList(String instrument) {
        String[] bandNames = InstrumentConsts.getInstance().getGeomBandNames(instrument);
        ArrayList<RasterDataNode> list = new ArrayList<RasterDataNode>();
        for (int i=0; i<bandNames.length; i++) {
            list.add(sourceProduct.getRasterDataNode(bandNames[i]));
        }
        return list;
    }

    private ArrayList<Band> getSpecBandList(String instrument) {
        String[] bandNames = InstrumentConsts.getInstance().getSpecBandNames(instrument);
        ArrayList<Band> list = new ArrayList<Band>();
        for (int i=0; i<bandNames.length; i++) {
            list.add(sourceProduct.getBand(bandNames[i]));
        }
        if (list.get(0).getSpectralWavelength() > 0) {
            Comparator<Band> byWavelength = new WavelengthComparator();
            Collections.sort(list, byWavelength);
        }
        return list;
    }

    private String getInstrument() {
        String inst = null;
        String[] supportedInstruments = InstrumentConsts.getInstance().getSupportedInstruments();
        for (String suppInstr : supportedInstruments) {
            String[] specBands = InstrumentConsts.getInstance().getSpecBandNames(suppInstr);
            if (specBands.length > 0 && sourceProduct.containsBand(specBands[0])) {
                inst = suppInstr;
                break;
            }
        }
        if (inst.equals(null)) throw new OperatorException("Product not supported.");
        return inst;
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

    private <T extends RasterDataNode> Map<String, Tile> getTiles(ArrayList<T> bandList, Rectangle rec) {
        Map<String, Tile> tileMap = new HashMap<String, Tile>(bandList.size());
        for (T b : bandList) {
            tileMap.put(b.getName(), getSourceTile(b, rec, borderExt, ProgressMonitor.NULL));
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

    private Rectangle addBorder(Rectangle rec, Rectangle borderRec) {
        return new Rectangle(rec.x + borderRec.x,
                             rec.y + borderRec.y,
                             rec.width + borderRec.width,
                             rec.height + borderRec.height);
    }

    private SourcePixelField readPixelData(Map<String, Tile> inputTiles, int iX, int iY, Rectangle window) {
        InputPixelData[] inPixel = new InputPixelData[window.width * window.height];
        boolean[] valid = new boolean[window.width * window.height];

        for (int j = 0; j < window.height; j++){
            for (int i = 0; i < window.width; i++){
                int tileX = i + iX + window.x;
                int tileY = j + iY + window.y;
                double noData = inputTiles.get(geomBandNames[0]).getRasterDataNode().getGeophysicalNoDataValue();
                float sza = inputTiles.get(geomBandNames[0]).getSampleFloat(tileX, tileY);
                float saa = inputTiles.get(geomBandNames[1]).getSampleFloat(tileX, tileY);
                float vza = inputTiles.get(geomBandNames[2]).getSampleFloat(tileX, tileY);
                float vaa = inputTiles.get(geomBandNames[3]).getSampleFloat(tileX, tileY);
                PixelGeometry geom = new PixelGeometry(sza, saa, vza, vaa);
                valid[i+j*window.width] = (sza != (float)noData);

                noData = inputTiles.get(surfPresName).getRasterDataNode().getGeophysicalNoDataValue();
                float pres = inputTiles.get(surfPresName).getSampleFloat(tileX, tileY);
                valid[i+j*window.width] = valid[i+j*window.width] && (pres != noData);

                noData = inputTiles.get(surfPresName).getRasterDataNode().getGeophysicalNoDataValue();
                float o3col = inputTiles.get(ozoneName).getSampleFloat(tileX, tileY) * 1000; // O3 in [DU]
                valid[i+j*window.width] = valid[i+j*window.width] && (o3col != noData);

                float[] spec = new float[specWvl.length];
                for (int k=0; k<specWvl.length; k++) {
                    Tile t = inputTiles.get(specBandNames[k]);
                    spec[k] = t.getSampleFloat(tileX, tileY);
                    valid[i+j*window.width] = valid[i+j*window.width] && (spec[k] != t.getRasterDataNode().getGeophysicalNoDataValue());
                }
                inPixel[i+j*window.width] = new InputPixelData(geom, pres, o3col, specWvl, spec);
            }
        }
        
        return new SourcePixelField(inPixel, valid, (-window.x)+(-window.y)*window.width);
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
