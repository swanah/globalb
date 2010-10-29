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
import java.util.List;
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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import javax.media.jai.BorderExtender;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoLut;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.globalbedo.sdr.lutUtils.Aardvarc4DLut;
import org.esa.beam.globalbedo.sdr.lutUtils.MomoSynLut;
import org.esa.beam.gpf.operators.standard.BandMathsOp;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.math.LookupTable;
import org.esa.beam.util.math.RsMathUtils;

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
    @Parameter(defaultValue="false")
    private boolean saveToaBands = false;
    @Parameter(defaultValue="false")
    private boolean saveSdrBands = false;
    @Parameter(defaultValue="false")
    private boolean saveModelBands = false;

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
    private ArrayList<TiePointGrid> geometryBandList;
    private float[] specWvl;
    private int nSpecBands;
    private float[] soilSurfSpec;
    private float[] vegSurfSpec;
    private double[] specWeights;
    private MomoLut momo;
    private Band validBand;
    private VirtualBand surfPresBand;
    private LookupTable[] synLut;
    private float[] lutAlb;
    private boolean useSynLut = false;
    private Aardvarc4DLut aardvarcLut;
    private boolean useAardvarcLut = false;
    private BorderExtender borderExt;
    private String validExpression;
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

        instrument = InstrumentConsts.getInstance().getInstrument(sourceProduct);

        specBandList = getSpecBandList(instrument);
        specWvl = getSpectralWvl(specBandList);
        nSpecBands = specBandList.size();
        specBandNames = InstrumentConsts.getInstance().getSpecBandNames(instrument);

        geometryBandList = getGeomBandList(instrument);
        geomBandNames = InstrumentConsts.getInstance().getGeomBandNames(instrument);

        surfPresName = InstrumentConsts.getInstance().getSurfPressureName(instrument);
        ozoneName    = InstrumentConsts.getInstance().getOzoneName(instrument);
        /*
        if (!sourceProduct.containsBand(surfPresName)){
            String presExpr = null;
            if (!(sourceProduct.containsBand("elevation") || sourceProduct.containsTiePointGrid("elevation"))){
                sourceProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(CreateElevationBandOp.class), GPF.NO_PARAMS, sourceProduct);
            }
            presExpr = "(1013.25 * exp(-elevation/8400))";
            surfPresBand = new VirtualBand(surfPresName, ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight, presExpr);
            surfPresBand.setDescription("estimated sea level pressure (p0=1013.25hPa, hScale=8.4km)");
            surfPresBand.setNoDataValue(0);
            surfPresBand.setNoDataValueUsed(true);
            surfPresBand.setUnit("hPa");
            sourceProduct.addBand(surfPresBand);
        }
         */
        if (!sourceProduct.containsBand(ozoneName)){
            String ozoneExpr = "0.350";
            VirtualBand ozoneBand = new VirtualBand(ozoneName, ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight, ozoneExpr);
            ozoneBand.setDescription("constant ozone band 350DU");
            ozoneBand.setNoDataValue(0);
            ozoneBand.setNoDataValueUsed(true);
            ozoneBand.setUnit("atm.cm");
            sourceProduct.addBand(ozoneBand);
        }

        validExpression = InstrumentConsts.getInstance().getValidExpression(instrument);
        final BandMathsOp validBandOp = BandMathsOp.createBooleanExpressionBand(validExpression, sourceProduct);
        validBand = validBandOp.getTargetProduct().getBandAt(0);

        readSurfaceSpectra(SurfaceSpecName);
        specWeights = InstrumentConsts.getInstance().getSpectralFitWeights(instrument);
        if (useSynLut){
            momo = null;
            aardvarcLut = null;
            MomoSynLut momoSynLut = new MomoSynLut("e:/model_data/momo/bin", 8, "MERIS", specWvl);
            lutAlb = momoSynLut.getAlbDim();
            synLut = momoSynLut.getLut();
        }
        else {
            String lutName = InstrumentConsts.getInstance().getLutName(instrument);
            int nLutBands = InstrumentConsts.getInstance().getnLutBands(instrument);
            momo = new MomoLut(lutName, nLutBands);
            lutAlb = null;
            synLut = null;
            if (useAardvarcLut){
                Guardian.assertEquals("instrument", instrument, "AATSR");
                lutName = "e:/projects/Synergy/wgrey/AATSR/src/aardvarc/aardvarc_v1/LUT6S/LUT6s_1_4D";
                aardvarcLut = new Aardvarc4DLut(lutName);
            }
        }

        createTargetProduct();
        borderExt = BorderExtender.createInstance(BorderExtender.BORDER_COPY);
        pixelWindow = new Rectangle(0, 0, 1,1);
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
        Tile scaleVegTile = null;
        Tile scaleSoilTile = null;
        Tile addConstTile = null;
        Tile[] sdrTiles = null;
        Tile[] modelTiles = null;
        if (saveModelBands) {
            scaleVegTile = targetTiles.get(targetProduct.getBand("scaleVeg"));
            scaleSoilTile = targetTiles.get(targetProduct.getBand("scaleSoil"));
            addConstTile = targetTiles.get(targetProduct.getBand("addConst"));
            modelTiles = new Tile[nSpecBands];
            for (int i=0; i<nSpecBands; i++) {
                modelTiles[i] = targetTiles.get(targetProduct.getBand(specBandNames[i]+"_model"));
            }
        }
        if (saveSdrBands) {
            sdrTiles = new Tile[nSpecBands];
            for (int i=0; i<nSpecBands; i++) {
                sdrTiles[i] = targetTiles.get(targetProduct.getBand(specBandNames[i]+"_sdr"));
            }
        }

        int x0 = (int) targetRectangle.getX();
        int y0 = (int) targetRectangle.getY();
        int width = (int) targetRectangle.getWidth() + x0 - 1;
        int height = (int) targetRectangle.getHeight() + y0 - 1;
        for (int iY=y0; iY <= height; iY++) {
            //System.out.printf(outS+" %f %% \r", (100.0*iY/height));
            for (int iX=x0; iX <= width; iX++) {
                //System.err.println(iX+" / "+iY);
                float aot = (float) aotTile.getRasterDataNode().getGeophysicalNoDataValue();
                float err = (float) errTile.getRasterDataNode().getGeophysicalNoDataValue();
                float scaleVeg = 0;
                float scaleSoil = 0;
                float addConst = 0;
                float[] model = null;
                if (saveModelBands){
                    scaleVeg = (float) scaleVegTile.getRasterDataNode().getGeophysicalNoDataValue();
                    scaleSoil = (float) scaleSoilTile.getRasterDataNode().getGeophysicalNoDataValue();
                    addConst = (float) addConstTile.getRasterDataNode().getGeophysicalNoDataValue();
                    model = new float[nSpecBands];
                    for (int i=0; i<nSpecBands; i++) {
                        model[i] = (float) modelTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                    }
                }
                float[] sdr = null;
                if (saveSdrBands){
                    sdr = new float[nSpecBands];
                    for (int i=0; i<nSpecBands; i++) {
                        sdr[i] = (float) sdrTiles[i].getRasterDataNode().getGeophysicalNoDataValue();
                    }
                }

                if (validTile.getSampleBoolean(iX, iY)) {
                    SourcePixelField srcPixelF = readPixelData(inputTiles, iX, iY, pixelWindow);
                    BrentFitFunction brentFitFunction = null;
                    emodSpecTau specBrentFct = null;
                    emodAngTau  angBrentFct = null;
                    emodSynTau  synBrentFct = null;
                    if (useSynLut){
                        specBrentFct = new emodSpecTau(srcPixelF, synLut, lutAlb, soilSurfSpec, vegSurfSpec, specWeights);
                        angBrentFct  = new emodAngTau(srcPixelF, synLut, lutAlb, specWeights);
                        synBrentFct  = new emodSynTau(srcPixelF, synLut, lutAlb, soilSurfSpec, vegSurfSpec, specWeights);

                    }
                    else if (useAardvarcLut){
                        angBrentFct  = new emodAngTau(srcPixelF, aardvarcLut, specWeights);
                        //brentFitFunction = new BrentFitFunction(BrentFitFunction.ANGULAR_MODEL, srcPixelF, aardvarcLut, specWeights);
                    }
                    else {
                        specBrentFct = new emodSpecTau(srcPixelF, momo, soilSurfSpec, vegSurfSpec, specWeights);
                        angBrentFct  = new emodAngTau(srcPixelF, momo, specWeights);
                        //synBrentFct  = new emodSynTau(srcPixelF, momo, soilSurfSpec, vegSurfSpec, specWeights);
                        //brentFitFunction = new BrentFitFunction(BrentFitFunction.ANGULAR_MODEL, srcPixelF, momo, specWeights);
                    }

                    RetrievalResults retrievalResult = null;
                    if (retrieveAOT) {
                        if (instrument.equals("AATSR")){
                            /*
                            if (iX==0 && iY==0){
                                System.out.printf("%f  %f  \n%f  %f\n", 90.0-inputTiles.get("sun_elev_nadir").getSampleFloat(0, 0),
                                        90.0-inputTiles.get("view_elev_nadir").getSampleFloat(0, 0),
                                        inputTiles.get("sun_azimuth_nadir").getSampleFloat(0, 0),
                                        inputTiles.get("view_azimuth_nadir").getSampleFloat(0, 0));
                                System.out.printf("%f  %f  %f\n", srcPixelF.inPixArr[srcPixelF.pixelIndex].geom.sza,
                                        srcPixelF.inPixArr[srcPixelF.pixelIndex].geom.vza,
                                        srcPixelF.inPixArr[srcPixelF.pixelIndex].geom.razi);
                                System.out.printf("fmin: %f\n", angBrentFct.f(0.2));
                                //float[][] sdrDiff = aardvarcLut.getSdrDiff(inPixData.inPixArr[inPixData.pixelIndex], 0.2f);
                                //for(int i=0; i<4; i++) System.out.printf("%f  %f  %f\n",specWvl[i],inPixData.inPixArr[inPixData.pixelIndex].toaReflec[i],sdrDiff[0][i]);
                            }
                            */
                            retrievalResult = new PointRetrieval(angBrentFct).runRetrieval(instrument);
                            //retrievalResult = new PointRetrieval(synBrentFct).runRetrieval(instrument);
                        }
                        else {
                            retrievalResult = new PointRetrieval(specBrentFct).runRetrieval(instrument);
                        }
                    }
                    /*
                    else if (instrument.equals("VGT")) {
                        retrievalResult = new PointRetrieval(specBrentFct).retrieveSDR(srcAotTile.getSampleFloat(iX,iY));
                    }
                     *
                     */
                    else {
                        if (instrument.equals("AATSR")){
                            retrievalResult = new PointRetrieval(angBrentFct).retrieveSDR(0.1);
                        }
                        else {
                            retrievalResult = new PointRetrieval(specBrentFct).retrieveSDR(0.1);
                        }
                    }
                    if (!retrievalResult.retrievalFailed){
                        aot = retrievalResult.optAOT;
                        err = retrievalResult.retrievalErr;
/*
                        if(instrument.equals("VGT")){
                            scaleVeg = (float) retrievalResult.pAtMin[0];
                            scaleSoil = (float) retrievalResult.pAtMin[1];
                            addConst = (float) retrievalResult.pAtMin[2];
                            sdr = retrievalResult.sdr;
                            model = retrievalResult.modelSpec;
                        }
*/
                    }
                }
                aotTile.setSample(iX, iY, aot);
                errTile.setSample(iX, iY, err);
                if (saveModelBands) {
                    scaleVegTile.setSample(iX, iY, scaleVeg);
                    scaleSoilTile.setSample(iX, iY, scaleSoil);
                    addConstTile.setSample(iX, iY, addConst);
                }
                for (int i=0; i<nSpecBands; i++) {
                    if (saveSdrBands)   sdrTiles[i].setSample(iX, iY, sdr[i]);
                    if (saveModelBands) modelTiles[i].setSample(iX, iY, model[i]);
                }
                if(pm.isCanceled()) return;
            }
            pm.worked(1);
            System.err.println(iY);
        }
        pm.done();
        System.err.println("tiledone");
    }

    private void createTargetProduct() {

        targetProduct = new Product(productName, productType, srcRasterWidth, srcRasterHeight);

        ProductUtils.copyMetadata(sourceProduct, targetProduct);
        ProductUtils.copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        ProductUtils.copyFlagBands(sourceProduct, targetProduct);
        //AerosolHelpers.addAerosolFlagBand(targetProduct, downscaledRasterWidth, downscaledRasterHeight);

        createTargetProductBands();

        ProductUtils.copyMasks(sourceProduct, targetProduct);
        ProductUtils.copyRoiMasks(sourceProduct, targetProduct);

        setTargetProduct(targetProduct);

    }

    private void createTargetProductBands() {

        List<String> specNameList = Arrays.asList(specBandNames);
        for (Band srcBand : sourceProduct.getBands()) {
            if (srcBand.isFlagBand()) {
                Band tarBand = targetProduct.getBand(srcBand.getName());
                tarBand.setSourceImage(srcBand.getSourceImage());
            }
            //else if (saveToaBands || !specNameList.contains(srcBand.getName())) {
            else if (saveToaBands && specNameList.contains(srcBand.getName())) {
                ProductUtils.copyBand(srcBand.getName(), sourceProduct, targetProduct);
                Band tarBand = targetProduct.getBand(srcBand.getName());
                tarBand.setSourceImage(srcBand.getSourceImage());
            }
        }

        Band targetBand = new Band("aot", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("best fitting aot Band");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        //targetBand.setValidPixelExpression(targetBand.getName() + ">= 0 AND " + targetBand.getName() + "<= 2");
        targetBand.setValidPixelExpression(validExpression);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        targetBand = new Band("aot_err", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
        targetBand.setDescription("");
        targetBand.setNoDataValue(-1);
        targetBand.setNoDataValueUsed(true);
        targetBand.setValidPixelExpression(validExpression);
        targetBand.setUnit("dl");
        targetProduct.addBand(targetBand);

        if (saveModelBands) {
            targetBand = new Band("scaleVeg", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setValidPixelExpression(validExpression);
            targetBand.setUnit("dl");
            targetProduct.addBand(targetBand);

            targetBand = new Band("scaleSoil", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setValidPixelExpression(validExpression);
            targetBand.setUnit("dl");
            targetProduct.addBand(targetBand);

            targetBand = new Band("addConst", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
            targetBand.setDescription("");
            targetBand.setNoDataValue(-1);
            targetBand.setNoDataValueUsed(true);
            targetBand.setValidPixelExpression(validExpression);
            targetBand.setUnit("dl");
            targetProduct.addBand(targetBand);

            for (int i = 0; i < nSpecBands; i++) {
                targetBand = new Band(specBandNames[i] + "_model", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
                targetBand.setDescription("");
                targetBand.setNoDataValue(-1);
                targetBand.setNoDataValueUsed(true);
                targetBand.setValidPixelExpression(validExpression);
                targetBand.setUnit("dl");
                targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
                targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
                targetProduct.addBand(targetBand);
            }
        }
        
        if (saveSdrBands) {
            for (int i = 0; i < nSpecBands; i++) {
                targetBand = new Band(specBandNames[i] + "_sdr", ProductData.TYPE_FLOAT32, srcRasterWidth, srcRasterHeight);
                targetBand.setDescription("");
                targetBand.setNoDataValue(-1);
                targetBand.setNoDataValueUsed(true);
                targetBand.setValidPixelExpression("aot_err<0.005");
                targetBand.setValidPixelExpression(validExpression);
                targetBand.setUnit("dl");
                targetBand.setSpectralBandwidth(sourceProduct.getBand(specBandNames[i]).getSpectralBandwidth());
                targetBand.setSpectralWavelength(sourceProduct.getBand(specBandNames[i]).getSpectralWavelength());
                targetProduct.addBand(targetBand);
            }
        }
    }

    private ArrayList<TiePointGrid> getGeomBandList(String instrument) {
        String[] bandNames = InstrumentConsts.getInstance().getGeomBandNames(instrument);
        ArrayList<TiePointGrid> list = new ArrayList<TiePointGrid>();
        for (int i=0; i<bandNames.length; i++) {
            list.add(sourceProduct.getTiePointGrid(bandNames[i]));
        }
        return list;
    }

    private ArrayList<Band> getSpecBandList(String instrument) {
        String[] bandNames = InstrumentConsts.getInstance().getSpecBandNames(instrument);
        ArrayList<Band> list = new ArrayList<Band>();
        for (int i=0; i<bandNames.length; i++) {
            list.add(sourceProduct.getBand(bandNames[i]));
        }
/*
        if (list.get(0).getSpectralWavelength() > 0) {
            Comparator<Band> byWavelength = new WavelengthComparator();
            Collections.sort(list, byWavelength);
        }
*/
        return list;
    }

    private float getRelativeAzi(float saa, float vaa) {
        float relAzi = Math.abs(saa - vaa);
        relAzi = (relAzi > 180.0f) ? 180 - (360 - relAzi) : 180 - relAzi;
        return relAzi;
    }

    private float[] getSpectralWvl(ArrayList<Band> specBandList) {
        int nBands = (instrument.equals("AATSR")) ? 4 : specBandList.size();
        float[] wvl = new float[nBands];
        for (int i=0; i<nBands; i++) {
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
                PixelGeometry geomFward = null;
                if (instrument.equals("AATSR")){
                    geom = new PixelGeometry(90f-sza, saa, 90f-vza, vaa);
                    sza = inputTiles.get(geomBandNames[4]).getSampleFloat(tileX, tileY);
                    saa = inputTiles.get(geomBandNames[5]).getSampleFloat(tileX, tileY);
                    vza = inputTiles.get(geomBandNames[6]).getSampleFloat(tileX, tileY);
                    vaa = inputTiles.get(geomBandNames[7]).getSampleFloat(tileX, tileY);
                    geomFward = new PixelGeometry(90f-sza, saa, 90f-vza, vaa);
                }
                valid[i+j*window.width] = (sza != (float)noData);

                noData = inputTiles.get(surfPresName).getRasterDataNode().getGeophysicalNoDataValue();
                float pres = inputTiles.get(surfPresName).getSampleFloat(tileX, tileY);
                valid[i+j*window.width] = valid[i+j*window.width] && (pres != noData);

                noData = inputTiles.get(surfPresName).getRasterDataNode().getGeophysicalNoDataValue();
                float o3col = inputTiles.get(ozoneName).getSampleFloat(tileX, tileY) * 1000; // O3 in [DU]
                float wvCol = 0.0f;
                valid[i+j*window.width] = valid[i+j*window.width] && (o3col != noData);

                float[] spec = new float[specWvl.length];
                float[] specFward = null;
                for (int k=0; k<specWvl.length; k++) {
                    Tile t = inputTiles.get(specBandNames[k]);
                    spec[k] = t.getSampleFloat(tileX, tileY);
                    valid[i+j*window.width] = valid[i+j*window.width] && (spec[k] != t.getRasterDataNode().getGeophysicalNoDataValue());
                }
                if (instrument.equals("AATSR")){
                    specFward = new float[specWvl.length];
                    for (int k=0; k<specWvl.length; k++) {
                        spec[k] = RsMathUtils.radianceToReflectance(spec[k], geom.sza, (float) Math.PI) * 0.01f;

                        Tile t = inputTiles.get(specBandNames[4+k]);
                        specFward[k] = t.getSampleFloat(tileX, tileY);
                        specFward[k] = RsMathUtils.radianceToReflectance(specFward[k], geomFward.sza, (float) Math.PI) * 0.01f;
                        valid[i+j*window.width] = valid[i+j*window.width] && (specFward[k] != t.getRasterDataNode().getGeophysicalNoDataValue());
                    }
                }
                inPixel[i+j*window.width] = new InputPixelData(geom, geomFward, pres, o3col, wvCol, specWvl, spec, specFward);
            }
        }
        
        return new SourcePixelField(inPixel, valid, (-window.x)+(-window.y)*window.width);
    }

    private Rectangle getPixelWindow(String pixelWindowStr) {
        String[] split = pixelWindowStr.split(",");
        Guardian.assertEquals("Number of Source Pixel Window parameters", split.length, 4);
        return new Rectangle(Integer.parseInt(split[0]), Integer.parseInt(split[1]),
                             Integer.parseInt(split[2]), Integer.parseInt(split[3]));
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
