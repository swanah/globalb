/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.awt.Dimension;
import java.awt.RenderingHints;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.dataio.envisat.EnvisatConstants;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

/**
 * Main Operator producing AOT for GlobAlbedo
 */
@OperatorMetadata(alias = "ga.MasterOp",
                  description = "",
                  authors = "Andreas Heckel",
                  version = "1.1",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class GaMasterOp  extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;
    @Parameter(defaultValue="false")
    private boolean copyToaRadBands;
    @Parameter(defaultValue="false")
    private boolean copyToaReflBands;
/*
    @Parameter(defaultValue="true")
    private boolean retrieveAOT = true;
    @Parameter(defaultValue="false")
    private boolean saveSdrBands = false;
    @Parameter(defaultValue="false")
    private boolean saveModelBands = false;
*/
    @Parameter(defaultValue="2")
    private int vegSpecId;
    private int scale;
    private String instrument;

    @Override
    public void initialize() throws OperatorException {
        scale = 9;
        //Dimension prefTileSize = sourceProduct.getPreferredTileSize();
        int xTileSize = (sourceProduct.getSceneRasterWidth() > 1000)? 1000 : sourceProduct.getSceneRasterWidth();
        int yTileSize = 9;
        Dimension prefTileSize = new Dimension(xTileSize, yTileSize);
        Dimension smallTileSize = new Dimension(prefTileSize.width/scale+1, prefTileSize.height/scale+1);

        final boolean isMerisProduct = sourceProduct.getProductType().equals(EnvisatConstants.MERIS_RR_L1B_PRODUCT_TYPE_NAME);
        final boolean isAatsrProduct = sourceProduct.getProductType().equals(EnvisatConstants.AATSR_L1B_TOA_PRODUCT_TYPE_NAME);
        final boolean isVgtProduct = sourceProduct.getProductType().startsWith("VGT PRODUCT FORMAT V1.");

        Guardian.assertTrue("is valid source product", (isMerisProduct ^ isAatsrProduct ^ isVgtProduct));

        Product reflProduct = null;
        if (isMerisProduct) {
            instrument = "MERIS";
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(MerisPrepOp.class), GPF.NO_PARAMS, sourceProduct);
        }
        else if (isAatsrProduct) {
            instrument = "AATSR";
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AatsrPrepOp.class), GPF.NO_PARAMS, sourceProduct);;
        }
        else if (isVgtProduct) {
            instrument = "VGT";
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(VgtPrepOp.class), GPF.NO_PARAMS, sourceProduct);
        }


        //Map<String, Object> sclParams = new HashMap<String, Object>(1);
        //sclParams.put("scale", scale);
        //Product downsclProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(DownSclOp.class), sclParams, reflProduct);

        Map<String, Object> aotParams = new HashMap<String, Object>(4);
        aotParams.put("vegSpecId", vegSpecId);
        aotParams.put("scale", scale);
        /*
        aotParams.put("retrieveAOT", retrieveAOT);
        aotParams.put("saveToaBands", saveToaBands);
        aotParams.put("saveSdrBands", saveSdrBands);
        aotParams.put("saveModelBands", saveModelBands);
         *
         */

        RenderingHints rh = new RenderingHints(GPF.KEY_TILE_SIZE, smallTileSize);
        Product aotDownsclProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AerosolOp2.class), aotParams, reflProduct, rh);

        Map<String, Product> fillSourceProds = new HashMap<String, Product>(2);
        fillSourceProds.put("aotProduct", aotDownsclProduct);
        fillSourceProds.put("fillProduct", aotDownsclProduct);
        Product fillAotProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(GapFillingOp.class), GPF.NO_PARAMS, fillSourceProds, rh);

        Map<String, Product> upsclProducts = new HashMap<String, Product>(2);
        upsclProducts.put("lowresProduct", fillAotProduct);
        upsclProducts.put("hiresProduct", reflProduct);
        Map<String, Object> sclParams = new HashMap<String, Object>(1);
        sclParams.put("scale", scale);
        rh = new RenderingHints(GPF.KEY_TILE_SIZE, prefTileSize);
        Product aotHiresProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(UpSclOp.class), sclParams, upsclProducts, rh);

        targetProduct = mergeToTargetProduct(reflProduct, aotHiresProduct);
        setTargetProduct(targetProduct);
    }

    private Product mergeToTargetProduct(Product reflProduct, Product aotHiresProduct) {
        String pname = sourceProduct.getName() + "_AOT";
        String ptype = sourceProduct.getProductType() + " GlobAlbedo AOT";
        int rasterWidth = sourceProduct.getSceneRasterWidth();
        int rasterHeight = sourceProduct.getSceneRasterHeight();
        Product tarP = new Product(pname, ptype, rasterWidth, rasterHeight);
        tarP.setStartTime(sourceProduct.getStartTime());
        tarP.setEndTime(sourceProduct.getEndTime());
        tarP.setPointingFactory(sourceProduct.getPointingFactory());
        ProductUtils.copyMetadata(aotHiresProduct, tarP);
        ProductUtils.copyTiePointGrids(sourceProduct, tarP);
        ProductUtils.copyGeoCoding(sourceProduct, tarP);
        ProductUtils.copyFlagBands(reflProduct, tarP);
        ProductUtils.copyFlagBands(aotHiresProduct, tarP);
        Band tarBand;
        String bname;
        if (copyToaRadBands){
            for (Band b : sourceProduct.getBands()){
                bname = b.getName();
                if (b.getSpectralWavelength() > 0){
                    tarBand = ProductUtils.copyBand(bname, sourceProduct, tarP);
                    tarBand.setSourceImage(b.getSourceImage());
                }
            }
        }
        for (Band b : reflProduct.getBands()){
            bname = b.getName();
            if (b.isFlagBand()){
                tarBand = tarP.getBand(bname);
                tarBand.setSourceImage(b.getSourceImage());
            }
            boolean copyBand = (copyToaReflBands && !tarP.containsBand(bname) && b.getSpectralWavelength() > 0);
            copyBand = copyBand || (instrument.equals("VGT") && InstrumentConsts.getInstance().isVgtAuxBand(b));
            copyBand = copyBand || (bname.equals("elevation"));

            if (copyBand){
                tarBand = ProductUtils.copyBand(bname, reflProduct, tarP);
                tarBand.setSourceImage(b.getSourceImage());
            }
        }
        for (Band b : aotHiresProduct.getBands()){
            bname = b.getName();
            if (b.isFlagBand()){
                tarBand = tarP.getBand(bname);
            }
            else {
                tarBand = ProductUtils.copyBand(bname, aotHiresProduct, tarP);
            }
            tarBand.setSourceImage(b.getSourceImage());
        }
        return tarP;
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
            super(GaMasterOp.class);
        }
    }
}
