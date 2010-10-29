/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.util.HashMap;
import java.util.Map;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.idepix.operators.CloudScreeningSelector;
import org.esa.beam.idepix.operators.ComputeChainOp;
import org.esa.beam.meris.brr.Rad2ReflOp;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

/**
 * Create Meris input product for Globalbedo aerosol retrieval and BBDR processor
 * @author akheckel
 */
@OperatorMetadata(alias = "MerisPrepOp",
                  description = "Create Meris product for input to Globalbedo aerosol retrieval and BBDR processor",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class MerisPrepOp extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;

    @Override
    public void initialize() throws OperatorException {
        InstrumentConsts instrC = InstrumentConsts.getInstance();
        final boolean needPixelClassif = (!sourceProduct.containsBand(instrC.getIdepixFlagBandName()));
        final boolean needElevation = (!sourceProduct.containsBand(instrC.getElevationBandName()));
        final boolean needSurfacePres = (!sourceProduct.containsBand(instrC.getSurfPressureName("MERIS")));

        // convert radiance bands to reflectance
        Product reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(Rad2ReflOp.class), GPF.NO_PARAMS, sourceProduct);

        // subset might have set ptype to null, thus:
        if (sourceProduct.getDescription() == null) sourceProduct.setDescription("MERIS Radiance product");

        // setup target product primarily as copy of sourceProduct
        final int rasterWidth = sourceProduct.getSceneRasterWidth();
        final int rasterHeight = sourceProduct.getSceneRasterHeight();
        targetProduct = new Product(sourceProduct.getName(),
                                    sourceProduct.getProductType(),
                                    rasterWidth, rasterHeight);
        targetProduct.setStartTime(sourceProduct.getStartTime());
        targetProduct.setEndTime(sourceProduct.getEndTime());
        targetProduct.setPointingFactory(sourceProduct.getPointingFactory());
        ProductUtils.copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        ProductUtils.copyFlagBands(sourceProduct, targetProduct);

        // create pixel calssification if missing in sourceProduct
        // and add flag band to targetProduct
        Product idepixProduct = null;
        if (needPixelClassif) {
            Map<String, Object> pixelClassParam = new HashMap<String, Object>(4);
            pixelClassParam.put("algorithm", CloudScreeningSelector.GlobAlbedo);
            pixelClassParam.put("gaCopyRadiances", false);
            pixelClassParam.put("gaCopyAnnotations", false);
            pixelClassParam.put("gaComputeFlagsOnly", true);
            idepixProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(ComputeChainOp.class), pixelClassParam, sourceProduct);
            ProductUtils.copyFlagBands(idepixProduct, targetProduct);
        }

        // create elevation product if band is missing in sourceProduct
        Product elevProduct = null;
        if (needElevation){
            elevProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(CreateElevationBandOp.class), GPF.NO_PARAMS, sourceProduct);
        }

        // create surface pressure estimate product if band is missing in sourceProduct
        VirtualBand surfPresBand = null;
        if (needSurfacePres){
            String presExpr = "(1013.25 * exp(-elevation/8400))";
            surfPresBand = new VirtualBand(instrC.getSurfPressureName("MERIS"),
                                                       ProductData.TYPE_FLOAT32,
                                                       rasterWidth, rasterHeight, presExpr);
            surfPresBand.setDescription("estimated sea level pressure (p0=1013.25hPa, hScale=8.4km)");
            surfPresBand.setNoDataValue(0);
            surfPresBand.setNoDataValueUsed(true);
            surfPresBand.setUnit("hPa");
        }



        // copy all non-radiance bands from sourceProduct and
        // copy reflectance bands from reflProduct
        Band tarBand;
        for (Band srcBand : sourceProduct.getBands()){
            String srcName = srcBand.getName();
            if (srcBand.isFlagBand()){
                tarBand = targetProduct.getBand(srcName);
                tarBand.setSourceImage(srcBand.getSourceImage());
            }
            else if (srcName.startsWith("radiance")){
                String reflName = "rho_toa_" + srcName.split("_")[1];
                String tarName = "reflectance_" + srcName.split("_")[1];
                tarBand = ProductUtils.copyBand(reflName, reflProduct, targetProduct);
                tarBand.setName(tarName);
                Band reflBand = reflProduct.getBand(reflName);
                tarBand.setSourceImage(reflBand.getSourceImage());
            }
            else {
                tarBand = ProductUtils.copyBand(srcName, sourceProduct, targetProduct);
                tarBand.setSourceImage(srcBand.getSourceImage());

            }
        }

        // add idepix flag band data if needed
        if (needPixelClassif){
            Guardian.assertNotNull("idepixProduct", idepixProduct);
            Band srcBand = idepixProduct.getBand(instrC.getIdepixFlagBandName());
            Guardian.assertNotNull("idepix Band", srcBand);
            tarBand = targetProduct.getBand(srcBand.getName());
            tarBand.setSourceImage(srcBand.getSourceImage());
        }

        // add elevation band if needed
        if (needElevation){
            Guardian.assertNotNull("elevProduct", elevProduct);
            Band srcBand = elevProduct.getBand(instrC.getElevationBandName());
            Guardian.assertNotNull("elevation band", srcBand);
            tarBand = ProductUtils.copyBand(srcBand.getName(), elevProduct, targetProduct);
            tarBand.setSourceImage(srcBand.getSourceImage());
        }

        // add vitrual surface pressure band if needed
        if (needSurfacePres) {
            targetProduct.addBand(surfPresBand);
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
            super(MerisPrepOp.class);
        }
    }
}
