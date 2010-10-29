/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.util.HashMap;
import java.util.Map;
import org.esa.beam.dataio.envisat.EnvisatConstants;
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

/**
 * Main Operator producing AOT for GlobAlbedo
 */
@OperatorMetadata(alias = "MasterOp",
                  description = "",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class MasterOp  extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;
/*
    @Parameter(defaultValue="true")
    private boolean retrieveAOT = true;
    @Parameter(defaultValue="false")
    private boolean saveToaBands = false;
    @Parameter(defaultValue="false")
    private boolean saveSdrBands = false;
    @Parameter(defaultValue="false")
    private boolean saveModelBands = false;
*/
    @Parameter(defaultValue="2")
    private int vegSpecId;
    private int scale;

    @Override
    public void initialize() throws OperatorException {

        final boolean isMerisProduct = sourceProduct.getProductType().equals(EnvisatConstants.MERIS_RR_L1B_PRODUCT_TYPE_NAME);
        final boolean isAatsrProduct = sourceProduct.getProductType().equals(EnvisatConstants.AATSR_L1B_TOA_PRODUCT_TYPE_NAME);
        final boolean isVgtProduct = sourceProduct.getProductType().startsWith("VGT PRODUCT FORMAT V1.");

        Guardian.assertTrue("is valid source product", (isMerisProduct ^ isAatsrProduct ^ isVgtProduct));

        Product reflProduct = null;
        if (isMerisProduct) {
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(MerisPrepOp.class), GPF.NO_PARAMS, sourceProduct);
        }
        else if (isAatsrProduct) {
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AatsrPrepOp.class), GPF.NO_PARAMS, sourceProduct);;
        }
        else if (isVgtProduct) {
            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(VgtPrepOp.class), GPF.NO_PARAMS, sourceProduct);
        }


        scale = 9;
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
        Product aotDownsclProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AerosolOp2.class), aotParams, reflProduct);
        //Product aotDownsclProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AerosolOp2.class), GPF.NO_PARAMS, elevProduct);

        Map<String, Product> upsclProducts = new HashMap<String, Product>(2);
        upsclProducts.put("lowresProduct", aotDownsclProduct);
        upsclProducts.put("hiresProduct", reflProduct);
        Map<String, Object> sclParams = new HashMap<String, Object>(1);
        sclParams.put("scale", scale);
        targetProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(UpSclOp.class), sclParams, upsclProducts);
        //targetProduct = aotDownsclProduct;
        //targetProduct = reflProduct;
        setTargetProduct(targetProduct);
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
            super(MasterOp.class);
        }
    }
}
