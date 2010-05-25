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
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.synergy.operators.CreateMerisOp;

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


    @Override
    public void initialize() throws OperatorException {
        Product reflProduct = null;
        if (sourceProduct.getProductType().equals(EnvisatConstants.MERIS_RR_L1B_PRODUCT_TYPE_NAME)) {
            // MERIS product
            Map<String, Object> merisParams = new HashMap<String, Object>(3);
            merisParams.put("copyToaRadiances", false);
            merisParams.put("copyCloudProbability", false);
            merisParams.put("copyCloudTopPreassureAndMask", false);
            merisParams.put("copyLandWaterReclass", false);

            reflProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(CreateMerisOp.class), merisParams, sourceProduct);
        }


        targetProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(AerosolOp.class), GPF.NO_PARAMS, reflProduct);
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
