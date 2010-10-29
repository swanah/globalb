/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import java.awt.Rectangle;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.gpf.operators.standard.SubsetOp;

/**
 *
 * @author akheckel
 */
@OperatorMetadata(alias = "SzaSubsetOp",
                  description = "upscaling product",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class SzaSubsetOp extends Operator{

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;

    @Override
    public void initialize() throws OperatorException {
        TiePointGrid szaTPG = sourceProduct.getTiePointGrid("sun_elev_nadir");
        int yStart=0;
        int yEnd=sourceProduct.getSceneRasterHeight()-1;
        for (int i=0; i<sourceProduct.getSceneRasterHeight(); i++){
            boolean foundStart = (szaTPG.getPixelFloat(0, yStart)>30);
            boolean foundEnd = (szaTPG.getPixelFloat(0, yEnd)>30);
            if (!foundStart) yStart++;
            if (!foundEnd) yEnd--;
            if (foundEnd && foundStart) break;
        }

        String[] tpNames = sourceProduct.getTiePointGridNames();
        String[] bandNames = sourceProduct.getBandNames();
        Map<String, Object> subsetParams = new HashMap<String, Object>(6);
        subsetParams.put("region", new Rectangle(0, yStart, sourceProduct.getSceneRasterWidth(), yEnd-yStart+1));
        subsetParams.put("subSamplingX", 1);
        subsetParams.put("subSamplingY", 1);
        subsetParams.put("tiePointGridNames", tpNames);
        subsetParams.put("bandNames", bandNames);
        subsetParams.put("copyMetadata", true);
        targetProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(SubsetOp.class), subsetParams, sourceProduct);
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
            super(SzaSubsetOp.class);
        }
    }
}
