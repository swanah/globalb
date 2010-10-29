/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import com.bc.ceres.core.ProgressMonitor;
import java.awt.Rectangle;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.gpf.operators.standard.BandMathsOp;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

/**
 *
 * @author akheckel
 */
@OperatorMetadata(alias = "DownSclOp",
                  description = "downscaling product",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class DownSclOp extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;

    @Parameter(defaultValue="9")
    private int scale;
    private int offset;
    private int targetWidth;
    private int targetHeight;
    private int sourceWidth;
    private int sourceHeight;
    private String instrument;
    private String validExpression;
    private Band validBand;

    public DownSclOp() {
    }

    @Override
    public void initialize() throws OperatorException {
        offset = scale/2;
        sourceWidth = sourceProduct.getSceneRasterWidth();
        sourceHeight = sourceProduct.getSceneRasterHeight();
        targetWidth = (int) ((sourceWidth + 0.5) / scale);
        targetHeight = (int) ((sourceHeight + 0.5) / scale);

        targetProduct = new Product(sourceProduct.getName(), sourceProduct.getProductType(), targetWidth, targetHeight);
        ProductUtils.copyMetadata(sourceProduct, targetProduct);
        copyTiePointGrids(sourceProduct, targetProduct);
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        ProductUtils.copyFlagCodings(sourceProduct, targetProduct);
        /*
         * I can't copy masks because resolution is changing
        ProductUtils.copyMasks(sourceProduct, targetProduct);
        ProductUtils.copyOverlayMasks(sourceProduct, targetProduct);
        ProductUtils.copyRoiMasks(sourceProduct, targetProduct);
        */
        copyBands(sourceProduct, targetProduct);


        targetProduct.setStartTime(sourceProduct.getStartTime());
        targetProduct.setEndTime(sourceProduct.getEndTime());
        targetProduct.setPointingFactory(sourceProduct.getPointingFactory());

        instrument = InstrumentConsts.getInstance().getInstrument(sourceProduct);
        validExpression = InstrumentConsts.getInstance().getValidExpression(instrument);
        final BandMathsOp validBandOp = BandMathsOp.createBooleanExpressionBand(validExpression, sourceProduct);
        validBand = validBandOp.getTargetProduct().getBandAt(0);
    }

    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        Rectangle tarRec = targetTile.getRectangle();
        int srcX0 = tarRec.x * scale;
        int srcY0 = tarRec.y * scale;
        int srcWidth = tarRec.width * scale;
        int srcHeight = tarRec.height * scale;
        if (srcX0+srcWidth > sourceWidth) srcWidth = (sourceWidth - srcX0);
        if (srcY0+srcHeight > sourceHeight) srcHeight = (sourceHeight - srcY0);
        Rectangle srcRec = new Rectangle(srcX0, srcY0, srcWidth, srcHeight);

        Tile srcTile = getSourceTile(sourceProduct.getBand(targetBand.getName()), srcRec, ProgressMonitor.NULL);
        Tile validTile = getSourceTile(validBand, srcRec, ProgressMonitor.NULL);
        double noData = targetBand.getGeophysicalNoDataValue();

        for (int ity = tarRec.y; ity < tarRec.y+tarRec.height; ity++){
            srcY0 = ity * scale;
            srcHeight = (ity+1)*scale;
            if (srcHeight > sourceHeight) srcHeight = sourceHeight - srcY0;
            for (int itx = tarRec.x; itx < tarRec.x+tarRec.width; itx++){
                srcX0 = itx * scale;
                srcWidth = (itx+1)*scale;
                if (srcWidth > sourceWidth) srcWidth = sourceWidth - srcX0;

                double sum = 0;
                int n = 0;
                for (int isy = srcY0; isy < srcHeight; isy++){
                    for (int isx = srcX0; isx < srcWidth; isx++){
                        double val = srcTile.getSampleDouble(isx, isy);
                        if (Double.compare(val, noData) != 0
                                && validTile.getSampleBoolean(isx, isy)){
                            sum += val;
                            n++;
                        }
                    }
                }

                double ave = (n>scale*scale*0.95)? sum/n : noData;
                targetTile.setSample(itx, ity, ave);
            }
        }
    }

    private void copyTiePointGrids(Product sourceProduct, Product targetProduct) {
        for (String tpgName : sourceProduct.getTiePointGridNames()){
            float[] tpgData = new float[targetHeight * targetWidth];
            TiePointGrid srcTpg = sourceProduct.getTiePointGrid(tpgName);
            for (int ity=0; ity<targetHeight; ity++){
                int isy = (ity) * scale + offset;
                //float isy = (ity) * scale;
                for (int itx=0; itx<targetWidth; itx++){
                    int isx = (itx) * scale + offset;
                    //float isx = (itx) * scale;
                    tpgData[ity*targetWidth+itx] = srcTpg.getPixelFloat(isx, isy);
                }
            }
            TiePointGrid tpg = new TiePointGrid(tpgName, targetWidth, targetHeight, 0.5f, 0.5f, 1, 1, tpgData);
            targetProduct.addTiePointGrid(tpg);
        }
    }

    private void copyBands(Product sourceProduct, Product targetProduct) {
        String[] bandNames = sourceProduct.getBandNames();
        for (String bn : bandNames){
            Band bSrc = sourceProduct.getBand(bn);
            Band b = copyBandScl(bn, sourceProduct, bn, targetProduct);
            if (bSrc.isFlagBand()){
                b.setSampleCoding(targetProduct.getFlagCodingGroup().get(bSrc.getFlagCoding().getName()));
            }
        }
    }

    /**
     * Copies the properties of the named band from the source product to the target product.<br>
     * <b>source and target product can have different image size</b>
     *
     * @param sourceBandName the name of the band to be copied.
     * @param sourceProduct  the source product.
     * @param targetBandName the name of the band copied.
     * @param targetProduct  the target product.
     *
     * @return the copy of the band, or <code>null</code> if the sourceProduct does not contain a band with the given name.
     */
    private Band copyBandScl(String sourceBandName, Product sourceProduct,
                                String targetBandName, Product targetProduct) {
        Guardian.assertNotNull("sourceProduct", sourceProduct);
        Guardian.assertNotNull("targetProduct", targetProduct);

        if (sourceBandName == null || sourceBandName.length() == 0) {
            return null;
        }
        final Band sourceBand = sourceProduct.getBand(sourceBandName);
        if (sourceBand == null) {
            return null;
        }
        Band targetBand = new Band(targetBandName,
                                   sourceBand.getDataType(),
                                   targetProduct.getSceneRasterWidth(),
                                   targetProduct.getSceneRasterHeight());
        ProductUtils.copyRasterDataNodeProperties(sourceBand, targetBand);
        targetProduct.addBand(targetBand);
        return targetBand;
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
            super(DownSclOp.class);
        }
    }
}
