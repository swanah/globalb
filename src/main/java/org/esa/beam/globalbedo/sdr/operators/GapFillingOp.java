/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import com.bc.ceres.core.ProgressMonitor;
import java.awt.Rectangle;
import javax.media.jai.BorderExtender;
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

/**
 *
 * @author akheckel
 */
@OperatorMetadata(alias = "ga.GapFillingOp",
                  description = "Fills Gaps in Grid",
                  authors = "Andreas Heckel",
                  version = "1.1",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class GapFillingOp extends Operator {

    @SourceProduct
    private Product aotProduct;
    @SourceProduct
    private Product fillProduct;
    @TargetProduct
    private Product targetProduct;
    private int box;
    private int off;
    private int rasterHeight;
    private int rasterWidth;
    private int climDist;
    private double climAot;
    private int climErr;


    @Override
    public void initialize() throws OperatorException {
        off = 25;
        box = off*2+1;
        climAot = 0.07;
        climErr = 2;
        climDist = 10;
        String pname = fillProduct.getName();
        String ptype = fillProduct.getProductType();
        rasterWidth = fillProduct.getSceneRasterWidth();
        rasterHeight = fillProduct.getSceneRasterHeight();
        targetProduct = new Product(pname, ptype, rasterWidth, rasterHeight);
        ProductUtils.copyBand("aot", fillProduct, targetProduct);
        ProductUtils.copyBand("aot_err", fillProduct, targetProduct);
        setTargetProduct(targetProduct);
    }

    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        String bandName = targetBand.getName();
        Rectangle tarRec = targetTile.getRectangle();
//System.err.println("gapfilling "+targetBand.getName()+" on tile: " + tarRec);
        Rectangle srcRec = new Rectangle(tarRec.x-off, tarRec.y-off, tarRec.width+box, tarRec.height+box);
        Tile aotTile = getSourceTile(aotProduct.getBand(bandName), srcRec, BorderExtender.createInstance(BorderExtender.BORDER_ZERO), ProgressMonitor.NULL);
        Tile fillTile = getSourceTile(fillProduct.getBand(bandName), srcRec, BorderExtender.createInstance(BorderExtender.BORDER_ZERO), ProgressMonitor.NULL);
        double noDataVal = fillTile.getRasterDataNode().getGeophysicalNoDataValue();
        double climVal = 0;
        if (bandName.equals("aot")) climVal = climAot;
        if (bandName.equals("aot_err")) climVal = climErr;
        float aotPixel;
        for (int y = tarRec.y; y < tarRec.y+tarRec.height; y++){
            for (int x = tarRec.x; x < tarRec.x+tarRec.width; x++){
                aotPixel = aotTile.getSampleFloat(x, y);
                if (Double.compare(noDataVal, aotPixel) != 0){
                    targetTile.setSample(x, y, aotPixel);
                }
                else {
                    targetTile.setSample(x, y, getAve(x, y, fillTile, noDataVal, climVal));
                }
            }
        }
    }

    private float getAve(int x, int y, Tile srcTile, double noDataValue, double climVal){
        double n0 = invDistanceWeight(climDist, 0, 4);
        double n = n0;
        double sum = climVal * n;
        float val;
        double weight;
        int ys = (y-off < 0) ? 0 : y-off;
        int xs = (x-off < 0) ? 0 : x-off;
        int ye = (y+off >= rasterHeight) ? rasterHeight-1 : y+off;
        int xe = (x+off >= rasterWidth) ? rasterWidth-1 : x+off;
        for (int j=ys; j<=ye; j++){
            for (int i=xs; i<=xe; i++){
                val = srcTile.getSampleFloat(i, j);
                if (Double.compare(noDataValue, val) != 0){
                    weight = invDistanceWeight(i-x, j-y, 4);
                    sum += val*weight;
                    n += weight;
                }
            }
        }
        return (float) ((n > 0) ? sum / n : noDataValue);
    }

    private double invDistanceWeight(int i, int j, int power) {
        return 1 / (Math.pow(Math.pow(i, 2) + Math.pow(j, 2), power/2) + 1.0E-5);
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
            super(GapFillingOp.class);
        }
    }
}
