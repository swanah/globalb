/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.esa.beam.globalbedo.sdr.operators;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.glevel.MultiLevelImage;
import java.awt.Rectangle;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferFloat;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.media.jai.BorderExtender;
import javax.media.jai.PlanarImage;
import javax.media.jai.RasterFactory;
import javax.media.jai.TiledImage;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.FlagCoding;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.ProductNodeGroup;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

/**
 *
 * @author akheckel
 */
@OperatorMetadata(alias = "UpSclOp",
                  description = "upscaling product",
                  authors = "Andreas Heckel",
                  version = "1.0",
                  copyright = "(C) 2010 by University Swansea (a.heckel@swansea.ac.uk)")
public class UpSclOp extends Operator {

    @SourceProduct
    private Product lowresProduct;
    @SourceProduct
    private Product hiresProduct;
    @TargetProduct
    private Product targetProduct;

    @Parameter(defaultValue="9")
    private int scale;
    private int offset;
    private int targetWidth;
    private int targetHeight;
    private int sourceRasterWidth;
    private int sourceRasterHeight;
    private String targetProductName;
    private String targetProductType;
    private ArrayList<String> hiresBandNames;
    private ArrayList<String> lowresBandNames;
    private Product fillLowresProduct;
    //private Map<String, Band> filledLowresBands;

    @Override
    public void initialize() throws OperatorException {
        sourceRasterWidth = lowresProduct.getSceneRasterWidth();
        sourceRasterHeight = lowresProduct.getSceneRasterHeight();
        targetWidth = hiresProduct.getSceneRasterWidth();
        targetHeight = hiresProduct.getSceneRasterHeight();

        offset = scale/2;

        targetProductName = lowresProduct.getName();
        targetProductType = hiresProduct.getProductType();
        targetProduct = new Product(targetProductName, targetProductType, targetWidth, targetHeight);
        targetProduct.setStartTime(hiresProduct.getStartTime());
        targetProduct.setEndTime(hiresProduct.getEndTime());
        targetProduct.setPointingFactory(hiresProduct.getPointingFactory());
        ProductUtils.copyMetadata(lowresProduct, targetProduct);
        ProductUtils.copyTiePointGrids(hiresProduct, targetProduct);
        ProductUtils.copyGeoCoding(hiresProduct, targetProduct);
        ProductUtils.copyFlagCodings(hiresProduct, targetProduct);
        addFlagCodings(lowresProduct, targetProduct);
        hiresBandNames = copyBands(hiresProduct, targetProduct);
        lowresBandNames = copyBands(lowresProduct, targetProduct);

        Map<String, Product> fillSourceProds = new HashMap<String, Product>(2);
        fillSourceProds.put("aotProduct", lowresProduct);
        fillSourceProds.put("fillProduct", lowresProduct);
        fillLowresProduct = GPF.createProduct(OperatorSpi.getOperatorAlias(GapFillingOp.class), GPF.NO_PARAMS, fillSourceProds);
        //filledLowresBands = generateGapFilledBands(lowresProduct, lowresBandNames);

        setTargetProduct(targetProduct);
    }

    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        Rectangle tarRec = targetTile.getRectangle();
        String targetBandName = targetBand.getName();
        Band sourceBand;
        Tile sourceTile;

        if (hiresBandNames.contains(targetBandName)){
            sourceBand = hiresProduct.getBand(targetBandName);
            sourceTile = getSourceTile(sourceBand, tarRec, ProgressMonitor.NULL);
            targetTile.setRawSamples(sourceTile.getRawSamples());
        }
        else if(lowresBandNames.contains(targetBandName)){
            final Rectangle srcRec = calcSourceRectangle(tarRec);

            if (fillLowresProduct.containsBand(targetBandName)){
                sourceBand = fillLowresProduct.getBand(targetBandName);
                sourceTile = getSourceTile(sourceBand, srcRec, ProgressMonitor.NULL);
                upscaleTileBilinear(sourceTile, targetTile, tarRec, ProgressMonitor.NULL);
            }
            else {
                sourceBand = lowresProduct.getBand(targetBandName);
                sourceTile = getSourceTile(sourceBand, srcRec, ProgressMonitor.NULL);
                upscaleTileCopy(sourceTile, targetTile, tarRec, ProgressMonitor.NULL);
            }
        }
    }

    private Rectangle calcSourceRectangle(Rectangle tarRec) {
        int srcX = (tarRec.x - offset) / scale;
        int srcY = (tarRec.y - offset) / scale;
        int srcWidth = tarRec.width / scale + 2;
        int srcHeight = tarRec.height / scale + 2;
        if (srcX >= sourceRasterWidth) {
            srcX = sourceRasterWidth - 2;
            srcWidth = 2;
        }
        if (srcY >= sourceRasterHeight) {
            srcY = sourceRasterHeight - 2;
            srcHeight = 2;
        }
        if (srcX + srcWidth > sourceRasterWidth) {
            srcWidth = sourceRasterWidth - srcX;
        }
        if (srcY + srcHeight > sourceRasterHeight) {
            srcHeight = sourceRasterHeight - srcY;
        }
        return new Rectangle(srcX, srcY, srcWidth, srcHeight);
    }

    private void addFlagCodings(Product sourceProduct, Product targetProduct) {
        Guardian.assertNotNull("source", sourceProduct);
        Guardian.assertNotNull("target", targetProduct);

        int numCodings = sourceProduct.getFlagCodingGroup().getNodeCount();
        for (int n = 0; n < numCodings; n++) {
            FlagCoding sourceFlagCoding = sourceProduct.getFlagCodingGroup().get(n);
            if ( ! targetProduct.getFlagCodingGroup().contains(sourceFlagCoding.getName())) {
                ProductUtils.copyFlagCoding(sourceFlagCoding, targetProduct);
            }
        }
    }

    private ArrayList<String> copyBands(Product sourceProduct, Product targetProduct) {
        Guardian.assertNotNull("source", sourceProduct);
        Guardian.assertNotNull("target", targetProduct);
        ArrayList<String> bandNames = new ArrayList<String>();
        Band sourceBand;
        Band targetBand;
        FlagCoding flgCoding;
        ProductNodeGroup<FlagCoding> targetFCG = targetProduct.getFlagCodingGroup();

        for (int iBand = 0; iBand < sourceProduct.getNumBands(); iBand++) {
            sourceBand = sourceProduct.getBandAt(iBand);
            //System.err.println(sourceBand.getName());
            boolean validBand = sourceBand.isFlagBand();
            validBand = validBand || (sourceBand.getProduct().equals(hiresProduct) && InstrumentConsts.getInstance().isToaBand(sourceBand));
            validBand = validBand || (sourceBand.getProduct().equals(hiresProduct) && InstrumentConsts.getInstance().isElevationBand(sourceBand));
            validBand = validBand || (sourceBand.getProduct().equals(lowresProduct) && sourceBand.getName().startsWith("aot"));
            if (!targetProduct.containsBand(sourceBand.getName()) && validBand) {
                targetBand = copyBandScl(sourceBand.getName(), sourceProduct, sourceBand.getName(), targetProduct);
                bandNames.add(targetBand.getName());
                if (sourceBand.isFlagBand()) {
                    flgCoding = sourceBand.getFlagCoding();
                    if (!targetFCG.contains(flgCoding.getName())) {
                        ProductUtils.copyFlagCoding(flgCoding, targetProduct);
                    }
                    targetBand.setSampleCoding(targetFCG.get(flgCoding.getName()));
                }
            }
        }
        return bandNames;
    }

    /**
     * Copies the named band from the source product to the target product.
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

    private Map<String, Band> generateGapFilledBands(Product sourceProduct, ArrayList<String> bandNames) {
        Map<String,Band> filledBands = new HashMap<String, Band>(bandNames.size());
        Band srcBand;
        Band filledBand = null;
        for (int iBand = 0; iBand < bandNames.size(); iBand++){
            srcBand = sourceProduct.getBand(bandNames.get(iBand));
            if (!srcBand.isFlagBand()){
                filledBand = fillGaps(srcBand);
                filledBands.put(bandNames.get(iBand), filledBand);
            }
        }
        return filledBands;
    }

    private Band fillGaps(Band srcBand) {
        int width = srcBand.getSceneRasterWidth();
        int height = srcBand.getSceneRasterHeight();
        float[] srcData = getDataFloat(srcBand);
        //System.err.println("filling...");
        float[] fillData = fillGaps(srcData, width, height, (float) srcBand.getNoDataValue());
        //System.err.println("...done");

        return createBandFromData(srcBand.getName(), width, height, fillData);
    }

    private float[] getDataFloat(Band band) {
        MultiLevelImage sourceImage = band.getSourceImage();
        DataBuffer db = sourceImage.getData().getDataBuffer();
        if (db.getDataType() != DataBuffer.TYPE_FLOAT) {
            throw new OperatorException("wrong datatype getting data array of band "+band.getName());
        }
        float[] srcData = ((DataBufferFloat) db).getData();
        return srcData;
    }

    private float[] fillGaps(float[] data, int width, int height, float noData) {
        float[] fillData1 = new float[height * width];
        float[] fillData2 = new float[height * width];
        System.arraycopy(data, 0, fillData1, 0, data.length);

        int searchSize = 10;
        double total_old = 0;
        double total_new = 0;
        for (int i=0; i<data.length; i++) {
            if (Float.compare(data[i],noData)!=0) total_new+=data[i];
        }
        double threshold = 5e-7*total_new;
        while (Math.abs(total_new - total_old) > threshold) {
            total_old = total_new;
            total_new = 0;
            if (searchSize > 1) {
                searchSize *= 0.8;
            }
            int[] sr = {searchSize, searchSize, searchSize, searchSize}; // yl,yu,xl,xu
            for (int iy = 0; iy < height; iy++) {
                sr[0] = (iy < searchSize) ? iy : searchSize;
                sr[1] = (iy + searchSize >= height) ? height - 1 - iy : searchSize;
                for (int ix = 0; ix < width; ix++) {
                    sr[2] = (ix < searchSize) ? ix : searchSize;
                    sr[3] = (ix + searchSize >= width) ? width - 1 - ix : searchSize;
                    int i = iy * width + ix;
                    if (Float.compare(data[i], noData) == 0) {
                        float sum = 0;
                        float n = 0;
                        for (int iiy = iy - sr[0]; iiy <= iy + sr[1]; iiy++) {
                            for (int iix = ix - sr[2]; iix <= ix + sr[3]; iix++) {
                                int j = iiy * width + iix;
                                float w = 1f;
                                if (Float.compare(data[j], noData) != 0) {
                                    w = 1.5f * (searchSize * 2 + 1) * (searchSize * 2 + 1);
                                }
                                if (Float.compare(fillData1[j], noData) != 0) {
                                    sum += fillData1[j] * w;
                                    n += w;
                                }
                            }
                        }
                        fillData2[i] = (n > 0) ? sum / n : data[i];
                    } else {
                        fillData2[i] = data[i];
                    }
                    if (Float.compare(fillData2[i], noData) != 0) total_new += fillData2[i];
                }
                sr[2] = searchSize;
                sr[3] = searchSize;
            }
            float[] swap = fillData1;
            fillData1 = fillData2;
            fillData2 = swap;
        }
        return fillData1;
    }

    private Band createBandFromData(String bandName, int width, int height, float[] data) {
        SampleModel sm = RasterFactory.createBandedSampleModel(DataBuffer.TYPE_FLOAT, width, height, 1);
        ColorModel cm = PlanarImage.createColorModel(sm);
        WritableRaster wr = RasterFactory.createWritableRaster(sm, null);
        TiledImage ti = new TiledImage(0, 0, width, height, 0, 0, sm, cm);
        wr.setSamples(0, 0, width, height, 0, data);
        ti.setData(wr);
        Band b = new Band(bandName , ProductData.TYPE_FLOAT32, width, height);
        b.setSourceImage(ti);
        return b;
    }

    private void upscaleTileBilinear(Tile srcTile, Tile tarTile, Rectangle tarRec, ProgressMonitor pm) {

        final int tarX = tarRec.x;
        final int tarY = tarRec.y;
        final int tarWidth = tarRec.width;
        final int tarHeight = tarRec.height;

        for (int iTarY = tarY; iTarY < tarY + tarHeight; iTarY++) {
            int iSrcY = (iTarY - offset) / scale;
            if (iSrcY >= srcTile.getMaxY()) iSrcY = srcTile.getMaxY() - 1;
            float yFac = (float) (iTarY - offset) / scale - iSrcY;
            for (int iTarX = tarX; iTarX < tarX + tarWidth; iTarX++) {
                checkForCancellation(pm);
                int iSrcX = (iTarX - offset) / scale;
                if (iSrcX >= srcTile.getMaxX()) iSrcX = srcTile.getMaxX() - 1;
                float xFrac = (float) (iTarX - offset) / scale - iSrcX;
                float erg = 0;
                try{
                    erg = (1.0f - xFrac) * (1.0f - yFac) * srcTile.getSampleFloat(iSrcX, iSrcY);
                    erg +=        (xFrac) * (1.0f - yFac) * srcTile.getSampleFloat(iSrcX+1, iSrcY);
                    erg += (1.0f - xFrac) *        (yFac) * srcTile.getSampleFloat(iSrcX, iSrcY+1);
                    erg +=        (xFrac) *        (yFac) * srcTile.getSampleFloat(iSrcX+1, iSrcY+1);
                } catch (Exception ex) {
                    System.err.println(iTarX+" / "+iTarY);
                    System.err.println(ex.getMessage());
                }
                tarTile.setSample(iTarX, iTarY, erg);
            }
        }
    }

    private void upscaleTileCopy(Tile srcTile, Tile tarTile, Rectangle tarRec, ProgressMonitor pm) {

        final int tarX = tarRec.x;
        final int tarY = tarRec.y;
        final int tarWidth = tarRec.width;
        final int tarHeight = tarRec.height;

        for (int iTarY = tarY; iTarY < tarY + tarHeight; iTarY++) {
            int iSrcY = iTarY / scale;
            if (iSrcY >= srcTile.getHeight()) iSrcY = srcTile.getHeight() - 1;
            for (int iTarX = tarX; iTarX < tarX + tarWidth; iTarX++) {
                if (pm.isCanceled()) {
                    break;
                }
                int iSrcX = iTarX / scale;
                if (iSrcY >= srcTile.getWidth()) iSrcY = srcTile.getWidth() - 1;
                float erg = srcTile.getSampleFloat(iSrcX, iSrcY);
                tarTile.setSample(iTarX, iTarY, erg);
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
            super(UpSclOp.class);
        }
    }
}
