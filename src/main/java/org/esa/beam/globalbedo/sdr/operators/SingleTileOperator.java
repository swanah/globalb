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
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.ProductUtils;

import java.awt.*;

// TODO - Rename this operator
// TODO - Adapt OperatorMetadata
// TODO - Adapt javadoc

/**
 * The sample operator implementation for an algorithm
 * that can output bands  of the target product independently of each other.
 */
@OperatorMetadata(alias="SingleTileOp",
                  description = "Computes <X> and <Y> from <Z>",
                  authors = "Me, myself and I",
                  version = "1.0",
                  copyright = "(C) 2010 by Brockmann Consult GmbH (beam@brockmann-consult.de)")
public class SingleTileOperator extends Operator {

    @SourceProduct
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct;
    @Parameter
    private String sourceBandName;
    @Parameter
    private String targetBandName1;
    @Parameter
    private String targetBandName2;

    private Band sourceBand;
    private Band targetBand1;
    private Band targetBand2;

    /**
     * Default constructor. The graph processing framework
     * requires that an operator has a default constructor.
     */
    public SingleTileOperator() {
    }

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link org.esa.beam.framework.datamodel.Product} annotated with the
     * {@link org.esa.beam.framework.gpf.annotations.TargetProduct TargetProduct} annotation or
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
        targetProduct = new Product("gpftestop1",
                                    "org.esa.beam",
                                    sourceProduct.getSceneRasterWidth(),
                                    sourceProduct.getSceneRasterHeight());
        sourceBand = sourceProduct.getBand(sourceBandName);
        targetBand1 = targetProduct.addBand(targetBandName1, sourceBand.getDataType());
        targetBand2 = targetProduct.addBand(targetBandName2, sourceBand.getDataType());
        ProductUtils.copyGeoCoding(sourceProduct, targetProduct);
        // Some target products may require more aid from ProductUtils methods...
    }

    /**
     * Called by the framework in order to compute a tile for the given target band.
     * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
     *
     * @param targetBand The target band.
     * @param targetTile The current tile associated with the target band to be computed.
     * @param pm         A progress monitor which should be used to determine computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException
     *          If an error occurs during computation of the target raster.
     */
    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        Rectangle targetTileRectangle = targetTile.getRectangle();
        Tile sourceRaster = getSourceTile(sourceBand, targetTileRectangle, pm);

        // Some algorithms may require the geoCoding...
        // GeoCoding geoCoding = sourceProduct.getGeoCoding();

        int x0 = targetTileRectangle.x;
        int y0 = targetTileRectangle.y;
        int w = targetTileRectangle.width;
        int h = targetTileRectangle.height;
        for (int y = y0; y < y0 + h; y++) {
            for (int x = x0; x < x0 + w; x++) {
                // Some algorithms may require a geoPos for given pixelPos...
                // GeoPos geoPos = geoCoding.getGeoPos(new PixelPos(x, y), null);

                double v = sourceRaster.getSampleDouble(x, y);
                if (targetTile.getRasterDataNode() == targetBand1) {
                    double v1 = 0.1 * v; // Place your transformation maths here
                    targetTile.setSample(x, y, v1);
                } else if (targetTile.getRasterDataNode() == targetBand2) {
                    double v2 = 0.2 * v; // Place your transformation maths here
                    targetTile.setSample(x, y, v2);
                }
            }
        }
    }

    /**
     * The SPI is used to register this operator in the graph processing framework
     * via the SPI configuration file
     * {@code META-INF/services/org.esa.beam.framework.gpf.OperatorSpi}.
     * This class may also serve as a factory for new operator instances.
     * @see OperatorSpi#createOperator()
     * @see OperatorSpi#createOperator(java.util.Map, java.util.Map)
     */
    public static class Spi extends OperatorSpi {
        public Spi() {
            super(SingleTileOperator.class);
        }
    }
}
