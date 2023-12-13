# This file is part of meas_extensions_shapeHSM.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = (
    "HigherOrderMomentsPlugin",
    "HigherOrderMomentsSourcePlugin",
    "HigherOrderMomentsPSFPlugin",
)

import lsst.geom as geom
import lsst.meas.base as measBase
import numpy as np
from lsst.pex.config import Field, FieldValidationError, ListField


class HigherOrderMomentsConfig(measBase.SingleFramePluginConfig):
    min_order = Field[int](
        doc="Minimum order of the higher order moments to compute",
        default=4,
    )

    max_order = Field[int](
        doc="Maximum order of the higher order moments to compute",
        default=4,
    )

    def validate(self):
        if self.min_order > self.max_order:
            raise FieldValidationError(
                self.min_order, self, "min_order must be less than or equal to max_order"
            )
        super().validate()


class HigherOrderMomentsPlugin(measBase.SingleFramePlugin):
    """Base plugin for higher moments measurement"""

    ConfigClass = HigherOrderMomentsConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Define flags for possible issues that might arise during measurement.
        flagDefs = measBase.FlagDefinitionList()

        self.FAILURE = flagDefs.addFailureFlag("General failure flag, set if anything went wrong")

        self.NO_PIXELS = flagDefs.add("flag_no_pixels", "No pixels to measure")
        self.NOT_CONTAINED = flagDefs.add(
            "flag_not_contained", "Center not contained in footprint bounding box"
        )

        self.GALSIM = flagDefs.add("flag_galsim", "GalSim failure")

        self.NO_PSF = flagDefs.add("flag_no_psf", "Exposure lacks PSF")

        # Embed the flag definitions in the schema using a flag handler.
        self.flagHandler = measBase.FlagHandler.addFields(schema, name, flagDefs)

        self.pqlist = self._get_pq_full()

    @classmethod
    def getExecutionOrder(cls):
        return cls.FLUX_ORDER

    def fail(self, record, error=None):
        # Docstring inherited.
        self.flagHandler.handleFailure(record)

    def getAllNames(self):
        for p, q in self._get_pq_full():
            if p + q >= self.config.min_order:
                yield f"{p}{q}"

    def _get_pq_full(self):
        """Get a list of the orders to measure as a tuple.

        Returns
        -------
        pqlist: `list` [`tuples`]
            A list of tuples of the form (p, q) where p and q denote the order
            in x and y direction.
        """
        pq_list = []

        for n in range(self.config.min_order, self.config.max_order + 1):
            p = 0
            q = n

            pq_list.append((p, q))

            while p < n:
                p += 1
                q -= 1
                pq_list.append((p, q))

        return pq_list

    def _calculate_higher_order_moments(self, mi, center, M, badpix=None, setMaskedPixelsToZero=False):
        """
        image : `~lsst.afw.image.Image`
            Image of the PSF
        center: `~lsst.geom.Point2D`
            First order moments of ``image``. This is used as the peak of the
            Gaussian weight image.
        M : `~numpy.ndarray`
            A 2x2 numpy array representing the second order moments of
            ``image``. This is used to generate the Gaussian weight image.
        badpix : `~numpy.ndarray`
            A 2D array having the same shape and orientation as ``image.array``
            that denotes which pixels are bad and should not be accounted for
            when computing the moments.
        setMaskedPixelsToZero: `bool`
            Whether to treat pixels corresponding to ``badpix`` should be set
            to zero, or replaced by a scaled version of the weight image.
        """

        results_list = []

        image_array = mi.array

        y, x = np.mgrid[: image_array.shape[0], : image_array.shape[1]]

        inv_M = np.linalg.inv(M)

        evalues, evectors = np.linalg.eig(inv_M)
        # Ensuring square root matrix exists

        sqrt_inv_M = evectors * np.sqrt(evalues) @ np.linalg.inv(evectors)

        bbox = mi.getBBox()
        pos = np.array([x - (center.getX() - bbox.getMinX()), y - (center.getY() - bbox.getMinY())])

        std_pos = np.einsum("ij,jqp->iqp", sqrt_inv_M, pos)
        weight = np.exp(-0.5 * np.einsum("ijk,ijk->jk", std_pos, std_pos))

        std_x, std_y = std_pos

        if badpix is not None and badpix.any():
            if setMaskedPixelsToZero:
                image_array[badpix] = 0.0
            else:
                image_array[badpix] = weight[badpix] * image_array[~badpix].sum() / weight[~badpix].sum()

        image_weight = weight * image_array
        normalization = np.sum(image_weight)

        for p, q in self.pqlist:
            this_moment = np.sum(std_x**p * std_y**q * image_weight) / normalization
            results_list.append(this_moment)

        return results_list


class HigherOrderMomentsSourceConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the source
    """

    badMaskPlanes = ListField[str](
        doc="Mask planes used to reject bad pixels.",
        default=["BAD", "SAT"],
    )

    setMaskedPixelsToZero = Field[bool](
        doc="Set masked pixels to zero? If False, they are replaced by the "
        "scaled version of the adaptive weights.",
        default=False,
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsSource")
class HigherOrderMomentsSourcePlugin(HigherOrderMomentsPlugin):
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsSourceConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)

        # Add column names for the sources

        for suffix in self.getAllNames():
            schema.addField(
                schema.join(name, suffix),
                type=float,
                doc=f"Higher order moments M_{suffix} for source",
            )

    def measure(self, record, exposure):
        try:
            center = geom.Point2D(
                record["ext_shapeHSM_HsmSourceMoments_x"],
                record["ext_shapeHSM_HsmSourceMoments_y"],
            )
            M = np.zeros((2, 2))
            M[0, 0] = record["ext_shapeHSM_HsmSourceMoments_xx"]
            M[1, 1] = record["ext_shapeHSM_HsmSourceMoments_yy"]
            M[0, 1] = M[1, 0] = record["ext_shapeHSM_HsmSourceMoments_xy"]
        except KeyError:
            raise measBase.FatalAlgorithmError("HSM moments algorithm was not run.")

        # get the bounding box of the source footprint
        bbox = record.getFootprint().getBBox()

        # self.logger.info(bbox)

        # Check that the bounding box has non-zero area.
        if bbox.getArea() == 0:
            raise measBase.MeasurementError(self.NO_PIXELS.doc, self.NO_PIXELS.number)

        # # Ensure that the centroid is within the bounding box.
        # if not bbox.contains(Point2I(center)):
        #     raise measBase.MeasurementError(self.NOT_CONTAINED.doc, self.NOT_CONTAINED.number)

        badpix = exposure.mask[bbox].array.copy()
        bitValue = exposure.mask.getPlaneBitMask(self.config.badMaskPlanes)
        badpix &= bitValue

        # imageArray[badpix.astype(bool)] = 0.0
        # np.savetxt('savefig.txt', imageArray)
        # self.logger.info('array' + str(imageArray))

        # get Galsim image from the image array

        # Measure all the moments together to save time
        try:
            hm_measurement = self._calculate_higher_order_moments(
                exposure.image[bbox],
                center,
                M,
                badpix.astype(bool),
                setMaskedPixelsToZero=self.config.setMaskedPixelsToZero,
            )
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments
        for i in range(len(hm_measurement)):
            (p, q) = self.pqlist[i]
            M_pq = hm_measurement[i]
            this_column_name = self.name + f"_{p}{q}"

            record.set(this_column_name, M_pq)


class HigherOrderMomentsPSFConfig(HigherOrderMomentsConfig):

    """
    Configuration for the higher order moments of the PSF
    """

    useSourceCentroidOffset = Field[bool](
        doc="Use source centroid offset?",
        default=False,
    )


@measBase.register("ext_shapeHSM_HigherOrderMomentsPSF")
class HigherOrderMomentsPSFPlugin(HigherOrderMomentsPlugin):
    """
    Plugin class for Higher Order Moments measurement of the source
    """

    ConfigClass = HigherOrderMomentsPSFConfig

    def __init__(self, config, name, schema, metadata, logName=None):
        super().__init__(config, name, schema, metadata, logName=logName)
        self.centroidExtractor = measBase.SafeCentroidExtractor(schema, name)

        # Add column names for the PSFs
        for suffix in self.getAllNames():
            schema.addField(
                schema.join(name, suffix),
                type=float,
                doc=f"Higher order moments M_{suffix} for PSF",
            )

    def measure(self, record, exposure):
        M = np.zeros((2, 2))
        try:
            M[0, 0] = record["ext_shapeHSM_HsmPsfMoments_xx"]
            M[1, 1] = record["ext_shapeHSM_HsmPsfMoments_yy"]
            M[0, 1] = M[1, 0] = record["ext_shapeHSM_HsmPsfMoments_xy"]
            center = geom.Point2D(
                record["ext_shapeHSM_HsmPsfMoments_x"],
                record["ext_shapeHSM_HsmPsfMoments_y"],
            )
            # center = self.centroidExtractor(record, self.flagHandler)

        except KeyError:
            raise measBase.FatalAlgorithmError("HSM PSF moments algorithm was not run.")

        # get the psf and psf image from the exposure
        psf = exposure.getPsf()

        # check if the psf is none:

        if not psf:
            raise measBase.FatalAlgorithmError(self.NO_PSF.doc, self.NO_PSF.number)

        centroid = self.centroidExtractor(record, self.flagHandler)
        if self.config.useSourceCentroidOffset:
            # 1. Using `computeImage()` returns an image in the same coordinate
            # system as the pixelized image.
            psfImage = psf.computeImage(centroid)
            # centroid = center
            centroid.x -= center.x
            centroid.y -= center.y
        else:
            psfImage = psf.computeKernelImage(center)
            # 2. Using `computeKernelImage()` to return an image does not
            # retain any information about the original bounding box of the
            # PSF. We therefore reset the origin to be the same as the
            # pixelized image.
            center0 = geom.Point2I(center)
            xy0 = geom.Point2I(center0.x + psfImage.getX0(), center0.y + psfImage.getY0())
            psfImage.setXY0(xy0)
            psfBBox = psfImage.getBBox()
            centroid = geom.Point2D(psfBBox.getMin() + psfBBox.getDimensions() // 2)

        # Measure all the moments together to save time
        try:
            hm_measurement = self._calculate_higher_order_moments(psfImage, centroid, M)
        except Exception as e:
            raise measBase.MeasurementError(e)

        # Record the moments
        for i in range(len(hm_measurement)):
            (p, q) = self.pqlist[i]
            M_pq = hm_measurement[i]

            this_column_name = self.name + f"_{p}{q}"
            record.set(this_column_name, M_pq)
