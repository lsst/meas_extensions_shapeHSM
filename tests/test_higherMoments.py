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

import unittest

import galsim
import lsst.afw.geom
import lsst.meas.base.tests
import lsst.meas.extensions.shapeHSM  # noqa: F401
import lsst.utils.tests as tests
import numpy as np
from lsst.meas.base import SingleFrameMeasurementConfig, SingleFrameMeasurementTask

MOMENTS_DECIMAL = 1e-5  # 3e-1 # Number of decimals for equality in moments


class HigherMomentsBaseTestCase(tests.TestCase):
    """Base test case to test higher order moments."""

    def setUp(self):
        """Create an exposure and run measurement on the source and the PSF"""
        super().setUp()

        # Initialize a config and activate the plugin
        sfmConfig = SingleFrameMeasurementConfig()
        sfmConfig.plugins.names |= [
            "ext_shapeHSM_HsmSourceMoments",
            "ext_shapeHSM_HsmSourceMomentsRound",
            "ext_shapeHSM_HsmPsfMoments",
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        ]
        # The min and max order determine the schema and cannot be changed
        # after the Task is created. So we set it generously here.
        for plugin_name in (
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        ):
            sfmConfig.plugins[plugin_name].max_order = 7
            sfmConfig.plugins[plugin_name].min_order = 0

        # Create a minimal schema (columns)
        self.schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

        # Create a task
        sfmTask = SingleFrameMeasurementTask(config=sfmConfig, schema=self.schema)

        dataset = self.createDataset()

        # Get the exposure and catalog.
        exposure, catalog = dataset.realize(0.0, sfmTask.schema, randomSeed=0)

        self.catalog = catalog
        self.exposure = exposure
        self.task = sfmTask

        self.addMaskBits()

        return self.catalog

    @staticmethod
    def addMaskBits():
        """Add mask bits to the exposure.

        This must go along with the createDataset method. This is a no-op for
        the base class and subclasses must set mask bits depending on the test.
        """
        pass

    @staticmethod
    def createDataset():
        # Create a simple, fake dataset
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        # Create a point source with Gaussian PSF
        dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

        # Create a galaxy with Gaussian PSF
        dataset.addSource(300000.0, lsst.geom.Point2D(76.3, 79.2), lsst.afw.geom.Quadrupole(2.0, 3.0, 0.5))
        return dataset

    def runMeasurement(self, **kwargs):
        """Run measurement on the source and the PSF"""
        self.task.run(self.catalog, self.exposure, **kwargs)

    def check_odd_moments(self, row, plugin_name, atol, orders=(3, 5)):
        for n in orders:
            for p in range(n + 1):
                with self.subTest((p, n - p)):
                    self.assertFloatsAlmostEqual(row[f"{plugin_name}_{p}{n-p}"], 0.0, atol=atol)

    def check_even_moments(self, row, plugin_name, atol):
        M_source_40 = row[f"{plugin_name}_40"]
        M_source_31 = row[f"{plugin_name}_31"]
        M_source_22 = row[f"{plugin_name}_22"]
        M_source_13 = row[f"{plugin_name}_13"]
        M_source_04 = row[f"{plugin_name}_04"]

        M_source_60 = row[f"{plugin_name}_60"]
        M_source_51 = row[f"{plugin_name}_51"]
        M_source_42 = row[f"{plugin_name}_42"]
        M_source_33 = row[f"{plugin_name}_33"]
        M_source_24 = row[f"{plugin_name}_24"]
        M_source_15 = row[f"{plugin_name}_15"]
        M_source_06 = row[f"{plugin_name}_06"]

        self.assertFloatsAlmostEqual(M_source_40, 0.75, atol=atol)
        self.assertFloatsAlmostEqual(M_source_31, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_22, 0.25, atol=atol)
        self.assertFloatsAlmostEqual(M_source_13, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_04, 0.75, atol=atol)

        self.assertFloatsAlmostEqual(M_source_60, 1.875, atol=atol)
        self.assertFloatsAlmostEqual(M_source_51, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_42, 0.375, atol=atol)
        self.assertFloatsAlmostEqual(M_source_33, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_24, 0.375, atol=atol)
        self.assertFloatsAlmostEqual(M_source_15, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_06, 1.875, atol=atol)

    def check(self, row, plugin_name, atol):
        self.check_odd_moments(row, plugin_name, atol)
        self.check_even_moments(row, plugin_name, atol)


class HigherOrderMomentsTestCase(HigherMomentsBaseTestCase):
    @lsst.utils.tests.methodParameters(
        plugin_name=(
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        )
    )
    def testHsmSourceMoments(self, plugin_name):
        """Test that we can instantiate and play with a measureShape"""

        self.runMeasurement()

        catalog1 = self.catalog

        for row in catalog1:
            self.check(row, plugin_name, atol=MOMENTS_DECIMAL)

    @lsst.utils.tests.methodParameters(useSourceCentroidOffset=(False, True))
    def test_hsm_psf_higher_moments(self, useSourceCentroidOffset):
        """Test that we can instantiate and play with a measureShape"""

        self.task.config.plugins[
            "ext_shapeHSM_HsmPsfMoments"
        ].useSourceCentroidOffset = useSourceCentroidOffset
        self.task.config.plugins[
            "ext_shapeHSM_HigherOrderMomentsPSF"
        ].useSourceCentroidOffset = useSourceCentroidOffset

        self.runMeasurement()

        catalog1 = self.catalog

        # useSourceCentroidOffset = False results in more accurate results.
        # Adjust the absolute tolerance accordingly.
        atol = 9e-3 if useSourceCentroidOffset else MOMENTS_DECIMAL

        for i, row in enumerate(catalog1):
            with self.subTest(i=i):
                self.check(row, "ext_shapeHSM_HigherOrderMomentsPSF", atol=atol)

    @lsst.utils.tests.methodParameters(useSourceCentroidOffset=(True, False))
    def test_hsm_psf_lower_moments(self, useSourceCentroidOffset):
        """Test that we can instantiate and play with a measureShape"""
        plugin_name = "ext_shapeHSM_HigherOrderMomentsPSF"
        self.task.config.plugins[
            "ext_shapeHSM_HsmPsfMoments"
        ].useSourceCentroidOffset = useSourceCentroidOffset
        self.task.config.plugins[
            "ext_shapeHSM_HigherOrderMomentsPSF"
        ].useSourceCentroidOffset = useSourceCentroidOffset

        self.runMeasurement()

        catalog1 = self.catalog

        # useSourceCentroidOffset = False results in more accurate results.
        # Adjust the absolute tolerance accordingly.
        atol = 5e-4 if useSourceCentroidOffset else MOMENTS_DECIMAL

        for i, row in enumerate(catalog1):
            with self.subTest(i=i):
                self.assertFloatsAlmostEqual(row[f"{plugin_name}_00"], 1.0, atol=atol)

                self.assertFloatsAlmostEqual(row[f"{plugin_name}_01"], 0.0, atol=atol)
                self.assertFloatsAlmostEqual(row[f"{plugin_name}_10"], 0.0, atol=atol)

                self.assertFloatsAlmostEqual(row[f"{plugin_name}_20"], 0.5, atol=atol)
                self.assertFloatsAlmostEqual(row[f"{plugin_name}_11"], 0.0, atol=atol)
                self.assertFloatsAlmostEqual(row[f"{plugin_name}_02"], 0.5, atol=atol)

    @lsst.utils.tests.methodParameters(
        target_plugin_name=(
            "truth",
            "base_SdssShape",
            "ext_shapeHSM_HsmSourceMomentsRound",
        )
    )
    def test_source_consistent_weight(self, target_plugin_name):
        """Test that when we get expected results when use a different set of
        consistent weights to measure the higher order moments of sources.
        """
        # Pause the execution of the measurement task before the higher order
        # moments plugins.

        pause_order = self.task.plugins["ext_shapeHSM_HigherOrderMomentsSource"].getExecutionOrder()
        self.runMeasurement(endOrder=pause_order)

        for suffix in (
            "x",
            "y",
            "xx",
            "yy",
            "xy",
        ):
            self.catalog[f"ext_shapeHSM_HsmSourceMoments_{suffix}"] = self.catalog[
                f"{target_plugin_name}_{suffix}"
            ]

        # Resume the execution of the measurement task.
        self.runMeasurement(beginOrder=pause_order)

        catalog1 = self.catalog
        atol = 1.2e-1 if target_plugin_name == "ext_shapeHSM_HsmSourceMomentsRound" else 6e-4
        plugin_name = "ext_shapeHSM_HigherOrderMomentsSource"

        for i, row in enumerate(catalog1):
            with self.subTest((plugin_name, i)):
                self.check(row, plugin_name, atol=atol)

    @lsst.utils.tests.methodParametersProduct(
        target_plugin_name=(
            "truth",
            "base_SdssShape_psf",
        ),
        useSourceCentroidOffset=(False, True),
    )
    def test_psf_consistent_weight(self, target_plugin_name, useSourceCentroidOffset):
        """Test that when we get expected results when use a different set of
        consistent weights to measure the higher order moments of PSFs.
        """
        self.task.config.plugins[
            "ext_shapeHSM_HigherOrderMomentsPSF"
        ].useSourceCentroidOffset = useSourceCentroidOffset

        # Pause the execution of the measurement task before the higher order
        # moments plugins.
        pause_order = self.task.plugins["ext_shapeHSM_HigherOrderMomentsPSF"].getExecutionOrder()
        self.runMeasurement(endOrder=pause_order)

        # Create a dictionary of PSF moments corresponding to the truth.
        # These are hardcoded in dataset.realize.
        truth_psf = {"xx": 4.0, "yy": 4.0, "xy": 0.0}

        for suffix in (
            "xx",
            "yy",
            "xy",
        ):
            if target_plugin_name == "truth":
                self.catalog[f"ext_shapeHSM_HsmPsfMoments_{suffix}"] = truth_psf[suffix]
            else:
                self.catalog[f"ext_shapeHSM_HsmPsfMoments_{suffix}"] = self.catalog[
                    f"{target_plugin_name}_{suffix}"
                ]

        # Resume the execution of the measurement task.
        self.runMeasurement(beginOrder=pause_order)

        catalog1 = self.catalog
        # useSourceCentroidOffset = False results in more accurate results.
        # Adjust the absolute tolerance accordingly.
        atol = 1.2e-2 if useSourceCentroidOffset else 8e-6
        plugin_name = "ext_shapeHSM_HigherOrderMomentsPSF"

        for i, row in enumerate(catalog1):
            with self.subTest((plugin_name, i)):
                self.check(row, plugin_name, atol=atol)


class HigherMomentTestCaseWithMask(HigherMomentsBaseTestCase):
    """A test case to measure higher order moments in the presence of masks.

    The tests serve checking the validity the algorithm on non-Gaussian
    profiles.
    """

    def addMaskBits(self):
        # Docstring inherited.
        for position in (
            lsst.geom.Point2I(48, 47),
            lsst.geom.Point2I(76, 79),
        ):
            self.exposure.mask[position] |= self.exposure.mask.getPlaneBitMask("BAD")
        for position in (
            lsst.geom.Point2D(49, 49),
            lsst.geom.Point2D(76, 79),
        ):
            self.exposure.mask[position] |= self.exposure.mask.getPlaneBitMask("SAT")

    def test_lower_order_moments(self, plugin_name="ext_shapeHSM_HigherOrderMomentsSource"):
        """Test that the lower order moments (2nd order or lower) is consistent
        even in the presence of masks.
        """
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = True

        self.runMeasurement()

        # catalog1 = self.runMeasurement(plugin_name)
        catalog1 = self.catalog

        for row in catalog1:
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_00"], 1.0, atol=MOMENTS_DECIMAL)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_01"], 0.0, atol=MOMENTS_DECIMAL)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_10"], 0.0, atol=MOMENTS_DECIMAL)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_20"], 0.5, atol=MOMENTS_DECIMAL)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_11"], 0.0, atol=MOMENTS_DECIMAL)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_02"], 0.5, atol=MOMENTS_DECIMAL)

        # for row in catalog1:
        #    self.check(row, plugin_name, atol=1.7e-1)

    def test_kurtosis(self):
        """Test the the kurtosis measurement against GalSim HSM implementation."""
        # GalSim does not set masked pixels to zero.
        # So we set them to zero as well for the comparison.
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = True

        self.runMeasurement()
        catalog1 = self.catalog

        delta_rho4s = []
        for i, row in enumerate(catalog1):
            bbox = row.getFootprint().getBBox()
            im = galsim.Image(self.exposure[bbox].image.array)
            badpix = self.exposure.mask[bbox].array.copy()
            bitValue = self.exposure.mask.getPlaneBitMask(["BAD", "SAT"])
            badpix &= bitValue
            badpix = galsim.Image(badpix, copy=False)
            shape = galsim.hsm.FindAdaptiveMom(im, badpix=badpix, strict=False)
            # r^4 = (x^2+y^2)^2 = x^4 + y^4 + 2x^2y^2
            rho4 = sum(
                (
                    row["ext_shapeHSM_HigherOrderMomentsSource_40"],
                    row["ext_shapeHSM_HigherOrderMomentsSource_04"],
                    row["ext_shapeHSM_HigherOrderMomentsSource_22"] * 2,
                )
            )
            delta_rho4s.append(abs(rho4 - 2.0))
            with self.subTest(i=i):
                self.assertFloatsAlmostEqual(shape.moments_rho4, rho4, atol=4e-7)

        # Check that at least one rho4 moment is non-trivial and differs from
        # the fiducial value of 2, by an amount much larger than the precision.
        self.assertTrue((np.array(delta_rho4s) > 1e-2).any(), "Unit test is too weak.")

    def testHsmSourceHigherMoments(self, plugin_name="ext_shapeHSM_HigherOrderMomentsSource"):
        """Test that we can instantiate and play with a measureShape"""

        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].badMaskPlanes = ["BAD", "SAT"]
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = False

        self.runMeasurement()

        catalog1 = self.catalog
        atol = 3e-1
        for row in catalog1:
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_00"], 1.0, atol=atol)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_01"], 0.0, atol=atol)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_10"], 0.0, atol=atol)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_20"], 0.5, atol=atol)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_11"], 0.0, atol=atol)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_02"], 0.5, atol=atol)

            self.check(row, plugin_name, atol=atol)


class HigherMomentTestCaseWithSymmetricMask(HigherMomentTestCaseWithMask):
    @staticmethod
    def createDataset():
        # Create a simple, fake dataset with centroids at integer or
        # half-integer positions to have a definite symmetry.
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        # Create a point source with Gaussian PSF
        dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

        # Create a galaxy with Gaussian PSF
        dataset.addSource(300000.0, lsst.geom.Point2D(76, 79), lsst.afw.geom.Quadrupole(2.0, 3.0, 0.5))
        return dataset

    def addMaskBits(self):
        # Docstring inherited.
        for position in (
            lsst.geom.Point2I(48, 48),
            lsst.geom.Point2I(73, 79),
        ):
            self.exposure.mask[position] |= self.exposure.mask.getPlaneBitMask("BAD")
        for position in (
            lsst.geom.Point2D(51, 51),
            lsst.geom.Point2D(79, 79),
        ):
            self.exposure.mask[position] |= self.exposure.mask.getPlaneBitMask("SAT")

    @lsst.utils.tests.methodParameters(plugin_name=("ext_shapeHSM_HigherOrderMomentsSource",))
    def testOddMoments(self, plugin_name):
        """Test that the odd order moments are close to expect values."""

        self.runMeasurement()

        catalog1 = self.catalog

        for row in catalog1:
            self.check_odd_moments(row, plugin_name, atol=MOMENTS_DECIMAL)
            self.check_even_moments(row, plugin_name, atol=3e-1)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
