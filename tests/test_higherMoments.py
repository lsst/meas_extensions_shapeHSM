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
    """A test case for shape measurement"""

    def setUp(self):
        super().setUp()
        """Create an exposure and run measurement on the source and the PSF"""

        # Initialize a config and activate the plugin
        sfmConfig = SingleFrameMeasurementConfig()
        sfmConfig.plugins.names |= [
            "ext_shapeHSM_HsmSourceMoments",
            "ext_shapeHSM_HsmPsfMoments",
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        ]
        for plugin_name in (
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        ):
            sfmConfig.plugins[plugin_name].max_order = 6
            sfmConfig.plugins[plugin_name].min_order = 0
        # sfmConfig.plugins['ext_shapeHSM_HigherOrderMomentsSource'].max_order = 6

        # Create a minimal schema (columns)
        self.schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

        # Create a task
        sfmTask = SingleFrameMeasurementTask(config=sfmConfig, schema=self.schema)

        # schema.getAliasMap().set("slot_Centroid", "ext_shapeHSM_HsmSourceMoments")
        # schema.getAliasMap().set("slot_Shape", "ext_shapeHSM_HsmSourceMoments")
        # schema.getAliasMap().set("slot_PsfShape", "ext_shapeHSM_HsmPsfMoments")

        dataset = self.createDataset()

        # Get the exposure and catalog.
        exposure, catalog = dataset.realize(0.0, sfmTask.schema, randomSeed=0)

        self.sfmConfig = sfmConfig
        self.catalog = catalog
        self.exposure = exposure
        self.task = sfmTask

        self.addMaskBits()

        # sfmTask.run(catalog, exposure)

        # self.catalog = catalog.asAstropy()
        return self.catalog

    @staticmethod
    def addMaskBits():
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
        # return self.catalog

    def check_odd_moments(self, row, plugin_name, atol):
        for n in (3, 5):
            for p in range(n + 1):
                with self.subTest((p, n - p)):
                    self.assertFloatsAlmostEqual(row[f"{plugin_name}_{p}{n-p}"], 0.0, atol=atol)

    def check_even_moments(self, row, plugin_name, atol):
        M_source_40 = row[f"{plugin_name}_40"]
        M_source_31 = row[f"{plugin_name}_31"]
        M_source_22 = row[f"{plugin_name}_22"]
        M_source_13 = row[f"{plugin_name}_13"]
        M_source_04 = row[f"{plugin_name}_04"]

        M_source_50 = row[f"{plugin_name}_50"]
        M_source_41 = row[f"{plugin_name}_41"]
        M_source_32 = row[f"{plugin_name}_32"]
        M_source_23 = row[f"{plugin_name}_23"]
        M_source_14 = row[f"{plugin_name}_14"]
        M_source_05 = row[f"{plugin_name}_05"]

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

        self.assertFloatsAlmostEqual(M_source_50, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_41, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_32, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_23, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_14, 0.0, atol=atol)
        self.assertFloatsAlmostEqual(M_source_05, 0.0, atol=atol)

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

        # catalog1 = self.runMeasurement(plugin_name)
        catalog1 = self.catalog

        for row in catalog1:
            self.check(row, plugin_name, atol=MOMENTS_DECIMAL)

            # M_psf_40 = row["ext_shapeHSM_HigherOrderMomentsPSF_40"]
            # M_psf_31 = row["ext_shapeHSM_HigherOrderMomentsPSF_31"]
            # M_psf_22 = row["ext_shapeHSM_HigherOrderMomentsPSF_22"]
            # M_psf_13 = row["ext_shapeHSM_HigherOrderMomentsPSF_13"]
            # M_psf_04 = row["ext_shapeHSM_HigherOrderMomentsPSF_04"]

            # M_psf_30 = row["ext_shapeHSM_HigherOrderMomentsPSF_30"]
            # M_psf_21 = row["ext_shapeHSM_HigherOrderMomentsPSF_21"]
            # M_psf_12 = row["ext_shapeHSM_HigherOrderMomentsPSF_12"]
            # M_psf_03 = row["ext_shapeHSM_HigherOrderMomentsPSF_03"]

            # M_psf_50 = row["ext_shapeHSM_HigherOrderMomentsPSF_50"]
            # M_psf_41 = row["ext_shapeHSM_HigherOrderMomentsPSF_41"]
            # M_psf_32 = row["ext_shapeHSM_HigherOrderMomentsPSF_32"]
            # M_psf_23 = row["ext_shapeHSM_HigherOrderMomentsPSF_23"]
            # M_psf_14 = row["ext_shapeHSM_HigherOrderMomentsPSF_14"]
            # M_psf_05 = row["ext_shapeHSM_HigherOrderMomentsPSF_05"]

            # M_psf_60 = row["ext_shapeHSM_HigherOrderMomentsPSF_60"]
            # M_psf_51 = row["ext_shapeHSM_HigherOrderMomentsPSF_51"]
            # M_psf_42 = row["ext_shapeHSM_HigherOrderMomentsPSF_42"]
            # M_psf_33 = row["ext_shapeHSM_HigherOrderMomentsPSF_33"]
            # M_psf_24 = row["ext_shapeHSM_HigherOrderMomentsPSF_24"]
            # M_psf_15 = row["ext_shapeHSM_HigherOrderMomentsPSF_15"]
            # M_psf_06 = row["ext_shapeHSM_HigherOrderMomentsPSF_06"]

            # self.assertAlmostEqual(M_psf_40, 0.75, MOMENTS_DECIMAL)
            # self.assertAlmostEqual(M_psf_31, 0.0, MOMENTS_DECIMAL)
            # self.assertAlmostEqual(M_psf_22, 0.25, MOMENTS_DECIMAL)
            # self.assertAlmostEqual(M_psf_13, 0.0, MOMENTS_DECIMAL)
            # self.assertAlmostEqual(M_psf_04, 0.75, MOMENTS_DECIMAL)
            # self.assertFloatsAlmostEqual(M_psf_30, 0.0,atol= 2e-3)
            # self.assertFloatsAlmostEqual(M_psf_21, 0.0,atol= 2e-3)
            # self.assertFloatsAlmostEqual(M_psf_12, 0.0,atol= 2e-3)
            # self.assertFloatsAlmostEqual(M_psf_03, 0.0,atol= 2e-3)

            # self.assertFloatsAlmostEqual(M_psf_50, 0.0, atol= 8e-3)
            # self.assertFloatsAlmostEqual(M_psf_41, 0.0, atol= 8e-3)
            # self.assertFloatsAlmostEqual(M_psf_32, 0.0, atol= 8e-3)
            # self.assertFloatsAlmostEqual(M_psf_23, 0.0, atol= 8e-3)
            # self.assertFloatsAlmostEqual(M_psf_14, 0.0, atol= 8e-3)
            # self.assertFloatsAlmostEqual(M_psf_05, 0.0, atol= 8e-3)

            # self.assertFloatsAlmostEqual(M_psf_60, 1.875, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_51, 0.0, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_42, 0.375, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_33, 0.0, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_24, 0.375, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_15, 0.0, atol= 5e-3)
            # self.assertFloatsAlmostEqual(M_psf_06, 1.875, atol= 5e-3)

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
        atol = 2e-2 if useSourceCentroidOffset else MOMENTS_DECIMAL
        for row in catalog1:
            with self.subTest():
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
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsPSF"].min_order = 0
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsPSF"].max_order = 2

        self.runMeasurement()

        catalog1 = self.catalog
        atol = 2e-2 if useSourceCentroidOffset else MOMENTS_DECIMAL
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
            "ext_shapeHSM_HsmSourceMomentsRound",
            "base_SdssShape",
        )
    )
    def test_consistent_weight(self, target_plugin_name):
        # Pause the execution of the measurement task before the higher order
        # moments plugins.

        if False:
            self.task.config.slots.centroid = "ext_shapeHSM_HsmSourceMoments"
            self.task.config.slots.shape = "ext_shapeHSM_HsmSourceMoments"
            self.task.config.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
        elif False:
            self.task.config.slots.centroid = "base_SdssShape"
            self.task.config.slots.shape = "base_SdssShape"
            self.task.config.slots.psfShape = "base_SdssShape"

        pause_order = self.task.plugins["ext_shapeHSM_HigherOrderMomentsSource"].getExecutionOrder()
        self.runMeasurement(endOrder=pause_order)

        for suffix in (
            "x",
            "y",
            "xx",
            "yy",
            "xy",
        ):
            self.catalog[f"ext_shapeHSM_HsmSourceMoments_{suffix}"] = self.catalog[f"base_SdssShape_{suffix}"]

        for suffix in (
            "xx",
            "yy",
            "xy",
        ):
            self.catalog[f"ext_shapeHSM_HsmPsfMoments_{suffix}"] = self.catalog[
                f"base_SdssShape_psf_{suffix}"
            ]

        # Resume the execution of the measurement task.
        self.runMeasurement(beginOrder=pause_order)

        # catalog1 = self.runMeasurement(plugin_name)
        catalog1 = self.catalog
        atol = 6e-4
        for plugin_name in (
            "ext_shapeHSM_HigherOrderMomentsSource",
            "ext_shapeHSM_HigherOrderMomentsPSF",
        ):
            for i, row in enumerate(catalog1):
                with self.subTest((plugin_name, i)):
                    self.check(row, plugin_name, atol=atol)


class HigherMomentTestCaseWithMask(HigherMomentsBaseTestCase):
    def addMaskBits(self):
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
        """Test that we can instantiate and play with a measureShape"""

        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].min_order = 0
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].max_order = 2
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
        # GalSim does not set masked pixels to zero.
        # So we set them to zero as well for the comparison.
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = True
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].min_order = 4
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].max_order = 4

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
                # Check that the rho4 moments are non-trivial and differ from 2.0

        self.assertTrue((np.array(delta_rho4s) > 10 * MOMENTS_DECIMAL).any())

    def testHsmSourceHigherMoments(self, plugin_name="ext_shapeHSM_HigherOrderMomentsSource"):
        """Test that we can instantiate and play with a measureShape"""

        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].badMaskPlanes = ["BAD", "SAT"]
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = False

        self.runMeasurement()

        # catalog1 = self.runMeasurement(plugin_name)
        catalog1 = self.catalog

        for row in catalog1:
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_00"], 1.0, atol=3e-1)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_01"], 0.0, atol=3e-1)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_10"], 0.0, atol=3e-1)

            self.assertFloatsAlmostEqual(row[f"{plugin_name}_20"], 0.5, atol=3e-1)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_11"], 0.0, atol=3e-1)
            self.assertFloatsAlmostEqual(row[f"{plugin_name}_02"], 0.5, atol=3e-1)

            self.check(row, plugin_name, atol=3e-1)


class NonGaussianTestCase(HigherMomentTestCaseWithMask):
    @staticmethod
    def createDataset():
        # Create a simple, fake dataset with centroids at integer or
        # half-integer positions.
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
        dataset = lsst.meas.base.tests.TestDataset(bbox)
        # Create a point source with Gaussian PSF
        dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

        # Create a galaxy with Gaussian PSF
        dataset.addSource(300000.0, lsst.geom.Point2D(76, 79), lsst.afw.geom.Quadrupole(2.0, 3.0, 0.5))
        return dataset

    def addMaskBits(self):
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
        """Test that we can instantiate and play with a measureShape"""

        self.runMeasurement()

        # catalog1 = self.runMeasurement(plugin_name)
        catalog1 = self.catalog

        for row in catalog1:
            self.check_odd_moments(row, plugin_name, atol=MOMENTS_DECIMAL)
            self.check_even_moments(row, plugin_name, atol=3e-1)

    def testAgainstPiff(self, plugin_name="ext_shapeHSM_HigherOrderMomentsSource"):
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].badMaskPlanes = ["BAD", "SAT"]
        self.task.config.plugins["ext_shapeHSM_HigherOrderMomentsSource"].setMaskedPixelsToZero = False

        self.runMeasurement()
        catalog1 = self.catalog


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
