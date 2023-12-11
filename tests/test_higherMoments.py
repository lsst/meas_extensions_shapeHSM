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


import sys
sys.path = ['/pscratch/sd/z/ztq1996/DM/meas_extensions_shapeHSM/python'] + sys.path

from lsst.meas.base import SingleFrameMeasurementConfig, SingleFrameMeasurementTask
import lsst.meas.extensions.shapeHSM as shapeHSM
import lsst.meas.base.tests
import lsst.afw.geom

import os
import numpy as np
import unittest
import lsst.utils.tests as tests
import itertools


MOMENTS_DECIMAL = 3  # Number of decimals for equality in moments



class HigherMomentsTestCase(tests.TestCase):
    """A test case for shape measurement"""


    def runMeasurement(self, testname, plugin_name):
        """Create an exposure and run measurement on the source and the PSF"""

        if testname == "point_source+gaussian_PSF":
            # Initialize a config and activate the plugin
            sfmConfig = SingleFrameMeasurementConfig()
            sfmConfig.plugins.names |= ["ext_shapeHSM_HsmSourceMoments",
                                        "ext_shapeHSM_HsmPsfMoments",
                plugin_name]
            sfmConfig.plugins[plugin_name].max_order = 6
            sfmConfig.plugins[plugin_name].min_order = 0
            # sfmConfig.plugins['ext_shapeHSM_HigherOrderMomentsSource'].max_order = 6

            if True:
                sfmConfig.slots.centroid = "ext_shapeHSM_HsmSourceMoments"
                sfmConfig.slots.shape = "ext_shapeHSM_HsmSourceMoments"
                sfmConfig.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
            else:
                sfmConfig.slots.centroid = "base_SdssShape"
                sfmConfig.slots.shape = "base_SdssShape"
                sfmConfig.slots.psfShape = "base_SdssShape"

            # Create a minimal schema (columns)
            schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

            # Create a task
            sfmTask = SingleFrameMeasurementTask(config=sfmConfig, schema=schema)

            #schema.getAliasMap().set("slot_Centroid", "ext_shapeHSM_HsmSourceMoments")
            # schema.getAliasMap().set("slot_Shape", "ext_shapeHSM_HsmSourceMoments")
            # schema.getAliasMap().set("slot_PsfShape", "ext_shapeHSM_HsmPsfMoments")

            # Create a simple, fake dataset
            bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(100, 100))
            dataset = lsst.meas.base.tests.TestDataset(bbox)
            # Create a point source with Gaussian PSF
            dataset.addSource(100000.0, lsst.geom.Point2D(49.5, 49.5))

            # Create a galaxy with Gaussian PSF
            dataset.addSource(300000.0, lsst.geom.Point2D(76.3, 79.2), lsst.afw.geom.Quadrupole(2.0, 3.0, 0.5))

            # Get the exposure and catalog.
            exposure, catalog = dataset.realize(0.0, sfmTask.schema, randomSeed=0)

        sfmTask.run(catalog, exposure)


        return catalog.asAstropy()

    @lsst.utils.tests.methodParameters(plugin_name=("ext_shapeHSM_HigherOrderMomentsSource", "ext_shapeHSM_HigherOrderMomentsPSF",))
    def testHsmSourceMoments(self, plugin_name):
        """Test that we can instantiate and play with a measureShape"""

        nFail = 0
        msg = ""

        for i in range(1):
            test1 = "point_source+gaussian_PSF"
            catalog1 = self.runMeasurement(test1, plugin_name)

        for row in catalog1:

            M_source_30 = row[f"{plugin_name}_30"]
            M_source_21 = row[f"{plugin_name}_21"]
            M_source_12 = row[f"{plugin_name}_12"]
            M_source_03 = row[f"{plugin_name}_03"]

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


            self.assertAlmostEqual(M_source_30, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_21, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_12, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_03, 0.0, MOMENTS_DECIMAL)


            self.assertAlmostEqual(M_source_40, 0.75, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_31, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_22, 0.25, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_13, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_04, 0.75, MOMENTS_DECIMAL)

            self.assertAlmostEqual(M_source_50, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_41, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_32, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_23, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_14, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_05, 0.0, MOMENTS_DECIMAL)

            self.assertAlmostEqual(M_source_60, 1.875, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_51, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_42, 0.375, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_33, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_24, 0.375, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_15, 0.0, MOMENTS_DECIMAL)
            self.assertAlmostEqual(M_source_06, 1.875, MOMENTS_DECIMAL)



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


        self.assertEqual(nFail, 0, "\n"+msg)

class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
