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

__all__ = ["configure_hsm"]

import lsst.meas.extensions.shapeHSM  # noqa: F401 needed for plugins below


def configure_hsm(config):
    """Enable HSM shape measurements for a single frame measurement task.

    Parameters
    ----------
    config : `lsst.meas.base.SingleFrameMeasurementConfig`
        Measurement config to set up to use HSM shape; modified in place.
    """
    config.plugins.names |= [
        "ext_shapeHSM_HsmShapeRegauss",
        "ext_shapeHSM_HsmSourceMoments",
        "ext_shapeHSM_HsmPsfMoments",
        "ext_shapeHSM_HsmSourceMomentsRound",
        "ext_shapeHSM_HigherOrderMomentsSource",
        "ext_shapeHSM_HigherOrderMomentsPSF",
        "ext_shapeHSM_HsmPsfMomentsDebiased",
    ]
    config.slots.shape = "ext_shapeHSM_HsmSourceMoments"
    config.slots.psfShape = "ext_shapeHSM_HsmPsfMoments"
    config.plugins["ext_shapeHSM_HsmShapeRegauss"].deblendNChild = "deblend_nChild"
