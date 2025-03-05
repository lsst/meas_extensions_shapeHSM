.. py:currentmodule:: lsst.meas.extensions.shapeHSM

HSMMoments Plugins
==================

The purpose of the `HsmMoments` algorithms is to measure shapes of both the sources (or objects) and PSFs.
The plugins wrap the implementation of the adaptive moments algorithm described in `Bernstein & Jarvis (2002) <http://adsabs.harvard.edu/abs/2002AJ....123..583B>`_ from `GalSim <https://github.com/GalSim-developers/GalSim>`_.
This implementation is known to be more robust than :py:class:`SdssShape <lsst.meas.base.SdssShapeAlgorithm>` algorithm.

The various plugins are specialized for sources (`HsmSourceMomentsPlugin` and its subclasses) and PSFs (`HsmPsfMomentsPlugin` and its subclasses).
These specializations allows an object to be measured under different configurations simultaneously.
For instance, `HsmSourceMomentsRoundPlugin` is specialized to use circular adaptive Gaussians to measure the moments of sources.
The `HsmPsfMomentsDebiasedPlugin` adds noise to the PSF image to degrade it to have the same signal-to-noise ratio (SNR) as the source image.
This makes the ellipticity calculated from this plugin have the same bias as the source ellipticity
The PSF moments from this plugin should be used when calculating ellipticity residuals so the bias is largely cancelled.
