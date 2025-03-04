.. py:currentmodule:: lsst.meas.extensions.shapeHSM

HSMShape Plugins
================

The purpose of the `HsmShape` algorithms is to perform PSF correction to measure galaxy shapes for the purpose of weak lensing.
The plugins can do several methods of PSF correction, described and tested on simulated data in `Hirata & Seljak (2003) <http://adsabs.harvard.edu/abs/2003MNRAS.343..459H>`_.

The PSF correction methods are

1. :py:class:`KSB <lsst.meas.extensions.shapeHSM.HsmShapeKsbPlugin>` --
   An implementation of the method from `Kaiser, Squires, & Broadhurst (1995) <http://adsabs.harvard.edu/abs/1995ApJ...449..460K>`_.

2. :py:class:`BJ <lsst.meas.extensions.shapeHSM.HsmShapeBjPlugin>` --
   Method from Appendix C of `Bernstein & Jarvis (2002) <http://adsabs.harvard.edu/abs/2002AJ....123..583B>`_ to apply a
   rounding kernel, and then correct the galaxy image moments using an
   equation that is valid in the limit that the galaxy and PSF are
   approximately Gaussian.

3. :py:class:`Linear <lsst.meas.extensions.shapeHSM.HsmShapeLinearPlugin>` --
   A variation of the BJ method described in Appendix B of `Hirata & Seljak (2003) <http://adsabs.harvard.edu/abs/2003MNRAS.343..459H>`_.
   This method accounts for the first order departure of both the galaxy and PSF from Gaussianity.

4. :py:class:`Re-Gaussianization <lsst.meas.extensions.shapeHSM.HsmShapeRegaussPlugin>` --
   Method described in section 2.4 of `Hirata & Seljak (2003) <http://adsabs.harvard.edu/abs/2003MNRAS.343..459H>`_ that includes a correction for the full departure of the PSF from Gaussianity.
   Note that this method has been tested more extensively in, e.g., `Mandelbaum et al. (2005) <http://adsabs.harvard.edu/abs/2005MNRAS.361.1287M>`_ and `Massey et al. (2007) <http://adsabs.harvard.edu/abs/2007MNRAS.376...13M>`_ for ground-based data only.

The shapes are defined using :math:`e = (1-q^2)/(1+q^2)` for an axis ratio :math:`0 \lt q \le 1`.
The output values from these algorithms may not satisfy conditions such as :math:`|e| \lt 1`, so they are reported as :math:`e_1` and :math:`e_2` values instead of quadrupole moments (i.e., :math:`xx`, :math:`yy` and :math:`xy`).
