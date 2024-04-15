.. py:currentmodule:: lsst.meas.extensions.shapeHSM

PSF-correction algorithms
=========================

The purpose of the `HsmShape` algorithms is to perform PSF correction to
measure galaxy shapes for the purpose of weak lensing.  The plugins can do
several methods of PSF correction, described and tested on simulated
data in Hirata & Seljak (2003),
http://adsabs.harvard.edu/abs/2003MNRAS.343..459H.

The PSF correction methods are

1. KSB --
   An implementation of the method from Kaiser, Squires, & Broadhurst (1995)
   described in http://adsabs.harvard.edu/abs/1995ApJ...449..460K.

2. BJ --
   Method from Appendix C of Bernstein & Jarvis (2002) to apply a
   rounding kernel, and then correct the galaxy image moments using an
   equation that is valid in the limit that the galaxy and PSF are
   approximately Gaussian.

3. Linear --
   A variation of the BJ method described in Appendix B of
   Hirata & Seljak (2003). This method accounts for the first
   order departure of both the galaxy and PSF from Gaussianity.

4. Re-Gaussianization --
   Method described in section 2.4 of Hirata & Seljak (2003) that includes a
   correction for the full departure of the PSF from Gaussianity.  Note that
   this method has been tested more extensively in, e.g.,
   http://adsabs.harvard.edu/abs/2005MNRAS.361.1287M
   and http://adsabs.harvard.edu/abs/2007MNRAS.376...13M
   for ground-based data only.

The shapes are defined using :math:`e = (1-q^2)/(1+q^2)` for an axis ratio
:math:`0 \lt q \le 1`.
The output values from these algorithms may not satisfy conditions such as
:math:`|e| \lt 1`, so they are reported back as :math:`e_1` and :math:`e_2`
values instead of quadrupole moments
(i.e., :math:`xx`, :math:`yy` and :math:`xy`).
