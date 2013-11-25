# check various properties of a 2-component model for the training galaxies, as needed for GREAT3

import galsim
import pyfits
import os
import numpy as np
import math

def training_galaxy_props_real(psf,
                               in_dir = '/home/rmandelb.proj/data-shared/great3_fit_data',
                               in_filename = 'real_galaxy_catalog_23.5.fits',
                               out_dir = '.',
                               out_filename = 'real_galaxy_catalog_23.5_real_props.fits',
                               pix_scale = 0.03,
                               size_factor = 0.6,
                               ps_size = 48,
                               do_orig = False,
                               n_use = None): # change to None!
    """
    A routine for estimating the properties of the training sample galaxies in different imaging
    conditions.
    
    Given input PSF and pixel scale, we make a real galaxy image in the same way as for GREAT3, and check the
    following:
    * did the original image have reasonable S/N within an elliptical Gaussian filter?  (Note:
      result does not depend on final PSF / pixel scale.)
    * what is the minimum noise variance in a whitened image?

    @params psf           GalSim object representing the PSF (without pixel convolution).
    @params in_dir        Directory containing the catalog of fit parameters for the sample.
    @params in_filename   Name of catalog of fit parameters for the sample.
    @params out_dir       Directory in which to put output file.
    @params out_filename  Name of output catalog
    @params pix_scale     Pixel scale for images
    @params size_factor   Multiplicative factor by which to modify galaxy sizes to represent a deeper
                          sample.
    @params ps_size       Number of pixels per side for postage stamp into which to draw image
    @params do_orig       Measure original PS?  Or just simulated.
    @params n_use         Number of galaxies to use; =None for using whole catalog
"""

    # define PSF
    pix = galsim.Pixel(pix_scale)
    epsf = galsim.Convolve(psf, pix)
    im_epsf = epsf.draw(dx=pix_scale)

    # set up RealGalaxyCatalog object
    rgc = galsim.RealGalaxyCatalog(in_filename, dir=in_dir)
    n = rgc.nobjects
    print "Read in ",n," from ",in_filename

    if n_use is not None:
        print "Using ",n_use
        n = n_use

    # loop over objects
    sn_ellip_gauss = np.zeros(n)
    min_var_white = np.zeros(n)
    for i in range(n):
        if i % 1000 == 0:
            print "...",i

        # Make the RealGalaxy object
        rg = galsim.RealGalaxy(rgc, index=i)

        # First do the test of S/N for original image:
        if do_orig:
            orig_var = rgc.variance[i] / 0.316 # fudge factor for correlated noise
            try:
                res = rg.original_image.draw(dx=0.03).FindAdaptiveMom()
                aperture_noise = np.sqrt(orig_var*2.*np.pi*(res.moments_sigma**2))
                sn_ellip_gauss[i] = res.moments_amp / aperture_noise
            except:
                sn_ellip_gauss[i] = -10.
        
        # Now make the simulated object, and check the minimum noise variance after whitening.
        # First we rescale the size of the object - could significantly change noise properties.
        rg.applyDilation(size_factor)
        # No need to apply shear/magnification since a small shear or magnification induces minimal
        # changes in noise properties.  So just convolve with target PSF and draw at target pixel
        # scale.
        obj = galsim.Convolve(rg, epsf)
        im = galsim.ImageF(ps_size, ps_size)
        try:
            im = obj.draw(dx = pix_scale)
            min_var_white[i] = obj.noise.applyWhiteningTo(im)
        except:
            min_var_white[i] = -10.

    # save results to file
    tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='sn_ellip_gauss',
                                                           format='D',
                                                           array=sn_ellip_gauss),
                                             pyfits.Column(name='min_var_white',
                                                           format='D',
                                                           array=min_var_white)]
                                            ))

    # a sanity check:
    print len(sn_ellip_gauss[sn_ellip_gauss<20.])," have S/N<20"
    
    # write outputs
    outfile = os.path.join(out_dir, out_filename)
    print "Writing to file ",outfile
    tbhdu.writeto(outfile, clobber=True)
    
