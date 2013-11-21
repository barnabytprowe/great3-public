# check various properties of a 2-component model for the training galaxies, as needed for GREAT3

import galsim
import pyfits
import os
import numpy as np
import math

def training_galaxy_props(psf,
                          in_dir = '.',
                          in_filename = 'real_galaxy_catalog_23.5_fits.fits',
                          out_dir = '.',
                          out_filename = 'real_galaxy_catalog_23.5_props.fits',
                          pix_scale = 0.03,
                          size_factor = 0.6,
                          ps_size = 48,
                          n_use = None):
    """
    A routine for estimating the properties of the training sample galaxies in different imaging
    conditions.
    
    Given input PSF and pixel scale, we make a galaxy image in the same way as for GREAT3, and check the
    following:
    * what noise variance would yield S/N=20?
    * what is the resolution given that PSF?
    * what fraction of the flux is contained in the fiducial postage stamp size?

    @params psf           GalSim object representing the PSF (without pixel convolution).
    @params in_dir        Directory containing the catalog of fit parameters for the sample.
    @params in_filename   Name of catalog of fit parameters for the sample.
    @params out_dir       Directory in which to put output file.
    @params out_filename  Name of output catalog
    @params pix_scale     Pixel scale for images
    @params size_factor   Multiplicative factor by which to modify galaxy sizes to represent a deeper
                          sample.
    @params ps_size       Number of pixels per side for postage stamp into which to draw image
    @params n_use         Number of galaxies to use; =None for using whole catalog
"""

    # define PSF
    pix = galsim.Pixel(pix_scale)
    epsf = galsim.Convolve(psf, pix)
    im_epsf = epsf.draw(dx=pix_scale)

    # read in catalog
    infile = os.path.join(in_dir, in_filename)
    dat = pyfits.getdata(infile)
    n = len(dat)
    print "Read in ",n," from ",infile

    if n_use is not None:
        dat = dat[0:n_use]
        print "Using ",n_use
        n = n_use

    # Create output arrays for the following quantities: bulge-to-total flux ratio B/T (bt),
    # flux_frac (fraction of the flux included in a GREAT3 postage stamp), resolution factor based
    # on adaptive moments, noise_var_snr_20 (noise variance to make the object have a S/N of 20) and
    # use_bulgefit (should we use 2-component fit, or 1-component Sersic fit).
    bt = np.zeros(n)
    do_meas = np.zeros(n)
    flux_frac = np.zeros(n)
    resolution = np.zeros(n)
    noise_var_snr_20 = np.zeros(n)
    use_bulgefit = np.ones(n)
    # loop over objects
    for i in range(n):
        if i % 1000 == 0:
            print "...",i
        params = dat[i].field('bulgefit')
        sparams = dat[i].field('sersicfit')
        bstat = dat[i].field('fit_status')[0]
        sstat = dat[i].field('fit_status')[4]
        dvc_btt = dat[i].field('fit_dvc_btt')
        bmad = dat[i].field('fit_mad_b')
        smad = dat[i].field('fit_mad_s')

        if bstat<1 or bstat>4 or dvc_btt<0.1 or dvc_btt>0.9 or np.isnan(dvc_btt) or params[9]<=0 or params[1]<=0 or params[11]<0.051 or params[3]<0.051 or smad<bmad:
            use_bulgefit[i] = 0
            # Then check if sersicfit is viable; if not, this object is a total failure:
            if sstat<1 or sstat>4 or sparams[1]<=0 or sparams[0]<=0:
                use_bulgefit[i] = -1
                do_meas[i] = -1
                resolution[i] = -10.
                noise_var_snr_20[i] = -10.
                flux_frac[i] = -10.
                bt[i] = -10.
                continue

        if use_bulgefit[i]:
            bulge_q = params[11]
            bulge_beta = params[15]*galsim.radians
            bulge_hlr = 0.03*size_factor*np.sqrt(bulge_q)*params[9]
            bulge_flux = 2.0*np.pi*3.607*(bulge_hlr**2)*params[8]

            disk_q = params[3]
            disk_beta = params[7]*galsim.radians
            disk_hlr = 0.03*size_factor*np.sqrt(disk_q)*params[1]
            disk_flux = 2.0*np.pi*1.901*(disk_hlr**2)*params[0]
            
            bfrac = bulge_flux/(bulge_flux+disk_flux)

            if bfrac < 0 or bfrac > 1 or np.isnan(bfrac):
                use_bulgefit[i] = -1
                do_meas[i] = -1
                resolution[i] = -10.
                noise_var_snr_20[i] = -10.
                flux_frac[i] = -10.
                bt[i] = -10.
                continue

            bt[i] = bfrac
            bulge = galsim.Sersic(4.0, half_light_radius = bulge_hlr, flux = bulge_flux)
            disk = galsim.Exponential(half_light_radius = disk_hlr, flux = disk_flux)
            if bulge_q > 0.:
                bulge.applyShear(q = bulge_q, beta = bulge_beta)
            if disk_q > 0.:
                disk.applyShear(q = disk_q, beta = disk_beta)
            gal = bulge+disk
            gal_flux = bulge_flux + disk_flux

        else:
            gal_n = sparams[2]
            if gal_n < 0.3: gal_n = 0.3
            tmp_ser = galsim.Sersic(gal_n, half_light_radius=1.)
            gal_bn = (1./tmp_ser.getScaleRadius())**(1./gal_n)
            prefactor = gal_n * math.gamma(2.*gal_n) * math.exp(gal_bn) / (gal_bn**(2.*gal_n))
            gal_q = sparams[3]
            gal_beta = sparams[7]*galsim.radians
            gal_hlr = 0.03*size_factor*np.sqrt(gal_q)*sparams[1]
            gal_flux = 2.*np.pi*prefactor*(gal_hlr**2)*params[0]

            gal = galsim.Sersic(gal_n, half_light_radius = gal_hlr, flux = gal_flux)
            if gal_q > 0.:
                gal.applyShear(q = gal_q, beta = gal_beta)

        obj = galsim.Convolve(gal, epsf, gsparams = galsim.GSParams(maximum_fft_size=15000))
        im = galsim.ImageF(ps_size, ps_size)
        try:
            im = obj.draw(dx = pix_scale)
        except:
            do_meas[i] = -0.5 # fit parameters make object impossible to draw
            resolution[i] = -10.
            noise_var_snr_20[i] = -10.
            flux_frac[i] = -10.
        flux_frac[i] = im.array.sum()/gal_flux
        noise_var_snr_20[i] = np.sum(im.array**2) / 20.**2

        # try getting resolution
        try:
            result = galsim.hsm.EstimateShear(im, im_epsf, guess_sig_gal=20)
            resolution[i] = result.resolution_factor
            do_meas[i] = 1. # made model and was able to measure shape
        except RuntimeError:
            resolution[i] = -10.
            do_meas[i] = 0. # made model and did other calculations, but could not measure shape

    # save results to file
    tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                           format='J',
                                                           array=dat.field('IDENT')),
                                             pyfits.Column(name='bulge_tot',
                                                           format='D',
                                                           array=bt),
                                             pyfits.Column(name='flux_frac',
                                                           format='D',
                                                           array=flux_frac),
                                             pyfits.Column(name='resolution',
                                                           format='D',
                                                           array=resolution),
                                             pyfits.Column(name='noise_var_snr_20',
                                                           format='D',
                                                           array=noise_var_snr_20),
                                             pyfits.Column(name='do_meas',
                                                           format='D',
                                                           array=do_meas),
                                             pyfits.Column(name='use_bulgefit',
                                                           format='D',
                                                           array=use_bulgefit)]
                                            ))
    
    fail_ind = np.where(do_meas < -0.5)[0]
    print len(fail_ind),' failures'
    bulgefit_ind = np.where(use_bulgefit == 1)[0]
    print len(bulgefit_ind),' use 2-component fits'

    # write outputs
    outfile = os.path.join(out_dir, out_filename)
    print "Writing to file ",outfile
    tbhdu.writeto(outfile, clobber=True)
    
