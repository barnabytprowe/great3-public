# check the shape of a 2-component model for the galaxies using adaptive moments

import galsim
import pyfits
import os
import numpy as np

# define filenames etc.
in_dir = '.'
in_filename = 'real_galaxy_catalog_23.5_fits.fits'
out_dir = '.'
out_filename = 'real_galaxy_23.5_shapes.fits'
pix_scale = 0.03

# define PSF - for later use, could turn into a list of different ones
airy = galsim.Airy(lam_over_diam=206265.*800.e-9/4.)
pix = galsim.Pixel(pix_scale)
psf = galsim.Convolve(airy, pix)
im_psf = psf.draw(dx=pix_scale)

# read in catalog
infile = os.path.join(in_dir, in_filename)
dat = pyfits.getdata(infile)
n = len(dat)
print "Read in ",n," from ",infile

# loop over objects
e1 = np.zeros(n)
e2 = np.zeros(n)
bt = np.zeros(n)
do_meas = np.zeros(n)
use_bulgefit = np.ones(n)
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

    # First decide if we can / should use bulgefit, otherwise sersicfit
    if bstat<1 or bstat>4 or dvc_btt<0.1 or dvc_btt>0.9 or np.isnan(dvc_btt) or params[9]<=0 or params[1]<=0 or params[11]<0.051 or params[3]<0.051 or smad<bmad:
        use_bulgefit[i] = 0
        # Then check if sersicfit is viable; if not, this object is a total failure:
        if sstat<1 or sstat>4 or sparams[1]<=0 or sparams[0]<=0:
            use_bulgefit[i] = -1
            do_meas[i] = -1
            e1[i] = -10.
            e2[i] = -10.
            bt[i] = -10.
            continue

    if use_bulgefit[i]:
        bulge_beta = params[15]*galsim.radians
        bulge_hlr = 0.03*np.sqrt(bulge_q)*params[9]
        bulge_flux = 2.0*np.pi*3.607*(bulge_hlr**2)*params[8]
        bulge_q = params[11]

        disk_beta = params[7]*galsim.radians
        disk_hlr = 0.03*np.sqrt(disk_q)*params[1]
        disk_flux = 2.0*np.pi*1.901*(disk_hlr**2)*params[0]
        disk_q = params[3]
        
        bfrac = bulge_flux/(bulge_flux+disk_flux)

        if bfrac < 0 or bfrac > 1 or np.isnan(bfrac):
            e1[i] = -10.
            e2[i] = -10.
            do_meas[i] = -1
            use_bulgefit[i] = -1
            continue

        bt[i] = bfrac
        bulge = galsim.Sersic(4.0, half_light_radius = bulge_hlr, flux = bulge_flux)
        disk = galsim.Exponential(half_light_radius = disk_hlr, flux = disk_flux)
        if bulge_q > 0.:
            bulge.applyShear(q = bulge_q, beta = bulge_beta)
        if disk_q > 0.:
            disk.applyShear(q = disk_q, beta = disk_beta)
        gal = bulge+disk
        obj = galsim.Convolve(gal, psf, gsparams = galsim.GSParams(maximum_fft_size=15000))
        try:
            im = obj.draw(dx = pix_scale)
            res = galsim.hsm.EstimateShear(im, im_psf, guess_sig_gal=20)
            e1[i] = res.corrected_e1
            e2[i] = res.corrected_e2
            do_meas[i] = 1.
        except RuntimeError:
            e1[i] = -10.0
            e2[i] = -10.0
            do_meas[i] = -1.
    else:
        e = galsim.Shear(q = sparams[3], beta = sparams[7]*galsim.radians)
        e1[i] = e.e1
        e2[i] = e.e2

# save results to file
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=dat.field('IDENT')),
                                         pyfits.Column(name='bulge_tot',
                                                       format='D',
                                                       array=bt),
                                         pyfits.Column(name='e1',
                                                       format='D',
                                                       array=e1),
                                         pyfits.Column(name='e2',
                                                       format='D',
                                                       array=e2),
                                         pyfits.Column(name='do_meas',
                                                       format='J',
                                                       array=do_meas),
                                         pyfits.Column(name='use_bulgefit',
                                                       format='J',
                                                       array=use_bulgefit)]
                                        ))

fail_ind = np.where(do_meas < -0.5)[0]
print len(fail_ind),' failures'

# write outputs
outfile = os.path.join(out_dir, out_filename)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)

