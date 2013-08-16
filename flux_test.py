import pyfits
import numpy as np
import math
import galsim
import os
import matplotlib.pyplot as plt
import scipy.stats as stats

# define dirs, filenames, etc.
dir = '/Users/rmandelb/great3/data-23.5'
fit_catalog = 'real_galaxy_catalog_23.5_fits.fits'
sel_catalog = 'real_galaxy_selection_info.fits'
out_catalog = 'real_galaxy_deltamag_info.fits'

# read in the fit and selection catalogs
fit_data = pyfits.getdata(os.path.join(dir, fit_catalog))
sel_data = pyfits.getdata(os.path.join(dir, sel_catalog))
use_bulgefit = sel_data['use_bulgefit'][:,0]
n = len(fit_data)
if len(sel_data) != n:
    raise RuntimeError("Catalog lengths do not match!")
print "Read in ",n," from fit and selection catalogs"

# for each object, calculate the flux (according to whether we use
# bulgefit or sersicfit)
claire_flux = np.zeros(n)
print "Getting fluxes..."
for i in range(n):
    if i % 1000 == 0:
        print "    Galaxy ",i
    if use_bulgefit[i] == 1.:
        params = fit_data[i].field('bulgefit')

        bulge_q = params[11]
        bulge_hlr = 0.03*np.sqrt(bulge_q)*params[9] # arcsec 
        # It is needed if we're drawing with 'flux' normalization convention.
        bulge_flux = 2.0*np.pi*3.607*(bulge_hlr**2)*params[8]/(0.03**2)

        disk_q = params[3]
        disk_hlr = 0.03*np.sqrt(disk_q)*params[1] # arcsec
        disk_flux = 2.0*np.pi*1.901*(disk_hlr**2)*params[0]/(0.03**2)

        claire_flux[i] = bulge_flux + disk_flux
    else:
        # Make a single Sersic model instead
        params = fit_data[i].field('sersicfit')
        gal_n = params[2]
        # Fudge this if it is at the edge.  Now that GalSim #325 and #449 allow Sersic n in
        # the range 0.3<=n<=6, the only problem is that Claire occasionally goes as low as
        # n=0.2.
        gal_q = params[3]
        gal_hlr = 0.03*np.sqrt(gal_q)*params[1]
        if gal_n < 0.3: gal_n = 0.3
        tmp_ser = galsim.Sersic(gal_n, half_light_radius=1.)
        gal_bn = (1./tmp_ser.getScaleRadius())**(1./gal_n)
        prefactor = gal_n * math.gamma(2.*gal_n) * math.exp(gal_bn) / (gal_bn**(2.*gal_n))
        claire_flux[i] = 2.*np.pi*prefactor*(gal_hlr**2)*params[0]/0.03**2
print len(claire_flux[claire_flux<=0])," have zero or negative fluxes; fudging!"
claire_flux[claire_flux<=0] = 1e-5
claire_mag = -2.5*np.log10(claire_flux)
catalog_mag = fit_data[0:n].field('mag_auto')

# do a linear regression - thsi should basically work if there aren't many outliers
(a, b, r, tt, stderr) = stats.linregress(catalog_mag, claire_mag)
print "Linear regression params: ",a,b
best_fit = a*catalog_mag + b
# now do an outlier rejection with delta mag ~ 0.5
fit_diff = claire_mag - best_fit
to_use = np.abs(fit_diff)<0.5
print "After outlier rejection, using ",np.sum(to_use.astype(int))," of ",len(to_use)
(a, b, r, tt, stderr) = stats.linregress(catalog_mag[to_use], claire_mag[to_use])
print "Slope and intercept are ",a,b
best_fit = a*catalog_mag + b

# compare -2.5 log10(flux) with quoted magnitude
# scatter plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(catalog_mag, claire_mag, label='Data')
ax.plot(catalog_mag, best_fit, label='Linear fit, a=%f, b=%f'%(a,b))
ax.set_xlabel('Catalog mag_auto')
ax.set_ylabel('Claire magnitude')
plt.legend(loc='upper left')
plt.savefig('claire_mags_scatter.png')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(catalog_mag, claire_mag, label='Data')
ax.plot(catalog_mag, best_fit, label='Linear fit, a=%f, b=%f'%(a,b))
ax.set_xlabel('Catalog mag_auto')
ax.set_ylabel('Claire magnitude')
ax.set_xlim(20,24)
ax.set_ylim(-8,0)
plt.legend(loc='upper left')
plt.savefig('claire_mags_scatter_zoom.png')

# histogram
fit_diff = claire_mag - best_fit
fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(fit_diff, 100, range=(-2,2))
ax.set_xlabel('Magnitude difference, Claire - linear fit')
plt.savefig('claire_mag_histdiff.png')
print "Number with difference >1 mag (either direction):",len(fit_diff[np.abs(fit_diff)>=1])
print "Number with difference >0.5 mag (either direction):",len(fit_diff[np.abs(fit_diff)>=0.5])
print "Number that are >1 mag fainter:",len(fit_diff[fit_diff>=1])
print "Difference for UFO that Barney identified,",fit_diff[fit_data.field('ident')==102152]

tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=fit_data.field('ident')),
                                         pyfits.Column(name='delta_mag',
                                                       format='D',
                                                       array=fit_diff)]
                                        ))
outfile = os.path.join(dir,out_catalog)
print "Writing to file ",outfile
tbhdu.writeto(outfile, clobber=True)
