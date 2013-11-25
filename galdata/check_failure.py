import pyfits
import numpy as np
import matplotlib.pyplot as plt

# define filenames etc.
fit_cat = '/Users/rmandelb/great3/data-23.5/real_galaxy_catalog_23.5_fits.fits'
shape_cat = '/Users/rmandelb/great3/data-23.5/real_galaxy_23.5_shapes.fits'
outfile = 'diagnostics/shape_failures.png'

# read in fit information, etc.
fit_dat = pyfits.getdata(fit_cat)
n_fit = len(fit_dat)
print 'Read in ',n_fit,' from ',fit_cat

# read in measurement information
shape_dat = pyfits.getdata(shape_cat)
n_shape = len(shape_dat)
print 'Read in ',n_shape,' from ',shape_cat
if n_fit != n_shape:
    raise RuntimeError("Inconsistent catalog sizes")
do_meas = shape_dat.field('do_meas')
do_meas_fail = do_meas[do_meas<-0.5]
print len(do_meas_fail),' failures'

# histograms: magnitude
mag = fit_dat.field('mag_auto')
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.hist(mag, normed=True, facecolor='blue', alpha=0.5, bins=0.5*np.arange(11)+18, label='All')
ax1.hist(mag[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.25*np.arange(26)+18, label='Failures')
ax1.set_xlabel('F814W magnitude')

# histograms: size
size = 0.03*fit_dat.field('flux_radius')
ax2 = fig.add_subplot(222)
ax2.hist(size, normed=True, facecolor='blue', alpha=0.5, bins=0.05*np.arange(40), label='All')
ax2.hist(size[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.05*np.arange(40), label='Failures')
ax2.set_xlabel('flux_radius')

# histograms: B/T, or bulge-to-total flux ratio
bt = shape_dat.field('bulge_tot')
ax3 = fig.add_subplot(223)
ax3.hist(bt, normed=True, facecolor='blue', alpha=0.5, bins=0.05*np.arange(21), label='All')
ax3.hist(bt[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.05*np.arange(21), label='Failures')
ax3.set_xlabel('B/T')

# histograms: n
n = fit_dat.field('sersicfit')[:,2]
ax3 = fig.add_subplot(224)
ax3.hist(n, normed=True, facecolor='blue', alpha=0.5, bins=0.5*np.arange(17), label='All')
ax3.hist(n[do_meas<-0.5], normed=True, facecolor='red', alpha=0.5, bins=0.5*np.arange(17), label='Failures')
ax3.set_xlabel('n')

plt.legend(loc='upper right')
fig.tight_layout()
plt.savefig(outfile)
