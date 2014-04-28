import space_optical_psf
import galsim
import numpy as np

rms = 1./13. # wavelentgh
tel_diam = 2.4 # meters
f_number = 8.1
dx = 1.e-6/tel_diam/f_number*206265 # arcsec

# read wfirst fits file
optical_psf_model = space_optical_psf.OpticalPSFModel(rms = rms)
xgrid = [optical_psf_model.xmin, 0.5*optical_psf_model.xmin,
         0., 0.5*optical_psf_model.xmax, optical_psf_model.xmax]
ygrid = [optical_psf_model.ymin, 0.5*optical_psf_model.ymin,
         0., 0.5*optical_psf_model.ymax, optical_psf_model.ymax]


# make image
stamp_size = 128
n_grid = len(xgrid)
image = galsim.ImageF(stamp_size * n_grid, stamp_size * n_grid)
additional_errors = np.array([1e-7, 0., 0., 0., 0., 0., 0., 0.])
for iy, y in enumerate(ygrid):
    for ix, x in enumerate(xgrid):
        print "generating PSF at (%f, %f)" % (x, y)
        optics = optical_psf_model.get_psf(x, y)

        b = galsim.BoundsI(ix*stamp_size+1 , (ix+1)*stamp_size-1, 
                           iy*stamp_size+1 , (iy+1)*stamp_size-1)
        sub_image = image[b]
        optics.draw(image = sub_image, dx = dx)
filename = "space_optical_psf_example.fits"
image.write(filename)
print "output is written in:", filename
