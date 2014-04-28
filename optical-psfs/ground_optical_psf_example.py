import numpy as np
import galsim
import ground_optical_psf

dz = -0.03
dx = -0.07
dy = 0.1
tx = 5.
ty = -3.


lam = 800            # nm 
diameter = 4.0       # meters
lam_over_diam = lam*1e-9/diameter*206265 # arcsec

# read ZEMAX files. Zernike coefficients were generated on 41x41 grids.
optical_psf_model = ground_optical_psf.OpticalPSFModel(lam = lam,
                                                       dz = dz, dx = dx, dy = dy,
                                                       tx = tx, ty = ty)
grid = [optical_psf_model.xmin, 0.5*optical_psf_model.xmin, 0.,
        0.5*optical_psf_model.xmax, optical_psf_model.xmax]

# make image
stamp_size = 80
n_grid = len(grid)
image = galsim.ImageF(stamp_size * n_grid, stamp_size * n_grid)
for iy, y in enumerate(grid):
    for ix, x in enumerate(grid):
        print "generating PSF at (%f, %f)" % (x, y)
        optics = optical_psf_model.get_psf(x, y)
        b = galsim.BoundsI(ix*stamp_size+1 , (ix+1)*stamp_size-1, 
                           iy*stamp_size+1 , (iy+1)*stamp_size-1)
        sub_image = image[b]
        optics.draw(image = sub_image, dx = lam_over_diam/2.)
filename = "ground_optical_psf_example.fits"
image.write(filename)
print "output is written in:", filename
