import os
import codecs
import numpy as np
import galsim

def loadZernikeCoefficients(filename_coeff):
    """
    Helper function to read Zernike coefficients from a ZEMAX output.
    """
    defocus = 0.
    a1 = 0.
    a2 = 0.
    c1 = 0.
    c2 = 0.
    spher = 0.
    t1 = 0.
    t2 = 0.
    with codecs.open(filename_coeff, "r", "utf-16") as f:
        for line in f:
            item = line.split()
            if len(item) == 0:
                continue
            elif item[0] == "Wavelength":
                wavelength = float(item[2]) # um
            # following Zernike coefficients are in wavelength
            elif item[0] == "Z":
                if item[1] == "4":
                    defocus = float(item[2])
                elif item[1] == "5":
                    a1 = float(item[2])
                elif item[1] == "6":
                    a2 = float(item[2])
                elif item[1] == "7":
                    c1 = float(item[2])
                elif item[1] == "8":
                    c2 = float(item[2])
                elif item[1] == "9":
                    t1 = float(item[2])
                elif item[1] == "10":
                    t2 = float(item[2])
                elif item[1] == "11":
                    spher = float(item[2])
    return wavelength, np.array([defocus, a1, a2, c1, c2, t1, t2, spher])

class OpticalPSFMisalignment:
    """Misalignment model based on code from Aaron Roodman (SLAC National Accelerator Laboratory).

    This class calculates how Zernike coefficients respond to defocus(dz), decenter(dx, dy),
    and tilt(tx, ty). This change is also a function of position in FOV(x, y).

    After investigating the response by using ZEMAX, we found that the response at each position
    can be approximately modeled by a linear function of (dz, dx, dy, tx, ty), and interpolated
    by a linear function of (x, y) across the field.
    Thus the response of a Zernike coefficient is written as

    delta_Z = A_ij * u_j + B_ij * u_j * x + C_ij * u_j * y      (1),

    where delta_Z = (delta_defocus, delta_astig1, delta_astig2, delta_coma1, delta_coma2, 
    delta_trefoil1, delta_trefoil2, delta_spher) and u = t(dz, dx, dy, tx, ty).

    We found the following responses of Zernike coefficients:
    - defocus depends only on dz and does not depend on (x,y)
    - astigmatism depends on (dx, dy, tx, ty), and also depends on (x,y)
    - coma depends on (dx, dy, tx, ty), but does not depend on (x,y)
    - trefoil and spherical do not depend on misalignments.
    In addition, we also found some symmetries, and Eq. (1) can be changed as
    1) separate defocus part as delta_defocus = a * dz
    2) define matrix M as
       M = [ 0  b  c  0]
           [ b  0  0 -c]
           [ 0  d  e  0]
           [ d  0  0 -e]
       , and vector v as v = M * u. And then
       delta_astig1 = v[0] * x + v[1] * y
       delta_astig2 = v[1] * x - v[0] * y
       delta_coma1 = v[2]
       delta_coma2 = v[3]

    These findings agree with Aaron's model.
    We implemented these responses based on the ZEMAX simulations we made.
    """

    def __init__(self, lam, dz, dx, dy, tx, ty):
        """
        Inputs
        - lam: wavelength [nm]
        - dz: defocus [mm]
        - dx: decenter along x axis [mm]
        - dy: decenter along y axis [mm]
        - tx: tilt about x [arcsec]
        - tx: tilt about y [arcsec]

        Notes of system convention
        (dx, dy, dz): right-handed system where z-axis is along the light coming into mirror.
        (tx, ty): + is the rotation by which a right screw goes into the direction of an axis.
        """
        self.lam_orig = 800.*1e-9
        self.lam = lam*1e-9 # m
        # transfer matrix from misalignment to Zernike coefficients.
        # These coefficients are based on wavelength of 800 nm, so that the results sould be
        # converted to wavelength one wants to get.
        self.M = np.matrix([[0.0e+00,  7.4e-02,  4.6e-03,  0.0e+00],
                            [7.4e-02,  0.0e+00,  0.0e+00, -4.6e-03], 
                            [0.0e+00,  2.4e-01,  1.8e-03,  0.0e+00],
                            [2.4e-01,  0.0e+00,  0.0e+00, -1.8e-03]])

        # defocus transformation coefficient
        self.a = -6.0e+00

        # storage arrays for defocus, a1, a2, c1, and c2.
        self.zernDelta = np.zeros(5)
        self.zernThetax = np.zeros(5)
        self.zernThetay = np.zeros(5)

        # defocus
        self.zernDelta[0] = self.a * dz

        # decenter and tilt
        misalignmentList = (dx, dy, tx, ty)
        misalignmentColVec = np.matrix(misalignmentList).transpose()
        zernikeVector = self.M*misalignmentColVec
        self.zernThetax[1] = zernikeVector[0][0]
        self.zernThetay[1] = zernikeVector[1][0]
        self.zernThetax[2] = zernikeVector[1][0]
        self.zernThetay[2] = -zernikeVector[0][0]
        self.zernDelta[3] = zernikeVector[2][0]
        self.zernDelta[4] = zernikeVector[3][0]


    def apply(self, x, y):
        """
        Inputs
        - x: x position on FOV [deg]
        - y: y position on FOV [deg]

        Outputs
        - delta_coefs: ndarray, residual of Zernike coefficients at (x, y) caused by misalignments
                       [wavelength]
        """
        dzernike = self.zernDelta + self.zernThetax * x + self.zernThetay * y
        ddefocus = dzernike[0]
        da1 = dzernike[1]
        da2 = dzernike[2]
        dc1 = dzernike[3]
        dc2 = dzernike[4]
        dt1 = 0.
        dt2 = 0.
        dspher = 0.
        delta_coefs = np.array([ddefocus, da1, da2, dc1, dc2,
                                dt1, dt2, dspher])*self.lam_orig/self.lam
        return delta_coefs

class OpticalPSFModel:
    """This class is used for obtaining an optical PSF at an arbitrary position (x,y).

    One should prepare ZEMAX output files in which Zernike coefficients are stored, and then make
    a list of positions where the ZEMAX output files are generated.  These positions should be on
    a grid, where the number of grids along the x-axis should be the same as that along the y-axis.

    First the class reads the ZEMAX files following the list of positions.
    The ZEMAX files should be created by using a ZEMAX analysis function "Zernike Standard
    Coefficients", stored in the same directory as the position list, and named as
    "%.4f_%.4f.txt % (x, y)".  An example of ZEMAX Macro to generate these files is given in
    "GREAT3_COEFFICIENTS.ZPL".

    Then this class interpolates based the Zernike coefficients across the FOV.
    One can add misalignments (defocus, decenter, tilt) for which the class uses
    OpticalPSFMisalignment. Based on these Zernike coefficients, this class makes an optical PSF at
    an arbitrary position.
    One can specify wavelength in nm, diameter in m, and obscuration in ratio between mirror
    diameter and obscuration diameter.
    One can get minimum/maximum x/y of the data used for the interpolation. They are stored as
    variables xmin, xmax, ymin, and ymax.

    A recommended setup based on a proto-type DES optics design, which we got by a courtesy of
    Aaron Roodman (SLAC National Accelerator Laboratory) and Steve Kent (Fermilab), is adopted as a
    default of this class. Since this is based on an actual telescope, the results could get
    unrealistic if the diameter is changed by much. According to Aaron, typical values of
    misalignments are
    - decenter(dx, dy): 0.15 mm
    - tilt(tx, ty): 10 arcsec
    . We found typical value of defocus is likely
    - defocus(dz): 0.06 mm
    , which corresponds to ~20% increase of rms size of optical PSF. 
    We also found that the default configuration of the DES optics is not optimized in terms of
    focus. Thus we added offset to dz (dz0). For the DES optics,
    - dz0 = 0.05 mm
    gives the minimum rms size.
    """
    def __init__(self, position_list_filename = 
                 "ground_optical_psf_zernike_coefficients_41x41/ZEMAXInput.dat",
                 lam = 800., diameter = 4.0, obscuration = 0.35,
                 nstruts = 0, strut_thick = 0.01, strut_angle = 0.*galsim.degrees, 
                 pad_factor = None,
                 dz = 0., dx = 0., dy = 0., tx = 0., ty = 0., dz0 = 0.05,
                 interpolant2d = None):
        """
        Inputs
        - position_list_filename: filename of position list used for reading ZEMAX output files.
        - lam: wavelength [nm]
        - diameter: diameter of telescope [m]
        - obscuration: central obscuration [ratio between mirror diameter and obscuration diameter]
        - nstruts: number of radial support struts to add to the central obscuration
        - strut_thick: thickness of support struts as a fraction of pupil diameter
        - strut_angle: angle made between the vertical and the strut starting closest to it,
                       defined to be positive in the counter-clockwise direction; must be a
                       galsim.Angle instance
        - pad_factor: optional padding specification if 1.5 is not good enough
        - dz: defocus [mm]
        - dx: decenter along x axis [mm]
        - dy: decenter along y axis [mm]
        - tx: tilt about x [arcsec]
        - tx: tilt about y [arcsec]
        - dz0: offset to defocus [mm]
        - interpolant2d: galsim._galsim.InterpolantXY
                         If None, galsim.InterpolantXY(galsim.Quintic())

        Notes of system convention
        (dx, dy, dz): right-handed system where z-axis is along the light coming into mirror.
        (tx, ty): + is the rotation by which a right screw goes into the direction of an axis.
        """
        self.lam = lam*1e-9 # meters
        self.lam_over_diam = self.lam/diameter*206265 # arcsec
        self.obscuration = obscuration
        self.nstruts = nstruts
        self.strut_thick = strut_thick
        self.strut_angle = strut_angle
        self.pad_factor = pad_factor

        # read position information from the list
        data = np.loadtxt(position_list_filename)
        x = data[:,0]
        y = data[:,1]
        n = int(np.sqrt(len(x)))
        d = y[1] - y[0]
        self.xmin = x.min()
        self.xmax = x.max()
        self.ymin = y.min()
        self.ymax = y.max()

        # read coefficients from ZEMAX file
        self.ncoefs = 8 # defocus, a1, a2, c1, c2, t1, t2, spher
        coefs = np.zeros((self.ncoefs, n, n))
        for i, (_x, _y) in enumerate(zip(x, y)):
            zernike_filename = os.path.join(os.path.dirname(position_list_filename),
                                            "%.4f_%.4f.txt" % (_x,_y))
            wavelength, coefs_tmp = loadZernikeCoefficients(zernike_filename)
            i_x = int(i/n)
            i_y = i%n
            # need to flip sign of x and y, since definition of axes in ZEMAX and GalSim seem
            #different.
            # We have to convert Zernike coefficients in units of wavelength we want. 
            coefs[:, n - i_y - 1, n - i_x - 1] = coefs_tmp*wavelength*1e-6/self.lam
        
        # get interpolated images
        self.interpolated_coefficients = list()
        for coef in coefs:
            im_coef = galsim.ImageViewD(coef)
            im_coef.setScale(d)
            if interpolant2d == None:
                interpolant2d = galsim.InterpolantXY(galsim.Quintic())
            self.interpolated_coefficients.append(galsim.InterpolatedImage(im_coef,
                                                                   x_interpolant = interpolant2d,
                                                                   normalization = "sb",
                                                                   calculate_stepk = False,
                                                                   calculate_maxk = False,
                                                                   ))
        # prepare for misalignment
        self.optical_psf_misalignment = OpticalPSFMisalignment(self.lam*1e9,
                                                               dz + dz0, dx, dy, tx, ty)

    def get_zernike_coefficients(self, x, y):
        """
        Inputs
        - x: x position on FOV [deg]
        - y: y position on FOV [deg]

        Outputs
        - delta_coefs: ndarray, Zernike coefficients at (x, y) including effects by misalignments
                       [wavelength]
        """
        if x < self.xmin or x > self.xmax or y < self.ymin or y > self.ymax:
            import warnings
            warnings.warn(
                    "Warning: position (%f,%f) not within the bounds "%(x,y) +
                    "of the gridded values.  Min, max x: (%f,%f) "%(self.xmin,self.xmax) +
                    "and min, max y: (%f, %f) "%(self.ymin,self.ymax))
        coefs = list()
        for interpolated_coefficient in self.interpolated_coefficients:
            coefs.append(interpolated_coefficient.xValue(galsim.PositionD(x, y)))
        coefs = np.array(coefs)
        delta_coefs = self.optical_psf_misalignment.apply(x, y)
        coefs += delta_coefs
        return coefs

    def get_psf(self, x, y, additional_coefs = None):
        """
        - x: x position on FOV [deg]
        - y: y position on FOV [deg]
        - additional_coefs: ndarray (Noll ordering: defocus, astig1, astig2, coma1, coma2,
                            trefoil1, trefoil2, spher), additional Zernike coefficients.
                            If None, no additional errors are added [wavelength]

        Outputs
        - optics: galsim.optics.OpticalPSF
        """
        if x < self.xmin or x > self.xmax or y < self.ymin or y > self.ymax:
            import warnings
            warnings.warn(
                    "Warning: position (%f,%f) not within the bounds "%(x,y) +
                    "of the gridded values.  Min, max x: (%f,%f) "%(self.xmin,self.xmax) +
                    "and min, max y: (%f, %f) "%(self.ymin,self.ymax))
        coefs = self.get_zernike_coefficients(x, y)
        if additional_coefs is not None:
            coefs += additional_coefs
        optics = galsim.OpticalPSF(lam_over_diam = self.lam_over_diam, 
                                   defocus = coefs[0],
                                   astig1 = coefs[1],
                                   astig2 = coefs[2],
                                   coma1 = coefs[3],
                                   coma2 = coefs[4],
                                   trefoil1 = coefs[5],
                                   trefoil2 = coefs[6],
                                   spher = coefs[7],
                                   obscuration = self.obscuration,
                                   nstruts = self.nstruts,
                                   strut_thick = self.strut_thick,
                                   strut_angle = self.strut_angle,
                                   pad_factor = self.pad_factor,
                                   suppress_warning = True)
        return optics
