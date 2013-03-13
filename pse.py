import numpy as np
from numpy import pi

class PowerSpectrumEstimator(object):
    """
    Class for estimating the shear power spectrum from gridded shears.

    This class stores all the data used in power spectrum estimation that is fixed with the geometry
    of the problem - the binning and spin weighting factors.

    The only public method is `estimate()`, which is called with 2D g1 and g2 arrays on a square
    grid.  It assumes the flat sky approximation (where ell and k are interchangeable), and rebins
    the observed ell modes into a user-defined number of logarithimic bins in ell.  Given that the
    grid parameters are precomputed and stored when the PowerSpectrumEstimator is initialized,
    computation of the PS for multiple sets of shears corresponding to the same grid setup can
    proceed more rapidly than if everything had to be recomputed each time.

    Some important notes:

    1) Power spectrum estimation requires a weight function which decides how the averaging
    is done across ell within each bin.  By default, that weighting is flat in ell, but
    this is easy to change (in the _bin_power function).  A keyword allows for weighting by the
    power itself, but requires the GalSim software package.

    2) This is the power spectrum of the gridded *data*, not the underlying field - we do not
    account for the effects of the finite grid (basically, ignoring all the reasons why power
    spectrum estimation of the hard).  Users must account for the contribution of noise in g1, g2
    and the masking.

    3) The binning is currently fixed as uniform in log(ell).

    4) The code for this class uses the notation of the GREAT10 handbook (Kitching et al. 2011,
    http://dx.doi.org/10.1214/11-AOAS484), equations 17-21.
    """
    def __init__(self, N=100, sky_size_deg=10., nbin=15):
        """Create a PSE object given some grid parameters.

        The grid parameters are:
        N: the number of pixels along each side;
        sky_size_deg: the total grid width (in one dimension) in degrees; and
        nbin: the number of evenly-spaced logarithmic ell bins to use.

        Each parameter has a reasonable default value.
        """

        # Set up the scales of the sky and pixels
        self.N = N
        self.sky_size = np.radians(sky_size_deg)
        self.dx = self.sky_size / N

        # Define the possible ell range, the bin edges and effective ell values.
        # This is necessary for binning the power spectrum in ell.
        lmin = 2*pi / self.sky_size
        lmax = np.sqrt(2.)*pi / self.dx # in 2 dimensions
        self.bin_edges = np.logspace(np.log10(lmin), np.log10(lmax), nbin+1)
        # By default, report an area-averaged value of ell, which should be fine if there is
        # no weighting (in which case it's recomputed) and if there are many ell modes in
        # each bin.  The latter assumption is most likely to break down at low ell.
        self.ell = (2./3.)*(self.bin_edges[1:]**3-self.bin_edges[:-1]**3) \
                                   / (self.bin_edges[1:]**2-self.bin_edges[:-1]**2)

        # Precompute and store two useful factors, both in the form of 2D grids in Fourier space.
        # These are the lengths of the wavevector |ell| for each point in the space, and the complex
        # valued spin-weighting that takes g_r, g_i -> E,B
        self.l_abs, self.eb_rot = self._generate_eb_rotation()

    def _generate_eb_rotation(self):
        # Set up the Fourier space grid lx, ly.
        ell = 2*pi*np.fft.fftfreq(self.N, self.dx)
        lx = np.vstack([ell for i in xrange(self.N)])
        ly = lx.T

        # Now compute the lengths and angles of the ell vectors.
        l_abs = (lx**2 + ly**2)**0.5
        l_ang = np.arctan2(ly, lx)

        self.lx = lx
        self.ly = ly
        self.l_ang = l_ang

        # Compute the spin-2 weights.
        rot = np.exp(-2j*l_ang)

        return l_abs, rot

    def _bin_power(self, C, ell_weight=None):
        # This little utility function bins a 2D C^{E/B, E/B}_{ell} based on |ell|.  The use of
        # histogram is a little hack, but is quite convenient since it means everything is done in C
        # so it is very fast. The first call to `histogram` just returns an array over the
        # logarithmic ell bins of
        # sum_{|ell| in bin} C_{ell_x,ell_y}
        # and the second call returns
        # sum_{|ell| in bin} 1.
        # Thus, the ratio is just the mean power in the bin.  If `ell_weight` is not None, then some
        # non-flat weighting scheme is used for averaging over the ell values within a bin.
        if ell_weight is not None:
            P,_ = np.histogram(self.l_abs, self.bin_edges, weights=C*ell_weight)
            count,_ = np.histogram(self.l_abs, self.bin_edges, weights=ell_weight)
        else:
            P,_ = np.histogram(self.l_abs, self.bin_edges, weights=C)
            count,_ = np.histogram(self.l_abs, self.bin_edges)
        return P/count

    def estimate(self, g1, g2, weight_EE=False, weight_BB=False, theory_func=None):
        """ Compute the EE, BB, and EB power spectra of two 2D arrays g1 and g2.

        In addition to the obvious arguments (the shear components as 2D NumPy arrays), this method
        can take three optional arguments:
        * weight_EE and weight_BB determine whether the E and B auto-power spectra are re-computed
          weighting by the power within the bin.
        * theory_func is some callable function that can be used to get an idealized value of power
          at each point on the grid, and then see what results it gives for our chosen ell binning.
        """
        # Check for square geometry is what we expect.
        assert g1.shape == g2.shape == (self.N, self.N)

        # Transform g1+j*g2 into Fourier space and rotate into E-B, then separate into E and B.
        EB = np.fft.ifft2(self.eb_rot * np.fft.fft2(g1 + 1j*g2))
        E = np.fft.fft2(EB.real)
        B = np.fft.fft2(EB.imag)

        # Use the internal function above to bin, and account for the normalization of the FFT.
        # Recall that power has units of angle^2, which is the reason why we need a self.dx^2 in the
        # equations below in addition to the standard 1/N^2 coming from the FFTs.
        C_EE = self._bin_power(E*np.conjugate(E))*(self.dx/self.N)**2
        C_BB = self._bin_power(B*np.conjugate(B))*(self.dx/self.N)**2
        C_EB = self._bin_power(E*np.conjugate(B))*(self.dx/self.N)**2

        if theory_func or weight_EE or weight_BB:
            import galsim

        if theory_func is not None:
            # theory_func needs to be a callable function; need to carefully sanity check this
            # later, but for now just assume.
            C_theory_ell = np.zeros_like(self.l_abs)
            C_theory_ell[self.l_abs>0] = theory_func(self.l_abs[self.l_abs>0])
            C_theory = self._bin_power(C_theory_ell)

        if weight_EE or weight_BB:
            # Need to interpolate C_EE to values of self.l_abs.  A bit of kludginess as we go off
            # the end of our final ell grid...
            new_ell = np.zeros(len(self.ell)+2)
            new_ell[1:len(self.ell)+1] = self.ell
            new_ell[len(self.ell)+1] = 10.*max(self.ell)
            C_theory = self._bin_power(C_theory_ell, ell_weight=C_theory_ell)

        if weight_EE:
            new_CEE = np.zeros_like(new_ell)
            new_CEE[1:len(self.ell)+1] = np.real(C_EE)
            new_CEE[len(self.ell)+1] = new_CEE[len(self.ell)]
            EE_table = galsim.LookupTable(new_ell, new_CEE)
            ell_weight = EE_table(self.l_abs)
            C_EE = self._bin_power(E*np.conjugate(E), ell_weight=ell_weight)*(self.dx/self.N)**2

        if weight_BB:
            new_CBB = np.zeros_like(new_ell)
            new_CBB[1:len(self.ell)+1] = np.real(C_BB)
            new_CBB[len(self.ell)+1] = new_CBB[len(self.ell)]
            BB_table = galsim.LookupTable(new_ell, new_CBB)
            ell_weight = BB_table(self.l_abs)
            C_BB = self._bin_power(B*np.conjugate(B), ell_weight=ell_weight)*(self.dx/self.N)**2

        # For convenience, return ell (copied in case the user messes with it) and the three power
        # spectra. If the user requested a binned theoretical spectrum, return that as well.
        if theory_func is None:
            return self.ell.copy(), np.real(C_EE), np.real(C_BB), np.real(C_EB)
        else:
            return self.ell.copy(), np.real(C_EE), np.real(C_BB), np.real(C_EB), np.real(C_theory)
