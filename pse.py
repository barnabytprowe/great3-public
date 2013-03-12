import numpy as np
from numpy import pi
import galsim

class PowerSpectrumEstimator(object):
	"""
	This class stores all the data used in power spectrum estimation that 
	is fixed with the geometry of the problem - the binning and spin weighting factors.

	The only public method is estimate, which you call with 2D g1 and g2 arrays on a square
	grid.  It assumes the flat sky approximation.

	Some important notes:
	1) Power spectrum estimation requires a weight function which decides how the averaging
	is done across ell within each bin.  That weighting is implicitly flat in ell here, but
	this is easy to change (in the _bin_power function)

	2) This is the power spectrum of the *data*, not the underlying field - we do 
	not account for the effects that the finite mask size has (basically, all the reasons
	that realistic PSE is hard are ignored).  In particular we have to account for contribution
	of noise on g1, g2 and the masking.

	3) The binning is currently fixed as uniform in log(ell).  It is easy to change this to
	either something user-defined or to add an option to get linear spacing instead.

	4) The code for this class uses the notation of the Great10 handbook, equations 17-21.
	"""
	def __init__(self, N, sky_size_deg, nbin):
		"""Create a PSE object with the number of pixels along each side, N, 
		the total sky width in degrees, and the number of evenly-log-spaced ell bins to use."""

		#Set up the scales of the sky and pixels
		self.N=N
		self.sky_size = np.radians(sky_size_deg)
		self.dx = self.sky_size / N
		
		#Set the possible ell range and the bin edges and centers
		#This is for binning the power spectrum in ell.
		lmin = 2*pi / self.sky_size
		lmax = np.sqrt(2.)*pi / self.dx # in 2 dimensions
		self.bin_edges = np.logspace(np.log10(lmin), np.log10(lmax), nbin+1)
                # By default, report an area-averaged value of ell, which should be fine if there is
		# no weighting (in which case it's recomputed) and if there are many ell modes in
		# each bin.  The latter assumption is most likely to break down at low ell.
		self.ell = (2./3.)*(self.bin_edges[1:]**3-self.bin_edges[:-1]**3) \
                                   / (self.bin_edges[1:]**2-self.bin_edges[:-1]**2)

		#Compute two useful factors, both in the form of 2D grids in Fourier space.
		#These are the lengths of the wavevector |ell| for each point in the space,
		#and the complex valued spin-weighting that takes g_r, g_i -> E,B
		self.l_abs, self.eb_rot = self._generate_eb_rotation()

	def _generate_eb_rotation(self):
		#Set up the Fourier space grid lx, ly
		#Restrict ourselves temporarily to nx=ny so that the grid is 
		#simpler
		ell=2*pi*np.fft.fftfreq(self.N, self.dx)
		lx=np.vstack([ell for i in xrange(self.N)])
		ly=lx.T

		# Now compute the lengths and angles of the ell vectors
		l_abs = (lx**2 + ly**2)**0.5
		l_ang = np.arctan2(ly, lx)

		self.lx=lx
		self.ly=ly
		self.l_ang=l_ang

		#and the spin-2 weights.  And return both.
		rot = np.exp(-2j*l_ang)
		return l_abs, rot

	def _bin_power(self, C, ell_weight=None):
		#This little utility function bins a 2D C^{E/B, E/B}_{ell}
		#based on |ell|.  The use of histogram is a little hack,
		#but is quite convenient since it means everything is done in C
		#so it is very fast. The first one just returns an array over the bins of
		# sum_{|ell| in bin } C_{ell_x,ell_y}
		#and the second 
		# sum_{|ell| in bin } 1
		#so the ratio is just the mean bin value.
                if ell_weight is not None:
                    P,_ = np.histogram(self.l_abs, self.bin_edges, weights=C*ell_weight)
                    count,_ = np.histogram(self.l_abs, self.bin_edges, weights=ell_weight)
                else:
                    P,_ = np.histogram(self.l_abs, self.bin_edges, weights=C)
                    count,_ = np.histogram(self.l_abs, self.bin_edges)
		return P*(self.dx/self.N)**2/count

	def estimate(self, g1, g2, weight_EE=False, weight_BB=False, weight_EB=False):
		""" Compute the EE,BB, and EB power spectra of two 2D arrays g1 and g2."""
		#Check geometry is what we expect.
		assert g1.shape == g2.shape == (self.N, self.N)

		#Transform g1+j*g2 into Fourier and rotate into E-B
		#then separate into E and B
		EB = np.fft.ifft2(self.eb_rot * np.fft.fft2(g1 + 1j*g2))
		E = np.fft.fft2(EB.real)
		B = np.fft.fft2(EB.imag)

		#Use the internal function above to bin,
		#and account for the normalization of the FFT
		C_EE = self._bin_power(E*np.conjugate(E)) 
		C_BB = self._bin_power(B*np.conjugate(B))
		C_EB = self._bin_power(E*np.conjugate(B))

                if weight_EE or weight_BB or weight_EB:
                    # need to interpolate C_EE to values of self.l_abs
                    new_ell = np.zeros(len(self.ell)+2)
                    new_ell[1:len(self.ell)+1] = self.ell
                    new_ell[len(self.ell)+1] = 10.*max(self.ell)

                if weight_EE:
                    new_CEE = np.zeros_like(new_ell)
                    new_CEE[1:len(self.ell)+1] = np.real(C_EE)
                    new_CEE[len(self.ell)+1] = max(np.real(C_EE))
                    EE_table = galsim.LookupTable(new_ell, new_CEE)
                    ell_weight = EE_table(self.l_abs)
                    C_EE = self._bin_power(E*np.conjugate(E), ell_weight=ell_weight)
                if weight_BB:
                    new_CBB = np.zeros_like(new_ell)
                    new_CBB[1:len(self.ell)+1] = np.real(C_BB)
                    new_CBB[len(self.ell)+1] = max(np.real(C_BB))
                    BB_table = galsim.LookupTable(new_ell, new_CBB)
                    ell_weight = BB_table(self.l_abs)
                    C_BB = self._bin_power(E*np.conjugate(E), ell_weight=ell_weight)
                if weight_BB:
                    new_CBB = np.zeros_like(new_ell)
                    new_CBB[1:len(self.ell)+1] = np.real(C_BB)
                    new_CBB[len(self.ell)+1] = max(np.real(C_BB))
                    BB_table = galsim.LookupTable(new_ell, new_CBB)
                    ell_weight = BB_table(self.l_abs)
                    C_BB = self._bin_power(E*np.conjugate(E), ell_weight=ell_weight)

		#For convenience return ell (copied in case the user messes with it)
		#and the three power spectra.
		return self.ell.copy(), np.real(C_EE), np.real(C_BB), np.real(C_EB)
