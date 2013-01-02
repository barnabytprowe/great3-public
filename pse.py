import numpy as np
from numpy import pi


class PowerSpectrumEstimator(object):
	"""This class stores all the data used in power spectrum estimation that 
	is fixed with the geometry of the problem - the binning and spin weighting factors.

	The only public method is estimate, which you call with 2D g1 and g2 arrays.
	"""
	def __init__(self, N, sky_size_deg, nbin):
		"""Create a PSE objection with the number of pixels along each side, N, 
		the total sky width in degrees, and the number of evenly-log-spaced ell bins to use."""
		self.nx=N
		self.ny=N
		self.sky_size = np.radians(sky_size_deg)
		self.dx = self.sky_size / nx
		self.lmin = 2*pi / sky_size
		self.lmax = 2*pi / dx
		self.bin_edges = np.logspace(np.log10(self.lmin), np.log10(self.lmax), nbin+1)
		self.ell = 0.5*(self.bin_edges[1:] + self.bin_edges[:-1])

		self.l_abs, self.eb_rot = self._generate_eb_rotation()

	def _generate_eb_rotation(self):
		lx=2*pi*np.fft.fftfreq(self.N, self.dx)
		LX=np.vstack([lx for i in xrange(self.N)])
		LY=LX.T
		L_ABS = (LX**2 + LY**2)**0.5
		L_ANG = np.arctan2(LY, LX)
		rotation = np.exp(-2j*L_ANG)
		return L_ABS, rotation

	def _bin_power(self, C):
		P,_ = np.histogram(self.l_abs, self.bin_edges, weights=C)
		count,_ = np.histogram(self.l_abs, self.bin_edges)
		return P/count

	def estimate(self, g1, g2):
		assert g1.shape == g2.shape == (self.N, self.N)

		EB = self.eb_rot * np.fft.fft2(g1 + 1j*g2)
		(E,B)  = (EB.real, EB.imag)

		C_EE = bin_power(E**2) / (nx*ny)
		C_BB = bin_power(B**2) / (nx*ny)
		C_EB = bin_power(E*B ) / (nx*ny)

		return self.ell.copy(), C_EE, C_BB, C_EB




