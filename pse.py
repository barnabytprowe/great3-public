import numpy as np
from numpy import pi


def generate_eb_rotation(n, dx):
	lx=2*pi*np.fft.fftfreq(n, dx)
	LX=np.vstack([lx for i in xrange(n)])
	LY=LX.T
	L_ABS = (LX**2 + LY**2)**0.5
	L_ANG = np.arctan2(LY, LX)
	rotation = np.exp(-2j*L_ANG)
	return L_ABS, rotation


class PowerSpectrumEstimator(object):
	def __init__(self, N, sky_size_deg, nbin):
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





def power_spectrum(g1, g2, sky_size_deg, nbin):
	#Get the shape of the array.
	#We require g1 and g2 to be the same shape
	#First, on the off-chance they are lists or something, convert 
	#g1 and g2 to arrays
	g1 = np.array(g1)
	g2 = np.array(g2)

	#Check that we have been passed 2D arrays
	assert g1.ndim == g2.ndim ==2, "g1 and g2 passed to power spectrum estimation must be 2D"

	#A few more checks - that the array shapes are the same and
	#that they are square (temporary restriction till I recode this)
	nx, ny = g1.shape
	assert g1.shape==g2.shape, "To make a power spectrum you need two arrays of the same shape, g1 and g2"
	assert nx==ny, "Code only written for nx==ny"

	#Calculate the sizes and wavenumber limits in the problem
	sky_size = np.radians(sky_size_deg)
	dx = sky_size / nx
	lmin = 2*pi / sky_size
	lmax = 2*pi / dx

	#The core of the 
	g = g1 + 1j*g2
	Ell, R = generate_eb_rotation(nx, dx)
	g_f = np.fft.fft2(g)
	EB = R * g_f
	(E,B)  = (EB.real, EB.imag)

	C_EE = E**2
	C_BB = B**2
	C_EB = E*B


	#TODO: add an option to do linspace
	bin_edges = np.logspace(np.log10(lmin), np.log10(lmax), nbin+1)
	# bin_edges = np.linspace(lmin, lmax, nbin+1)
	# print bin_edges
	def bin_power(C):
		P,_ = np.histogram(Ell, bin_edges, weights=C)
		count,_ = np.histogram(Ell, bin_edges)
		# print count
		return P/count


	C_EE = bin_power(C_EE) / (nx*ny)
	C_BB = bin_power(C_BB) / (nx*ny)
	C_EB = bin_power(C_EB) / (nx*ny)

	ell = 0.5*(bin_edges[1:] + bin_edges[:-1])

	return ell, C_EE, C_BB, C_EB

#imshow(C_EE)
