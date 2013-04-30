import numpy as np
import matplotlib.pyplot as plt

def run_corr2_ascii(x,y,e1,e2,min_sep,max_sep,nbins):
    f = open('temp.cat','w')
    for (i,j), value in np.ndenumerate(x):
        xij = value
        yij = y[i,j]
        g1ij = e1[i,j]
        g2ij = e2[i,j]
        f.write('%e  %e  %e  %e\n'%(xij,yij,g1ij,g2ij))
    f.close()
    print min_sep, max_sep, nbins
    subprocess.Popen(['/Users/rmandelb/svn/mjarvis-read-only/corr2','corr2.params',
                      'min_sep=%f'%min_sep,'max_sep=%f'%max_sep,'nbins=%f'%nbins]).wait()
    results = np.loadtxt('temp.e2')
    os.remove('temp.cat')
    os.remove('temp.e2')
    return results

fig = plt.figure()
ax = fig.add_subplot(111)
# define values for theta_0 (which had been pre-determined by me for mathematica calculations and
# used to make filenames, meaning this part of the code can't be changed without changing that other
# stuff too).
min_theta0 = 360. # 0.1 degree
max_theta0 = 7200. # 2 degrees
ntheta0 = 11
theta0 = np.logspace(np.log10(min_theta0),np.log10(max_theta0),ntheta0)

# loop over values of theta_0, read in P(k), plot
for t0 in theta0:
    strval = str(int(round(t0)))
    infile = 'pk_math/Pk'+strval+'.dat'
    pk_dat = np.loadtxt(infile).transpose()
    k = pk_dat[0]
    pk = 1.e-4*pk_dat[1] # put in the normalization here
    ax.plot(k, pk, label='theta0='+strval)
# and just for fun, let's plot a cosmological PS.  For that, we have to convert the iCosmo one to
# arcsec instead of radians.
pk_file = '/Users/rmandelb/great3/GalSim/examples/data/cosmo-fid.zmed1.00.out'
pk_data = np.loadtxt(pk_file).transpose()
rad_to_arcsec = 3600.*180./np.pi
k_arcsec = pk_data[0]/rad_to_arcsec
pk_arcsec = pk_data[1]*(rad_to_arcsec**2)
ax.plot(k_arcsec, pk_arcsec, label='<gg>, z=1')
ell_20 = 20./rad_to_arcsec
ell_2000 = 2000./rad_to_arcsec
ax.plot((ell_20, ell_20), (1.e-6,1.e5))
ax.plot((ell_2000, ell_2000), (1.e-6,1.e5))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('k [1/arcsec]')
ax.set_ylabel('P(k) [arcsec^2]')
plt.legend()
print "Saving plot to file"
plt.savefig('diagnostics/test_pk_math.jpg')

## Code below hasn't been modified to work for the mathematica P(k)!!
# the shear part is expensive, so let's not do it every time, particularly if I just want to mess
# with the P(k) calculation.
do_shear = False
# also make an option to just use some of the grid, to save time
use_subgrid = True
if do_shear:
    # Take the outputs and run Mike's corr2 program to check that the outputs correspond to
    # inputs.  For now, let's just use the theta0=100 arcsec case.  Use 0.5*P but give that to both P_E
    # and P_B.
    print "getting shears using GalSim"
    # some kludginess to avoid neg numbers, but I did check that this is not very bad
    new_p = 0.5*p1-np.min(0.5*p1)+1.e-12
    tab_pk = galsim.LookupTable(test_k, new_p, x_log=True, f_log=True)
    ps = galsim.PowerSpectrum(tab_pk, tab_pk, units=galsim.arcsec)
    # Note, typically this returns shear, not e.  For nearly round things, e~2*shear.  But for the
    # atmospheric PSF stuff, I gave the code amplitudes corresponding to some ellipticity variance, so
    # the results should correspond to e1/e2 in amplitude as well.  The alternative approach would have
    # been to use shear variances (factor of 4 lower) and then to take the e that come out and multiply
    # by 2.  This seemed silly so I didn't do it.  However, it does mean that the kappa fluctuations are
    # too high by a factor of 4, so I must reduce them.
    e1, e2, kappa = ps.buildGrid(grid_spacing=dtheta_arcsec*subsample,
                                 ngrid=ngrid/subsample,
                                 get_convergence=True)
    kappa /= 4.
    if not use_subgrid:
        grid_range = (dtheta_arcsec*subsample) * np.arange(ngrid/subsample)
        x, y = np.meshgrid(grid_range, grid_range)
        results = run_corr2_ascii(x,y,e1,e2,dtheta_arcsec*subsample,np.sqrt(2)*dtheta_arcsec*ngrid,ngrid/subsample)
        theta_corr2 = results[:,1]
        xip_corr2 = results[:,2]
        xim_corr2 = results[:,3]    
    else:
        n_iter = 0
        grid_range = (dtheta_arcsec*subsample) * np.arange(ngrid/(subsample*oversample))
        x, y = np.meshgrid(grid_range, grid_range)
        for x_ind in range(oversample):
            x_start = x_ind*ngrid/(subsample*oversample)
            x_end = (x_ind+1)*ngrid/(subsample*oversample)
            for y_ind in range(oversample):
                print "  Subsample ",n_iter
                y_start = y_ind*ngrid/(subsample*oversample)
                y_end = (y_ind+1)*ngrid/(subsample*oversample)
                e1_sub = e1[x_start:x_end]
                e2_sub = e2[y_start:y_end]
                results = run_corr2_ascii(x,y,e1_sub,e2_sub,dtheta_arcsec*subsample,np.sqrt(2)*dtheta_arcsec*ngrid/oversample,ngrid/(subsample*oversample))
                if n_iter == 0:
                    theta_corr2 = results[:,1]
                    xip_corr2 = results[:,2]
                    xim_corr2 = results[:,3]
                else:
                    xip_corr2 += results[:,2]
                    xim_corr2 += results[:,3]
                n_iter += 1
        xip_corr2 /= oversample**2
        xim_corr2 /= oversample**2

    print "Plotting ellip corr func"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(theta_corr2, xip_corr2, label='observed xi+')
    ax.plot(theta_corr2, xim_corr2, label='observed xi-')
    ax.plot(theta_corr2, 1.e-4/(1+theta_corr2/100.), label='input xi+')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('theta [arcsec]')
    ax.set_ylabel('xi')
    plt.legend()
    print "Saving plot to file"
    plt.savefig('diagnostics/test_atmos_pk_xi.jpg')
