# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
# endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""Script for making catalogs of galaxy fit data corresponding to a real galaxy training set used by
GalSim.  It has to collect information from several large files."""
import pyfits
import numpy as np

# Define filenames, etc.
galsim_catfile = 'real_galaxy_catalog_23.5.fits'
fit_catfiles = ['BRIGHTtotalRAW00000.26113.fits',
                'totalRAW00000.29949.fits.gz']
n_catfiles = len(fit_catfiles)
cosmos_catfile = 'lensing14.fits.gz'
out_fitfile = 'real_galaxy_catalog_23.5_fits.fits'
out_catfile = 'real_galaxy_catalog_23.5.fits'

# Read in real galaxy catalog.
galsim_cat = pyfits.getdata(galsim_catfile)
n_galsim_cat = len(galsim_cat)
print 'Read in ',n_galsim_cat,' from GalSim catalog ',galsim_catfile
galsim_ident = galsim_cat.field('ident')
# Fields: ('IDENT', 'RA', 'DEC', 'MAG', 'BAND', 'WEIGHT', 'GAL_FILENAME', 'PSF_FILENAME', 'GAL_HDU',
#          'PSF_HDU', 'PIXEL_SCALE', 'NOISE_MEAN', 'NOISE_VARIANCE')

# Read in the full COSMOS catalog.
cosmos_cat = pyfits.getdata(cosmos_catfile)
n_cosmos_cat = len(cosmos_cat)
print 'Read in ',n_cosmos_cat,' from COSMOS catalog ',cosmos_catfile
# Fields: ('IDENT', 'MAG_AUTO', 'FLUX_AUTO', 'MAGERR_AUTO', 'FLUX_RADIUS', 'FLUXERR_AUTO',
# 'KRON_RADIUS', 'MU_MAX', 'MU_CLASS', 'CLEAN', 'GOOD', 'FLAGS', 'SN', 'SN_NON_CORR', 'FWHM_IMAGE',
# 'ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
# 'PETRO_RADIUS', 'RRG_XX', 'RRG_YY', 'XXC', 'YYC', 'XYC', 'D', 'E1_R', 'E2_R', 'E1_RU', 'E2_RU',
# 'GAMMA1', 'GAMMA2', 'FOCUS_MODEL', 'IXX', 'IYY', 'IXY', 'WEIGHT_FUNCT_RADIUS', 'VAR_E1', 'VAR_E2',
# 'BOX', 'SPECZ', 'SPECZ_MARA', 'SPECZ_CLASS', 'SPECZ_ORIGIN', 'GOOD_SPECZ', 'SPECZ_BL_AGN',
# 'SPECZ_SELECTION', 'MIPS_Z', 'MIPS_LOG_L', 'MIPS_MASS', 'ZEST_TYPE', 'ZEST_BULGE',
# 'ZEST_IRREGULARITY', 'ZEST_ELONGATION', 'ZEST_GINI', 'ZEST_M20', 'ZEST_CONCENTRATION',
# 'ZEST_ASYMMETRY', 'BULGE', 'KT', 'OLD_ZPHOT', 'OLD_GOOD_ZPHOT', 'HL_KPC', 'MARA_AGN',
# 'MARA_AGN_ZPHOT', 'MARA_AGN_ZPHOT_LOW68', 'MARA_AGN_ZPHOT_HIGH68', 'KNUD_AGN', 'G1_TS', 'G2_TS',
# 'WEIGHT_TS', 'CHANDRA_GOOD', 'CHANDRA_AGN', 'CHANDRA_LX_HARD', 'CHANDRA_LX_SOFT',
# 'CHANDRA_LX_FULL', 'CHANDRA_ZETA', 'CHANDRA_ZSPEC', 'CHANDRA_CLASSZSPEC', 'CHANDRA_MODEL',
# 'CHANDRA_XMM_ID', 'XMM_GOOD', 'XMM_AGN', 'XMM_LX_HARD', 'XMM_LX_SOFT', 'XMM_LX_FULL', 'XMM_ZETA',
# 'XMM_ZSPEC', 'XMM_CLASSZSPEC', 'XMM_MODEL', 'XMM_CHANDRA_ID', 'EZE_AGN_SPECZ', 'EZE_AGN_PHOTOZ',
# 'EZE_LX', 'EZE_HR', 'EZE_SPECZ', 'EZE_PHOTOZ', 'K_CFHT', 'MATCH_CFHT', 'ERR_K_CFHT',
# 'KEVIN_MSTAR', 'KEVIN_MSTAR2', 'KEVIN_MASSERR', 'OLIV_MSTAR', 'MVIR', 'COLOR', 'TYPE2_ZPHOT_MARA',
# 'PETER_PASSIVE', 'PETER_ANGLE_PA', 'PETER_ELLIP', 'PHOTOZ_ORDER', 'PHOTOZ_NON_COMB',
# 'PHOTOZ_NON_COMB_LOW_68', 'PHOTOZ_NON_COMB_HIGH_68', 'PBZK', 'PBZK_ZPHOT', 'PBZK_MK', 'PBZK_MASS',
# 'SIGNALTONOISERATIO', 'QUASIPETROSIANAREAFRACTION', 'QUASIPETROSIANFRACTION', 'AXISRATIO', 'GINI',
# 'CONCENTRATION', 'BOB_E', 'BOB_GOOD', 'BOB_S0', 'FLUX_GIM2D', 'R_GIM2D', 'ELL_GIM2D', 'PA_GIM2D',
# 'DX_GIM2D', 'DY_GIM2D', 'SERSIC_N_GIM2D', 'R_0P5_GIM2D', 'CHI_GIM2D', 'CECILE_SL_Z',
# 'CECILE_SL_SAT', 'CECILE_SL', 'CECILE_SL_FLAG1', 'CECILE_SL_FLAG2', 'ISOLATED', 'BCG_SCALE',
# 'BCG_R200', 'ALL_P_MEM', 'ALL_GROUP_ID', 'N_GROUP_OVERLAP', 'BEST_P_MEM', 'BEST_GROUP_ID',
# 'ZPHOT', 'TYPE', 'ZPDF', 'PHOTZ_LOW_68', 'PHOTZ_HIGH_68', 'CHI', 'MODD', 'EBV', 'NBFILT',
# 'ZMINCHI2', 'ZL68_MINCHI2', 'ZU68_MINCHI2', 'ZP2', 'CHI2', 'NUV', 'U', 'SUBARU_R', 'SUBARU_I',
# 'J_WFCAM', 'K_WIRCAM', 'M36', 'DNUV', 'DU', 'DJ_WFCAM', 'DK_WIRCAM', 'DM36', 'AUTO_OFFSET',
# 'AUTO_FLAG', 'MNUV', 'MU', 'MB', 'MV', 'MG', 'MR', 'MI', 'MJ', 'MK', 'MNUV_MR', 'SFR_MED',
# 'STR_INF', 'SFR_SUP', 'SSFR_MED', 'SSFR_INF', 'SSFR_SUP', 'MATCH_S', 'MASK_S', 'GOOD_ZPHOT_LENS',
# 'GOOD_ZPHOT_SOURCE')
# That's a lot of info, so let's just pick out the things we care about: galaxy identifier, apparent
# magnitude, size, photo-z.
cos_ident = cosmos_cat.field('ident')
cos_mag_auto = cosmos_cat.field('mag_auto')
cos_flux_rad = cosmos_cat.field('flux_radius')
cos_zphot = cosmos_cat.field('zphot')

# Read in catalogs with fit parameters from Lackner & Gunn.
print "Reading in catalogs of fit parameters"
n_fit_tot = 0
for i_cat in range(n_catfiles):
    # Get this catalog
    dat = pyfits.getdata(fit_catfiles[i_cat])
    n = len(dat)
    print "Read in ",n," fit results from file ",fit_catfiles[i_cat]
    # Just extract the columns we want, and append to previous if i_cat!=0.
    if i_cat == 0:
        fit_ident = dat.field('ident')
        fit_sersicfit = dat.field('sersicfit')
        fit_bulgefit = dat.field('bulgefit')
        fit_status = dat.field('mpfit_status')
        fit_mag_auto = dat.field('mag_auto')
        fit_mad_s = dat.field('mad_sersic_mask')
        fit_mad_b = dat.field('mad_dvcb_mask')
        fit_dvc_btt = dat.field('dvc_btt')
    if i_cat > 0:
        fit_ident = np.append(fit_ident, dat.field('galid'))
        fit_sersicfit = np.append(fit_sersicfit, dat.field('sersicfit'), axis=0)
        fit_bulgefit = np.append(fit_bulgefit, dat.field('bulgefit'), axis=0)
        fit_status = np.append(fit_status, dat.field('mpfit_status'), axis=0)
        fit_mag_auto = np.append(fit_mag_auto, np.zeros_like(dat.field('galid')), axis=0)
        fit_mad_s = np.append(fit_mad_s, dat.field('mad_sersic_mask'), axis=0)
        fit_mad_b = np.append(fit_mad_b, dat.field('mad_dvcb_mask'), axis=0)
        fit_dvc_btt = np.append(fit_dvc_btt, dat.field('dvc_btt'), axis=0)

    # Increment counter.
    n_fit_tot += n

    # Unfortunately, the files do not have the same column names.   Here are their contents -
    # Fields in first file: ('IDENT', 'MAG_AUTO', 'FLUX_AUTO', 'MAGERR_AUTO', 'FLUX_RADIUS',
    # 'FLUXERR_AUTO', 'KRON_RADIUS', 'MU_MAX', 'MU_CLASS', 'CLEAN', 'GOOD', 'FLAGS', 'SN',
    # 'SN_NON_CORR', 'FWHM_IMAGE', 'ALPHA_J2000', 'DELTA_J2000', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE',
    # 'B_IMAGE', 'THETA_IMAGE', 'PETRO_RADIUS', 'D', 'E1_R', 'E2_R', 'E1_RU', 'E2_RU', 'GAMMA1',
    # 'GAMMA2', 'FOCUS_MODEL', 'IXX', 'IYY', 'IXY', 'WEIGHT_FUNCT_RADIUS', 'VAR_E1', 'VAR_E2',
    # 'BOX', 'SPECZ', 'SPECZ_MARA', 'SPECZ_CLASS', 'SPECZ_ORIGIN', 'GOOD_SPECZ', 'SPECZ_BL_AGN',
    # 'FORS2_OBJECT_FLAG', 'MIPS_Z', 'MIPS_LOG_L', 'MIPS_MASS', 'ZEST_TYPE', 'ZEST_BULGE',
    # 'ZEST_IRREGULARITY', 'ZEST_ELONGATION', 'ZEST_GINI', 'ZEST_M20', 'ZEST_CONCENTRATION',
    # 'ZEST_ASYMMETRY', 'BULGE', 'KT', 'OLD_ZPHOT', 'OLD_GOOD_ZPHOT', 'HL_KPC', 'CHANDRA_GOOD',
    # 'CHANDRA_AGN', 'CHANDRA_LX_HARD', 'CHANDRA_LX_SOFT', 'CHANDRA_LX_FULL', 'CHANDRA_ZETA',
    # 'CHANDRA_ZSPEC', 'CHANDRA_CLASSZSPEC', 'CHANDRA_MODEL', 'CHANDRA_TYPE', 'CHANDRA_LUSSO_MASS',
    # 'XMM_GOOD', 'XMM_AGN', 'XMM_LX_HARD', 'XMM_LX_SOFT', 'XMM_ZETA', 'XMM_ZSPEC',
    # 'XMM_CLASSZSPEC', 'XMM_MODEL', 'XMM_TYPE', 'XMM_LUSSO_MASS', 'AGN_GOOD', 'AGN_Z', 'AGN_TYPE',
    # 'AGN_LX', 'AGN_LX_SOFT', 'AGN_LX_HARD', 'AGN_LUSSO_MASS', 'BOSS_LRG', 'K_CFHT', 'MATCH_CFHT',
    # 'ERR_K_CFHT', 'KEVIN_MSTAR', 'KEVIN_MSTAR2', 'KEVIN_MASSERR', 'KEVIN_QUENCH_FLAG', 'MVIR',
    # 'TYPE2_ZPHOT_MARA', 'PHOTOZ_ORDER', 'PHOTOZ_NON_COMB', 'PHOTOZ_NON_COMB_LOW_68',
    # 'PHOTOZ_NON_COMB_HIGH_68', 'FLUX_GIM2D', 'R_GIM2D', 'ELL_GIM2D', 'PA_GIM2D', 'DX_GIM2D',
    # 'DY_GIM2D', 'SERSIC_N_GIM2D', 'R_0P5_GIM2D', 'CHI_GIM2D', 'CECILE_SL_Z', 'CECILE_SL_SAT',
    # 'CECILE_SL', 'CECILE_SL_FLAG1', 'CECILE_SL_FLAG2', 'GROUP_PROJECTION_MMGG',
    # 'GROUP_PROJECTION_MMGG_SPECZ', 'MMGG_SCALE', 'P_MEM_BEST', 'GROUP_ID_BEST', 'GROUP_FLAG_BEST',
    # 'P_MEM_ALL', 'GROUP_ID_ALL', 'GROUP_FLAG_ALL', 'DIST_BCG_R200', 'MMGG_SCALE_SPECZ',
    # 'P_MEM_BEST_SPECZ', 'GROUP_ID_BEST_SPECZ', 'GROUP_FLAG_BEST_SPECZ', 'P_MEM_ALL_SPECZ',
    # 'GROUP_ID_ALL_SPECZ', 'GROUP_FLAG_ALL_SPECZ', 'DIST_BCG_R200_SPECZ', 'ZPHOT', 'TYPE', 'ZPDF',
    # 'PHOTZ_LOW_68', 'PHOTZ_HIGH_68', 'CHI', 'MODD', 'EBV', 'NBFILT', 'ZMINCHI2', 'ZL68_MINCHI2',
    # 'ZU68_MINCHI2', 'ZP2', 'CHI2', 'NUV', 'U', 'B', 'SUBARU_R', 'SUBARU_I', 'J_WFCAM', 'K_WIRCAM',
    # 'M36', 'DNUV', 'DU', 'DJ_WFCAM', 'DK_WIRCAM', 'DM36', 'AUTO_OFFSET', 'AUTO_FLAG', 'MNUV',
    # 'MU', 'MB', 'MV', 'MG', 'MR', 'MI', 'MJ', 'MK', 'MNUV_MR', 'SFR_MED', 'STR_INF', 'SFR_SUP',
    # 'SSFR_MED', 'SSFR_INF', 'SSFR_SUP', 'MATCH_S', 'MASK_S', 'GOOD_ZPHOT_LENS',
    # 'GOOD_ZPHOT_SOURCE', 'RA', 'DEC', 'GALID', 'BULGEFIT', 'DISKFIT', 'SERSICFIT', 'CHISQ_BULGE',
    # 'CHISQ_DISK', 'CHISQ_SERSIC', 'COVAR_BULGE', 'COVAR_DISK', 'COVAR_SERSIC', 'PERR_BULGE',
    # 'PERR_DISK', 'PERR_SERSIC', 'MPFIT_STATUS', 'DOF_BULGE', 'DOF_DISK', 'DOF_SERSIC', 'DOF_DVC',
    # 'DOF_EXP', 'EXPFIT', 'DVCFIT', 'CHISQ_EXP', 'CHISQ_DVC', 'PERR_EXP', 'PERR_DVC', 'COVAR_EXP',
    # 'COVAR_DVC', 'FRACDEV', 'XCROP', 'YCROP', 'XLEN', 'YLEN', 'DVC_BTT', 'EXP_BTT', 'MAD_SKY',
    # 'MAD_SERSIC', 'MAD_SERSIC_MASK', 'MAD_DVCB', 'MAD_DVCB_MASK', 'MAD_EXPB', 'MAD_EXPB_MASK',
    # 'MAD_EXP', 'MAD_EXP_MASK', 'MAD_DVC', 'MAD_DVC_MASK', 'CHISQ_BULGE_MASK', 'CHISQ_DISK_MASK',
    # 'CHISQ_EXP_MASK', 'CHISQ_SERSIC_MASK', 'CHISQ_DVC_MASK', 'DOF_BULGE_MASK', 'DOF_DISK_MASK',
    # 'DOF_EXP_MASK', 'DOF_SERSIC_MASK', 'DOF_DVC_MASK', 'SN_REFF_SERSIC', 'SKY_SERSIC',
    # 'SKY_SERSIC_ERR', 'SKY_SERSIC_COVAR', 'DVC_BTT_ERR', 'EXP_BTT_ERR')

print "Read in ",n_fit_tot," from ",n_catfiles," fit files"

print "Making correspondence between IDENT values for all inputs"
cos_ind = np.zeros_like(galsim_ident)
fit_ind = np.zeros_like(galsim_ident)
cos_ident_list = list(cos_ident)
fit_ident_list = list(fit_ident)
n_fail_cos = 0
n_fail_fit = 0
for i in range(n_galsim_cat):
    if i % 1000 == 0:
        print "... object ",i
    if galsim_ident[i] in cos_ident_list:
        cos_ind[i] = cos_ident_list.index(galsim_ident[i])
    else:
        cos_ind[i] = -1
        n_fail_cos += 1
    if galsim_ident[i] in fit_ident_list:
        fit_ind[i] = fit_ident_list.index(galsim_ident[i])
    else:
        fit_ind[i] = -1
        n_fail_fit += 1
print "Number of match failures for COSMOS, fits: ",n_fail_cos, n_fail_fit

print "Rearranging arrays into proper order"
use_ind = (fit_ind >= 0) & (cos_ind >= 0)
out_ident = galsim_ident[use_ind]
print "Actually using ",len(out_ident)
out_mag_auto = cos_mag_auto[cos_ind[use_ind]]
out_flux_rad = cos_flux_rad[cos_ind[use_ind]]
out_zphot = cos_zphot[cos_ind[use_ind]]
test_mag_auto = fit_mag_auto[fit_ind[use_ind]]
print 'Mag auto test:'
print out_mag_auto[0:9]
print test_mag_auto[0:9]

# Rearrange the FIT arrays with fit quantities in the same order as galsim_ident.
out_sersicfit = fit_sersicfit[fit_ind[use_ind],:]
out_bulgefit = fit_bulgefit[fit_ind[use_ind],:]
out_fit_status = fit_status[fit_ind[use_ind],:]
out_fit_mad_s = fit_mad_s[fit_ind[use_ind],:]
out_fit_mad_b = fit_mad_b[fit_ind[use_ind],:]
out_fit_dvc_btt = fit_dvc_btt[fit_ind[use_ind],:]

# Make output data structure with IDENT, photo-z, magnitude, flux_radius, SERSICFIT, BULGEFIT, fit
# status.  SERSICFIT and BULGEFIT are actually arrays of fit parameters from single Sersic fits and
# two-component fits, respectively.
tbhdu = pyfits.new_table(pyfits.ColDefs([pyfits.Column(name='IDENT',
                                                       format='J',
                                                       array=out_ident),
                                         pyfits.Column(name='mag_auto',
                                                       format='D',
                                                       array=out_mag_auto),
                                         pyfits.Column(name='flux_radius',
                                                       format='D',
                                                       array=out_flux_rad),
                                         pyfits.Column(name='zphot',
                                                       format='D',
                                                       array=out_zphot),
                                         pyfits.Column(name='sersicfit',
                                                       format='8D',
                                                       array=out_sersicfit),
                                         pyfits.Column(name='bulgefit',
                                                       format='16D',
                                                       array=out_bulgefit),
                                         pyfits.Column(name='fit_status',
                                                       format='5J',
                                                       array=out_fit_status),
                                         pyfits.Column(name='fit_mad_s',
                                                       format='D',
                                                       array=out_fit_mad_s),
                                         pyfits.Column(name='fit_mad_b',
                                                       format='D',
                                                       array=out_fit_mad_b),
                                         pyfits.Column(name='fit_dvc_btt',
                                                       format='D',
                                                       array=out_fit_dvc_btt)]
                                        ))

# Write outputs.
print "Writing to file ",out_fitfile
tbhdu.writeto(out_fitfile, clobber=True)

# Write new subset of catalog file.
print "Re-writing to file ",out_catfile
galsim_cat = pyfits.BinTableHDU(galsim_cat[use_ind])
galsim_cat.writeto(out_catfile, clobber=True)
