PRO getpk

nb=200

fid = set_fiducial(cosmo_in={H:0.70d, omega_m:0.25d, omega_l:0.75d, w0:-1.d, n:0.96d, sigma8:0.8d, curv:0}, calc_in={fit_nl:2, fit_TK:1, verbose:2, speed:0, n_lbin:nb}, expt_in={sv1_N_ZBIN:1, sv1_Z_MED:0.75})
sv = mk_survey(fid, 'sv1')
cosmo = mk_cosmo(fid)
cl = mk_cl_tomo(fid, cosmo, sv)
forprint,textout='tables/cosmo-fid.zmed0.75.out',cl.l,cl.cl,format='(2E13.6)',/nocomment

fid = set_fiducial(cosmo_in={H:0.70d, omega_m:0.25d, omega_l:0.75d, w0:-1.d, n:0.96d, sigma8:0.8d, curv:0}, calc_in={fit_nl:2, fit_TK:1, verbose:2, speed:0, n_lbin:nb}, expt_in={sv1_N_ZBIN:1, sv1_Z_MED:1.00})
sv = mk_survey(fid, 'sv1')
cosmo = mk_cosmo(fid)
cl = mk_cl_tomo(fid, cosmo, sv)
forprint,textout='tables/cosmo-fid.zmed1.00.out',cl.l,cl.cl,format='(2E13.6)',/nocomment

fid = set_fiducial(cosmo_in={H:0.70d, omega_m:0.25d, omega_l:0.75d, w0:-1.d, n:0.96d, sigma8:0.8d, curv:0}, calc_in={fit_nl:2, fit_TK:1, verbose:2, speed:0, n_lbin:nb}, expt_in={sv1_N_ZBIN:1, sv1_Z_MED:1.25})
sv = mk_survey(fid, 'sv1')
cosmo = mk_cosmo(fid)
cl = mk_cl_tomo(fid, cosmo, sv)
forprint,textout='tables/cosmo-fid.zmed1.25.out',cl.l,cl.cl,format='(2E13.6)',/nocomment

END
