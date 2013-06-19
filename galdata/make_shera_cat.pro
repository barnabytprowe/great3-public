PRO make_shera_cat

; set up input and output file names
outfile = 'shera_catalog_23.5.fits'
in1 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.dat'
in4 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.weightfac.norm.out'
in5 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.useflags.out'

; read in input data: decide which fields we need!
readcol,in1,ident,f814w,ra,dec,exti,format='L,X,X,X,F,D,D,X,X,X,X,X,X,X,X,X,F'
;; note: f814w is not extinction corrected
f814wi = f814w-0.03
f814wcorr = f814w-exti

readcol,in4,wt,wtmdcstat,wtmd,wtcstat,format='F,F,F,F'

readcol,in5,mdflag,format='X,X,I'

ngal = n_elements(ra)
print,'Read in info about ',ngal,' galaxies with postage stamps'

; set up structure for data: use something like
datastr = {IDENT: 0L, $
           RA: 0.D, $
           DEC: 0.D, $
           F814W: 0., $
           PS_WT: 0., $
           USE_FLAG: 0L}
data =  replicate(datastr, ngal)

data.ident = ident
data.ra = ra
data.dec = dec
data.f814w = f814w
data.ps_wt = wt
data.use_flag = mdflag

; write to FITS file
print,'Writing to file ',outfile
mwrfits,data,outfile,/create

END
