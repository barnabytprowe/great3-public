PRO makecatalog_many_var,label,start,n,per,varfile,seed

;; program used by Rachel to generate the catalogs for the first N
;; real galaxies.  Assign different files to each "per" number of
;; of galaxies.  Read in a precomputed mean and variance for
;; for each file

; define input catalog, output file
outcat = './real_galaxy_catalog'+label+'.fits'
tmpoutcat = 'foo.fits'

; read input / output catalog
mycat = mrdfits(outcat, 1)
nlines = n_elements(mycat)
print,'Read in ',nlines,' from file ',outcat

; 
OPENR,1,varfile
H=FLTARR(2,nlines)
READF,1,H
CLOSE,1
mean = H(0,*)
var = H(1,*)

; make new data structure
datastr = {ident: 0L, $
           RA: 0.D, $
           DEC: 0.D, $
           mag: 0.D, $
           band: '', $
           weight: 0.D, $
           gal_filename: '', $
           PSF_filename: '', $
           gal_hdu: 0L, $
           PSF_hdu: 0L, $
           pixel_scale: 0.D, $; arcsec
           noise_mean: 0.D, $
           noise_variance: 0.D, $
           noise_filename: ''}
data = replicate(datastr, nlines)
           
; populate the data structure from the input catalog
data.ident = mycat.ident
data.ra = mycat.ra
data.dec = mycat.dec
data.mag = mycat.mag
data.band = mycat.band
data.weight = mycat.weight
data.pixel_scale = mycat.pixel_scale
data.gal_filename = mycat.gal_filename
data.PSF_filename = mycat.PSF_filename
data.gal_hdu = mycat.gal_hdu
data.PSF_hdu = mycat.PSF_hdu

for i=0L,nlines-1 do begin
   data[i].noise_mean = mean[i]
   data[i].noise_variance = var[i]
   data[i].noise_filename = 'acs_I_unrot_sci_20_cf.fits'
endfor

;data.gal_hdu = indices
;data.PSF_hdu = indices
           
; write output catalog
print,'Writing to file ',tmpoutcat
mwrfits,data,tmpoutcat,/create

END
