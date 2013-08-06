PRO makecatalog_many,label,start,n,per,seed

;; program used by Rachel to generate the catalogs for the first N
;; real galaxies.  Different filenames are assigned based on
;; the "per" value.

; define input catalog, output file
incatw = './shera_catalog_23.5.fits'
outcat = './real_galaxy_catalog'+label+'.fits'

; how many galaxies to include?
nuse = n

; read input catalog: note, ordering is NOT random, so let's reorder randomly
mycat = mrdfits(incatw, 1)
nlines = n_elements(mycat)
print,'Read in ',nlines,' from file ',incatw
in_seed=seed
sort_ind = randomu(seed,nlines)
mycat = mycat[sort(sort_ind)]
print,'Randomly reordered using seed',in_seed

; select galaxies, <= nuse depending on how many galaxies are in input catalog
mycatuse = mycat[start:start+(nuse-1 < nlines)]
ngal = n_elements(mycatuse)
print,'Using ',ngal

; renormalize the weights from the input catalog, but let's not
; be too crazy here: use max weight of 5.
mycatuse.ps_wt = mycatuse.ps_wt < 5.
mycatuse.ps_wt = mycatuse.ps_wt/max(mycatuse.ps_wt)

indices = lindgen(ngal)

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
           pixel_scale: 0.D} ; arcsec
data = replicate(datastr, ngal)
           
; populate the data structure from the input catalog
data.ident = mycatuse.ident
data.ra = mycatuse.ra
data.dec = mycatuse.dec
data.mag = mycatuse.F814W
data.band[*] = 'F814W'
data.weight = mycatuse.ps_wt
data.pixel_scale = 0.03

for i=0L,nuse-1 do begin
   nfile=i/per+1
   data[i].gal_filename = 'real_galaxy_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
   data[i].PSF_filename = 'real_galaxy_PSF_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
   data[i].gal_hdu=i-(nfile-1)*per
   data[i].PSF_hdu=i-(nfile-1)*per
endfor
   
; write output catalog
print,'Writing to file ',outcat
mwrfits,data,outcat,/create

END
