;# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
;# All rights reserved.
;#
;# Redistribution and use in source and binary forms, with or without modification, are permitted
;# provided that the following conditions are met:
;#
;# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
;# and the following disclaimer.
;#
;# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
;# conditions and the following disclaimer in the documentation and/or other materials provided with
;# the distribution.
;#
;# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
;# endorse or promote products derived from this software without specific prior written permission.
;#
;# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
;# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
;# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
;# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
;# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
;# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;
; This script is the first script driven by `run_many.pro`.  It takes the SHERA catalog, and puts it
; into the format needed by GalSim.  Different filenames are assigned based on the 'per' value
; (indicating how many galaxies to put into each file).
;
PRO makecatalog_many,label,start,n,per,seed

; Define input catalog, output file.  The latter is a file that was made public on the GREAT3 data
; download page ( http://great3.projects.phys.ucl.ac.uk/leaderboard/data ); however, before that,
; another column was added in a later step of the processing.
incatw = './shera_catalog_23.5.fits'
outcat = './real_galaxy_catalog'+label+'.fits'

; How many galaxies to include?
nuse = n

; Read input catalog.
mycat = mrdfits(incatw, 1)
nlines = n_elements(mycat)
print,'Read in ',nlines,' from file ',incatw
; Note, ordering is NOT random, so let's reorder randomly.
in_seed=seed
sort_ind = randomu(seed,nlines)
mycat = mycat[sort(sort_ind)]
print,'Randomly reordered using seed',in_seed

; Select galaxies, <= nuse depending on how many galaxies are in input catalog.
mycatuse = mycat[start:start+(nuse-1 < nlines)]
ngal = n_elements(mycatuse)
print,'Using ',ngal

; Renormalize the weights from the input catalog, but enforce a maximum weight of 5 (rarely
; important but can affect a few of the very largest galaxies).
mycatuse.ps_wt = mycatuse.ps_wt < 5.
mycatuse.ps_wt = mycatuse.ps_wt/max(mycatuse.ps_wt)

indices = lindgen(ngal)

; Make new data structure.
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
           
; Populate the data structure from the input catalog.
data.ident = mycatuse.ident
data.ra = mycatuse.ra
data.dec = mycatuse.dec
data.mag = mycatuse.F814W
data.band[*] = 'F814W'
data.weight = mycatuse.ps_wt
data.pixel_scale = 0.03

; Set up the filenames where the postage stamp images are to be found, and put them in the catalog.
for i=0L,nuse-1 do begin
   nfile=i/per+1
   data[i].gal_filename = 'real_galaxy_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
   data[i].PSF_filename = 'real_galaxy_PSF_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
   data[i].gal_hdu=i-(nfile-1)*per
   data[i].PSF_hdu=i-(nfile-1)*per
endfor
   
; Write output catalog.
print,'Writing to file ',outcat
mwrfits,data,outcat,/create

END
