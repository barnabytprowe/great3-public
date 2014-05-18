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
; This script is the third script driven by `run_many.pro`.  It adds the noise variances to the
; catalogs.
;
PRO makecatalog_many_var,label,start,n,per,varfile,seed

; Define input (and output) catalog filenames, and a temp file for intermediate outputs.  This
; should be checked manually before copying it over the input file.
outcat = './real_galaxy_catalog'+label+'.fits'
tmpoutcat = 'foo.fits'

; Read input / output catalog.
mycat = mrdfits(outcat, 1)
nlines = n_elements(mycat)
print,'Read in ',nlines,' from file ',outcat

; Read file with variances.
OPENR,1,varfile
H=FLTARR(2,nlines)
READF,1,H
CLOSE,1
mean = H(0,*)
var = H(1,*)

; Make new data structure with room for noise variances.
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
           
; Populate the data structure from the input catalogs.
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

; Write output catalog.
print,'Writing to file ',tmpoutcat
mwrfits,data,tmpoutcat,/create

END
