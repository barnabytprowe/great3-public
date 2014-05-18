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
; This script was used to put together a catalog listing all galaxies in the F814W<23.5 sample, in a
; format that can be used by GalSim.  It takes some inputs that are not in this repository, so it
; cannot be run here, but is made available to illustrate the steps in the process of preparing
; the sample.
;
PRO make_shera_cat

; Set up input and output file names.
outfile = 'shera_catalog_23.5.fits'
in1 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.dat'
in4 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.weightfac.norm.out'
in5 = '/scr/chimera0/rmandelb/cosmos/Aug2010/images.23.5.useflags.out'

; Read in input data: decide which fields we need!
readcol,in1,ident,f814w,ra,dec,exti,format='L,X,X,X,F,D,D,X,X,X,X,X,X,X,X,X,F'
; Note: F814W (COSMOS MAG_AUTO) is not extinction corrected. Below, we define a version that is
; corrected to i band, and a version that is extinction corrected.
f814wi = f814w-0.03
f814wcorr = f814w-exti

; This file includes weights for the probability of the galaxy making it into the catalog given the
; bias against large galaxies (since we throw out things within a certain distance of CCD edges).
readcol,in4,wt,wtmdcstat,wtmd,wtcstat,format='F,F,F,F'
; This file includes flags for issues with the images.
readcol,in5,mdflag,format='X,X,I'

ngal = n_elements(ra)
print,'Read in info about ',ngal,' galaxies with postage stamps'

; Set up structure for data.
datastr = {IDENT: 0L, $
           RA: 0.D, $
           DEC: 0.D, $
           F814W: 0., $
           PS_WT: 0., $
           USE_FLAG: 0L}
data =  replicate(datastr, ngal)

; Populate all the fields.
data.ident = ident
data.ra = ra
data.dec = dec
data.f814w = f814w
data.ps_wt = wt
data.use_flag = mdflag

; Write to FITS file.
print,'Writing to file ',outfile
mwrfits,data,outfile,/create

END
