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
; This script was used to run was used to drive several scripts, in sequence, to generate the
; postage stamps and a GalSim-style catalog representing the galaxy sample.  It was used with the
; following calling sequence in IDL:
;
;    run_many,'_23.5',0,56062,1000,10
;
; This calling sequence tells the script about the size of the galaxy sample, how many postage
; stamps to include in a given image file (1000), and a random seed to use when randomly reordering
; the sample (10).
;
PRO run_many,label,start,n,per,seed

; Make an initial catalog in the format needed for GalSim.
makecatalog_many,label,start,n,per,seed

; Create postage stamps for all the galaxies, and calculate their noise properties.
makestamps_many,label,per,'real_galaxy_var.txt'

; Take the catalog from the initial step, and create a new catalog with noise properties included.
; There's some amount of repetition here, but this step is not the rate-limiting step anyway.
makecatalog_many_var,label,start,n,per,'real_galaxy_var.txt',seed

end

