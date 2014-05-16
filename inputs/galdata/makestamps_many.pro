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
; This script is the second script driven by `run_many.pro`.  It creates postage stamps for all the
; galaxies in the proper format, and calculates their noise properties (variance).
;
PRO makestamps_many,label,per,varfile

; Define filenames etc.
catfile = './real_galaxy_catalog'+label+'.fits' ; input catalog from makecatalog.pro
listfile = './imgs.23.5.in'                     ; list of files and their locations

centroidval = 0

; Read in catalog file.
cat = mrdfits(catfile,1)
ncat = n_elements(cat)
print,'Read in ',ncat,' from ',catfile

; Read in image list with locations of all files.
readcol,listfile,imgfilename, imgpref,ident,format='A,A,L'
nlist = n_elements(imgfilename)
print,'Read in ',nlist,' from file ',listfile

; Open to delete previous filename.
openw,1,varfile
close,1

; Loop over the objects in the catalog.
for i=0L,ncat-1 do begin
   ; Find filenames for each ident, since my list of image file locations
   ; is not necessarily in the same order as the catalog.
   wthis = where(ident eq cat[i].ident,c)
   ; If we don't have a match in the catalog, then something is very wrong.
   if (c ne 1) then begin
      print,'Error: wrong number of matches made with ',cat[i].ident 
   endif else begin
      
      ; Read in image, PSF, segmentation map for this particular galaxy.
      ; There are a bunch of things to try, because the files might be zipped or not.
      inimgfile = imgpref[wthis[0]]+'_processed.fits'
      inimgfile2 = imgpref[wthis[0]]+'_masknoise.fits.gz'
      inpsffile = imgpref[wthis[0]]+'.psf.fits'
      insegfile = imgpref[wthis[0]]+'_seg.fits'
      if file_test(inimgfile) eq 1 then begin
         img = mrdfits(inimgfile,0)
         print,'Read from file ',inimgfile
      endif else begin
         if file_test(inimgfile2) eq 1 then begin
            img = mrdfits(inimgfile2,0)
            print,'Read from file ',inimgfile2
         endif else begin
            inimgfile = inimgfile+'.gz'
            print,'Read from file ',inimgfile
         endelse
      endelse

      if file_test(inpsffile) eq 1 then begin
         psf = mrdfits(inpsffile,0)
      endif else begin
         inpsffile = inpsffile+'.gz'
         psf = mrdfits(inpsffile,0)
      endelse
      print,'Read from file ',inpsffile

      if file_test(insegfile) eq 1 then begin
         seg = mrdfits(insegfile,0)
      endif else begin
         insegfile = insegfile+'.gz'
         seg = mrdfits(insegfile,0)
      endelse
      print,'Read from file ',insegfile

      ; Trim the postage stamp based on the segmentation image.  We do not
      ; need to save huge images, because GalSim will automatically pad
      ; them for us to the proper size needed for accurate DFT rendering.
      imgsize = size(seg)
      ximgsize = imgsize[1]
      yimgsize = imgsize[2]
      myindx = seg[floor(0.5*ximgsize), floor(0.5*yimgsize)]
      wimg = where(seg eq myindx,countin,complement=wnot)
      nimg = where(seg eq 0,countin)
      meanclip,img[nimg],mean,sig
      var= sig*sig
      

      ; We find the right region in x, y, and expand its size by factor of
      ; 2.00 (unless that goes beyond the whole postage stamp)
      tmpsize = size(img)
      centroidval = 0.5*tmpsize[1]
      arrind = array_indices(seg, wimg)
      sortxind = sort(arrind[0,*])
      minx = arrind[0,sortxind[0.05*n_elements(sortxind)]]
      maxx = arrind[0,sortxind[0.95*n_elements(sortxind)]]
      xplussize = maxx-centroidval
      xminussize = centroidval-minx
      xsize = max(xminussize, xplussize)
      minxc = round(centroidval - 2.00*xsize) > 0
      maxxc = round(centroidval + 2.00*xsize) < tmpsize[1]-1
   
      ; We check for the case when too few pixels are found in the x direction.  This could
      ; be due to centroiding issues.
      if( (maxxc-minxc) le 5L) then begin
         midx = (maxx+minx)/2
         xplussize = maxx-midx
         xminussize = midx-minx
         xsize = max(xminussize, xplussize)
         minx = round(midx - 2.0*xsize) > 0
         maxx = round(midx + 2.0*xsize) < tmpsize[1]-1
      endif else begin
         minx = minxc
         maxx = maxxc
      endelse
      
      sortyind = sort(arrind[1,*])
      miny = arrind[1,sortyind[0.05*n_elements(sortyind)]]
      maxy = arrind[1,sortyind[0.95*n_elements(sortyind)]]
      yplussize = maxy-centroidval
      yminussize = centroidval-miny
      ysize = max(yminussize, yplussize)
      
   
      minyc = round(centroidval - 2.00*ysize) > 0
      maxyc = round(centroidval + 2.00*ysize) < tmpsize[1]-1

      ; We check for the case when too few pixels are found in the y direction.  This could
      ; be due to centroiding issues.
      if( (maxyc-minyc) le 5L) then begin
         midy = (maxy+miny)/2
         yplussize = maxy-midy
         yminussize = midy-miny
         ysize = max(yminussize, yplussize)
         miny = round(midy - 2.0*ysize) > 0
         maxy = round(midy + 2.0*ysize) < tmpsize[1]-1
      endif else begin
         miny = minyc
         maxy = maxyc
      endelse


      ; Check if it's square.  If not, then make it square by
      ; adopting the size in the dimension that is larger.
      ddy = maxy + 1 - miny
      ddx = maxx + 1 - minx
      if (ddy gt ddx) then begin
         delta_size = ddy - ddx
         delta_size_1 = delta_size / 2
         delta_size_2 = delta_size - delta_size_1
         minx = minx - delta_size_1 > 0
         maxx = maxx + delta_size_2 < tmpsize[1]-1
      endif
      if (ddx gt ddy) then begin
         delta_size = ddx - ddy
         delta_size_1 = delta_size / 2
         delta_size_2 = delta_size - delta_size_1
         miny = miny - delta_size_1 > 0
         maxy = maxy + delta_size_2 < tmpsize[1]-1
      endif

      ; For the PSF, cut off the postage stamp where flux is <2000x peak.
      tmpsize = size(psf)
      centroidval = 0.5*tmpsize[1]
      wkeep = where(psf ge 0.0005*max(psf),nkeep)
      arrind = array_indices(psf,wkeep)
      minxp = min(arrind[0,*])
      maxxp = max(arrind[0,*])
      xplussize = maxxp-centroidval
      xminussize = centroidval-minxp
      xsize = max(xminussize, xplussize)
      minxp = round(centroidval - xsize) > 0
      maxxp = round(centroidval + xsize) < tmpsize[1]-1
      minyp = min(arrind[1,*])
      maxyp = max(arrind[1,*])
      yplussize = maxyp-centroidval
      yminussize = centroidval-minyp
      ysize = max(yminussize, yplussize)
      minyp = round(centroidval - ysize) > 0
      maxyp = round(centroidval + ysize) < tmpsize[1]-1
      
      ; Write to output image files.
      nfile=i/per+1
      outgalfile = 'real_galaxy_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
      outpsffile = 'real_galaxy_PSF_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
      
      
      print,'Writing to files!'
      mwrfits,img[minx:maxx,miny:maxy],outgalfile
      mwrfits,psf[minxp:maxxp,minyp:maxyp],outpsffile
      
      ; Save noise variance to text file.
      openw,1,varfile,/APPEND
      printf,1,mean,sig*sig
      close,1
                
    endelse
endfor

END
