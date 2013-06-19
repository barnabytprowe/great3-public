PRO makestamps_many,label,per,varfile

;; IDL script used by Rachel to make the postage stamps for the 
;; galaxies from the catalog provided.  Also need when to create
;; undesirable in the long term (how much padding is needed, etc.) but
;; is good enough for a start.

; define filenames etc.
catfile = './real_galaxy_catalog'+label+'.fits' ; input catalog from makecatalog.pro
listfile = './imgs.23.5.in' ; list of files and their locations

; 
centroidval = 0

; read in catalog file
cat = mrdfits(catfile,1)
ncat = n_elements(cat)
print,'Read in ',ncat,' from ',catfile

; read in image list with locations of all files
readcol,listfile,imgfilename, imgpref,ident,format='A,A,L'
nlist = n_elements(imgfilename)
print,'Read in ',nlist,' from file ',listfile

; Open to delete previous filename
openw,1,varfile
close,1

; loop over the number of objects in the catalog
for i=0L,ncat-1 do begin
; find filenames for each ident, since my list of image file locations
; is not necessarily in the same order as the catalog
   wthis = where(ident eq cat[i].ident,c)
   if (c ne 1) then begin
      print,'Error: wrong number of matches made with ',cat[i].ident 
   endif else begin
      
; read in image, PSF, seg
      inimgfile = imgpref[wthis[0]]+'_processed.fits'
      inimgfile2 = imgpref[wthis[0]]+_'masknoise.fits.gz'
      inpsffile = imgpref[wthis[0]]+'.psf.fits'
      insegfile = imgpref[wthis[0]]+'_seg.fits'
      if file_test(inimgfile) eq 1 then begin
         img = mrdfits(inimgfile,0)
         print,'Read from file ',inimgfile
      endif else if file_test(inimgfile2) begin
         img = mrdfits(inimgfile2,0)
         print,'Read from file ',inimgfile2
      endif else begin
         inimgfile = inimgfile+'.gz'
         print,'Read from file ',inimgfile
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

; trim the postage stamp based on the segmentation image
                ;; first, find pixels belonging to central object
   imgsize = size(seg)
   ximgsize = imgsize[1]
   yimgsize = imgsize[2]
   myindx = seg[floor(0.5*ximgsize), floor(0.5*yimgsize)]
   wimg = where(seg eq myindx,countin,complement=wnot)
   nimg = where(seg eq 0,countin)
   meanclip,img[nimg],mean,sig
   var= sig*sig


   ;; find that region in x, y, and expand its size by factor of
   ;; 1.25 (unless that goes beyond the whole postage stamp)
   tmpsize = size(img)
   centroidval = 0.5*tmpsize[1]
   arrind = array_indices(seg, wimg)
   sortxind = sort(arrind[0,*])
   minx = arrind[0,sortxind[0.05*n_elements(sortxind)]]
   maxx = arrind[0,sortxind[0.95*n_elements(sortxind)]]
   xplussize = maxx-centroidval
   xminussize = centroidval-minx
   xsize = max(xminussize, xplussize)
   
   
   minxc = round(centroidval - 1.25*xsize) > 0
   maxxc = round(centroidval + 1.25*xsize) < tmpsize[1]-1

   ;; check for case when too few pixels are found.  Could be due
   ;; to centroid being displaced.  just use midpoint of image
   
   if( (maxxc-minxc) le 5L) then begin
      midx = (maxx+minx)/2
      xplussize = maxx-midx
      xminussize = midx-minx
      xsize = max(xminussize, xplussize)
      minx = round(midx - 1.5*xsize) > 0
      maxx = round(midx + 1.5*xsize) < tmpsize[1]-1
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
   
   
   minyc = round(centroidval - 1.25*ysize) > 0
   maxyc = round(centroidval + 1.25*ysize) < tmpsize[1]-1

   ;; check for case when too few pixels are found.  Could be due
   ;; to centroid being displaced.  just use midpoint

   if( (maxyc-minyc) le 5L) then begin
      midy = (maxy+miny)/2
      yplussize = maxy-midy
      yminussize = midy-miny
      ysize = max(yminussize, yplussize)
      miny = round(midy - 1.5*ysize) > 0
      maxy = round(midy + 1.5*ysize) < tmpsize[1]-1
   endif else begin
      miny = minyc
      maxy = maxyc
   endelse


   ; for PSF, cut it off where flux is <1000x peak
   tmpsize = size(psf)
   centroidval = 0.5*tmpsize[1]
   wkeep = where(psf ge 0.001*max(psf),nkeep)
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
        
   ;; write files
   nfile=i/per+1
   outgalfile = 'real_galaxy_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'
   outpsffile = 'real_galaxy_PSF_images'+label+'_n'+string(nfile,format='(I0)')+'.fits'


   print,'Writing to files!'
   mwrfits,img[minx:maxx,miny:maxy],outgalfile
   mwrfits,psf[minxp:maxxp,minyp:maxyp],outpsffile
   

       
   
   openw,1,varfile,/APPEND
   printf,1,mean,sig*sig
   close,1
                
    endelse
endfor

END
