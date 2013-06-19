PRO run_many,label,start,n,per,seed

;; Make an initial catalog with filenames,properties, etc.  

makecatalog_many,label,start,n,per,seed

;; Create postage stamps and calculate noise properties

makestamps_many,label,per,'real_galaxy_var.txt'

;; create new catalog with noise properties included
;; mostly repeat of previous step, but doesn't take very long
makecatalog_many_var,label,start,n,per,'real_galaxy_var.txt',seed

end

