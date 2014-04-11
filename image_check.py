#!/usr/bin/env python

import os
import sys
import glob
import pyfits
import numpy
# These two lines allow you to run this script on a terminal with no Display defined, 
# eg through a PBS queue
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import constants
sys.path.append("..")
try:
    import great3sims.mapper
    from great3sims.constants import n_subfields, n_deep_subfields
except:
    n_subfields = 220
    n_deep_subfields = 20

nproc = 8
threshold = 3 # Report deviations >= this many sigma for the means, medians, etc
n_histogram_bins = 40
n_curves_per_histogram_plot = 20 
histogram_root_filename = 'pixel_histograms/histogram'
# If you change histogram_root_filename to point to a different directory, you can delete
# the next two lines.
if not os.path.exists('pixel_histograms'):
    os.makedirs('pixel_histograms')

if nproc>1:
    # See if we can import the multiprocessing module
    try:
        import multiprocessing
    except:
        print "Cannot import multiprocessing; setting nproc = 1 and running serially."
        nproc = 1
# Luminosity-matched colors
colors = ['#7c7c7c', '#db00db', '#6767ff', '#008989', 
          '#009000', '#808000', '#d35400', '#f80000']
len_colors = len(colors)
n_main_fields = n_subfields - n_deep_subfields

def draw_histogram(x,cnum=0):
    """
    A function to draw data from a numpy.histogram() call as a step-style histogram.
    This function doesn't do any matplotlib things EXCEPT the plot call--writing must
    be done in the calling function.
    """
    # Stop numpy printing errors when it tries to take log(0)
    numpy.seterr(all='ignore') 
    # numpy.histogram returns a tuple (frequencies, bin_edges)
    frequencies = numpy.log(x[0])
    bins=x[1]
    line_x = bins[:-1]
    # Offset the curves so they're easier to see
    line_y = 0.1*(cnum%n_curves_per_histogram_plot)+frequencies
    plt.step(line_x,line_y,color=colors[cnum%len_colors],where='post')
    # Return to previous numpy floating-point error behavior, ie to warn but not stop
    numpy.seterr(all='warn')
    return

# So we can pool.map without worrying about kwargs--kind of a hack, but easy
expected_size = 4800

def image_good_sizeonly(filename):
    """Run an image check on the given `filename`.

    @return a tuple of 8 values: 3 logical values for correct size, no NaNs, and no zero-valued
    pixels; then zeros to fill out the array so the return signature is the same as image_good_full.
    """
    image_file = pyfits.open(filename)
    image = image_file[0].data
    size = image.shape
    # Do basic FITS-file pixel level checks
    if size[0]==size[1] and size[0]==expected_size:
        good_size = True
    else:
        good_size = False
    if numpy.any(numpy.isnan(image)):
        good_pixels = False
    else:
        good_pixels = True
    if numpy.any(image==0):
        no_zero_pixels = False
    else:
        no_zero_pixels = True
    image_file.close()
    return good_size, good_pixels, no_zero_pixels, 0, 0, 0, 0, 0
     

def image_good_full(filename):
    """Run an image check on the given `filename` and make pixel histograms.

    @return a tuple of 8 values: 3 logical values for correct size, no NaNs, and no zero-valued
    pixels; the mean, standard deviation, median, and median absolute deviation; and a Numpy array
    of histogram values.
    """
    image_file = pyfits.open(filename)
    image = image_file[0].data
    size = image.shape
    # Do basic FITS-file pixel level checks
    if size[0]==size[1] and size[0]==expected_size:
        good_size = True
    else:
        good_size = False
    if numpy.any(numpy.isnan(image)):
        good_pixels = False
    else:
        good_pixels = True
    if numpy.any(image==0):
        no_zero_pixels = False
    else:
        no_zero_pixels = True
    # Do a histogram.  So many of the values are near 0 that doing log(pixel value) is more
    # easily readable; since some values are <1, do log(1+pixel value).
    histogram = numpy.histogram(numpy.log(1.+image),bins=n_histogram_bins)
    # Compute some summary statistics.
    mean = numpy.mean(image)
    stddev = numpy.std(image)
    median = numpy.median(image)
    mad = numpy.median(numpy.abs(image-median))
    image_file.close()
    return good_size, good_pixels, no_zero_pixels, mean, stddev, median, mad, histogram

def check_all(root_dir, experiments=constants.experiments, obs_types=constants.obs_types,
              shear_types=constants.shear_types, full_stats = True, full_stats_filename = None):
    """Check all image files are good in branches for the experiments, obs_types and shear_types
    specified,  within the specified root directory.

    Assumes all image filenames begin with "image" and end in suffix .fits.
    
    full_stats = True computes mean, stddev, median, and MAD.  If you want that information written
    to a file, give the kwarg full_stats_filename; if full_stats = False this kwarg is ignored.
    
    @return good  If all checks pass, returns a list of all the filenames that were found by
                  this check and passed the correct-size checks.  If any checks failed, prints all
                  failed filenames and raises a TypeError exception.
    """
    if full_stats:
        image_function = image_good_full
        if full_stats_filename:
            fsf = open(full_stats_filename,'w')
    else:
        image_function = image_good_sizeonly
    # Set storage lists
    good_sizes = []
    bad_sizes = []
    good_NaNs = []
    bad_NaNs = []
    good_zeros = []
    bad_zeros = []
    found_exps = []
    found_obs = []
    found_shears = []
    for experiment in experiments:
        
        for obs_type in obs_types:
            
            for shear_type in shear_types:

                # Check for all files matching image*.fits.
                # The try-except structure is here because msimet was using rmandelb's files,
                # rather than untarring everything, and didn't have permission to make new 
                # directories, which caused the mapper to fail when branches didn't exist.
                try:
                    mapper = great3sims.mapper.Mapper(root_dir, experiment, obs_type, shear_type)
                    fitsfiles = glob.glob(os.path.join(mapper.full_dir, "image*.fits"))
                except:
                    fitsfiles = glob.glob(os.path.join(root_dir, experiment, obs_type, shear_type, "image*.fits"))
                fitsfiles.sort() # So they're in numerical order
                # If any fits files found, check
                if len(fitsfiles) > 0:
                    # Reset for every branch since these may be different branch to branch.
                    means = [] 
                    stddevs = []
                    medians = []
                    mads = []
                    histograms = []
                    if experiment not in found_exps: found_exps.append(experiment)
                    if obs_type not in found_obs: found_obs.append(obs_type)
                    if shear_type not in found_shears: found_shears.append(shear_type)
                    # Expected image size for the different branches
                    global expected_size
                    if obs_type=="space" and not (experiment=="multiepoch" or experiment=="full"):
                        expected_size=9600
                    else:
                        expected_size=4800
                    if nproc>1:
                        pool = multiprocessing.Pool(nproc)  
                        results = pool.map(image_function,fitsfiles)
                    else:
                        results = [image_function(fitsfile) for fitsfile in fitsfiles]
                    # Turn the list of results into usefully-named quantities by
                    # transposing the list of tuples
                    (query_sizes, query_NaNs, query_zeros, 
                        means, stddevs, medians, mads, histograms) = zip(*results)
                    for i,good_size in enumerate(query_sizes):
                        if good_size:
                            good_sizes.append(fitsfiles[i])
                        else:
                            bad_sizes.append(fitsfiles[i])
                    for i,good_NaN in enumerate(query_NaNs):
                        if good_NaN:
                            good_NaNs.append(fitsfiles[i])
                        else:
                            bad_NaNs.append(fitsfiles[i])
                    for i,no_zero_pixels in enumerate(query_zeros):
                        if no_zero_pixels:
                            good_zeros.append(fitsfiles[i])
                        else:
                            bad_zeros.append(fitsfiles[i])
                            
                    if full_stats:
                        # Now check for outliers in the mean, stddev, median, or MAD
                        # and print a warning (but keep analyzing) for any found
                        qlist = [(means,"mean"),(stddevs,"stddev"),(medians,"median"),(mads,"MAD")]
                        for long_quantity, qname in qlist:
                            # Cut to just the main files--the deep subfields are often different,
                            # and too few to make good outlier calculations; make sure to examine
                            # the histograms for any weirdness.
                            quantity = long_quantity[:n_main_fields]
                            mean = numpy.mean(quantity)
                            stddev = numpy.std(quantity)
                            # [0] since numpy.where returns a tuple with an element for each
                            # dimension, but we only have (and care about) one.
                            outliers = numpy.where(
                                numpy.abs((quantity-mean)/stddev)>threshold)[0]
                            if len(outliers)>0:
                                for outlier in outliers:
                                    print "WARNING: Outlier (>"+str(threshold)+' sigma) for',
                                    print 'quantity', qname,':', fitsfiles[outlier],
                                    print 'value =', quantity[outlier], 'for mean',mean,
                                    print 'and standard deviation', stddev
                        histogram_filename_base = (histogram_root_filename+'.'+experiment+'.'+
                                                   obs_type+'.'+shear_type)
                        # Draw some png files with a huge number of histograms...
                        nhists=0
                        for i,hist in enumerate(histograms):
                            draw_histogram(hist,cnum=i)
                            if (i+1)%n_curves_per_histogram_plot==0:
                                plt.title('Pixel histogram for'+experiment+'-'+obs_type+
                                          '-'+shear_type)
                                plt.xlabel('log(1+pixel value)')
                                plt.ylabel('log(Number of pixels)')
                                plt.savefig(histogram_filename_base+'-'+str(nhists)+'.png')
                                plt.clf()
                                nhists+=1
                        # In case len(histograms) doesn't divide evenly into n_curves_per
                        if (i+1)%n_curves_per_histogram_plot!=0:
                            plt.title('Pixel histogram for '+experiment+'-'+obs_type+'-'+shear_type)
                            plt.xlabel('log(1+pixel value)')
                            plt.ylabel('log(Number of pixels)')
                            plt.savefig(histogram_filename_base+'-'+str(nhists)+'.png')
                            plt.clf()
                        if full_stats_filename:
                            fsf.write('# '+experiment+'-'+obs_type+'-'+shear_type+'\n')
                            fsf.write('# mean stddev median MAD\n')
                            for i in range(len(fitsfiles)):
                                fsf.write(fitsfiles[i]+' '+str(means[i])+' '+str(stddevs[i])+' '+
                                          str(medians[i])+' '+str(mads[i])+'\n')

    if full_stats and full_stats_filename:
        fsf.close()

    print "Ran image analysis on all image files in "+root_dir
    print "Found images for experiments "+str(found_exps)
    print "Found images for obs_types "+str(found_obs)
    print "Found images for shear_types "+str(found_shears)
    if len(bad_sizes) > 0:
        print "The following images had incorrect sizes:\n"
        for filename in bad_sizes:
            print filename
    if len(bad_NaNs) > 0:    
        print "The following images had NaNs:\n"
        for filename in bad_NaNs:
            print filename
    if len(bad_zeros) > 0:    
        print "The following images had pixels == 0:\n"
        for filename in bad_zeros:
            print filename
    if bad_sizes or bad_NaNs or bad_zeros:
        raise TypeError('Some failures found.  See stdout for more details.')
    else:
        print "All images analyzed successfully!"
        retlist = good_sizes
    return retlist

if __name__ == "__main__":

    good_public = check_all(constants.public_dir, full_stats_filename="stats.dat")
