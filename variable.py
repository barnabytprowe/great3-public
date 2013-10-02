#!/usr/bin/env python
import os
import sys
import numpy as np
import constants
sys.path.append("..")
sys.path.append(os.path.join("..", "server", "great3"))
import great3sims.mapper
import evaluate
import test_evaluate

def get_field(subfield_index, nsubfields_per_field=evaluate.NSUBFIELDS_PER_FIELD):
    return subfield_index / nsubfields_per_field

def get_subfield_within_field(subfield_index, nsubfields_per_field=evaluate.NSUBFIELDS_PER_FIELD):
    return subfield_index % nsubfields_per_field

def make_fits_cats(idarray, g1array, g2array, dir="./cats", prefix="vartest"):
    """Make a whole series of FITS catalogues of the sort expected by Melanie's presubmission.py
    script.
    
    Writes a whole series of catalogues named <prefix>-000.fits, <prefix>-001.fits etc. to the
    folder specified as `dir`.
    """
    import pyfits
    if g1array.shape != g2array.shape:
        raise ValueError("Input g1true and g2true must be same shape")
    if g1array.shape != (
        evaluate.NGALS_PER_SUBFIELD, evaluate.NSUBFIELDS_PER_FIELD, evaluate.NFIELDS):
        raise ValueError(
            "Input truth shear arrays do not match shape expected for GREAT3 simulations")
    # Loop over subfields making the catalogues
    for subfield_index in range(evaluate.NSUBFIELDS):

       field_index = get_field(subfield_index)
       subfield_within_field_index = get_subfield_within_field(subfield_index)
       g1 = g1array[:, subfield_within_field_index, field_index]
       g2 = g2array[:, subfield_within_field_index, field_index]
       identifier = idarray[:, subfield_within_field_index, field_index]
       col0 = pyfits.Column(name='ID', format='9A', array=identifier)    
       col1 = pyfits.Column(name='g1', format='E', array=g1)
       col2 = pyfits.Column(name='g2', format='E', array=g2)
       col3 = pyfits.Column(name='w', format='E', array=np.ones_like(g1))
       cols = pyfits.ColDefs([col0, col1, col2, col3])
       tbhdu = pyfits.new_table(cols)
       prhdu = pyfits.PrimaryHDU()
       thdulist = pyfits.HDUList([prhdu, tbhdu])
       outfile = os.path.join(dir, prefix+("-%03d" % subfield_index)+".fits")
       print "Saving FITS catalogue to "+outfile
       thdulist.writeto(outfile, clobber=True)
    return

if __name__ == "__main__":

    idt, xt, yt, g1t, g2t = test_evaluate.get_variable_gtrue("control", "space")
    ret = make_fits_cats(idt, g1t, g2t)
