#!/usr/bin/env python

"""Script for calculating the GREAT08 Q metric (Q08) from a constant shear submission.

The constant shear submission format takes the form of a list of average g1, g2 values, calculated
on an image by image basis.  The input to this code is a catalogue with three columns as follows:

# image_ID g1 g2

By comparing this information to tables with the true input shear, we calculate Q08.
"""

import os
import sys
import numpy as np

