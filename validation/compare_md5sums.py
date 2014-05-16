#!/usr/bin/env python

# Copyright (c) 2014, the GREAT3 executive committee (http://www.great3challenge.info/?q=contacts)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions
# and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of
# conditions and the following disclaimer in the documentation and/or other materials provided with
# the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to
# endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""Compare two md5sums files output by get_md5sums.py, checking consistency.

Reads in two full md5sum output files each with entries in the format

    <MD5 checksum> <filename>

and compares the MD5 column for consistency. The filenames column, which may be different on
different systems, is ignored for the purpose of this comparison.
"""

import os

def compare_md5sum_files(file1, file2):
    """Check two get_md5sum.py output files for consistency.
    """
    import numpy as np
    # Load the data
    with open(file1, "rb") as funit1:
        data1 = funit1.readlines()
    with open(file2, "rb") as funit2:
        data2 = funit2.readlines()
    # Check record length, totally basic!
    if len(data1) != len(data2):
        print (
            "FAIL: file "+str(file1)+" has "+str(len(data1))+" records but file "+str(file2)+
            " has "+str(len(data2))+" records")
        exit(1)
    # Then get the md5sum bits separated from the filename bits
    md1list = [ditem.split()[0] for ditem in data1]
    md2list = [ditem.split()[0] for ditem in data2]
    file1list = [(ditem.split()[1]).rstrip("\n") for ditem in data1]
    file2list = [(ditem.split()[1]).rstrip("\n") for ditem in data2]
    # Then check for inequality, if found look where things went wrong
    if md1list != md2list:
        import sys
        for md1, md2, filename1, filename2 in zip(md1list, md2list, file1list, file2list):

            if md1 != md2:
                print "FAIL: md5sums differ!"
                print "Mismatch:"
                print filename1+" ("+md1+")"
                print filename2+" ("+md2+")"
                print

        sys.exit(1)
    else:
        print "All md5sums in "+file1+" and "+file2+" check out"
    return

if __name__ == "__main__":

    from sys import argv
    if len(argv) != 3:
        print "usage: compare_md5sums.py file1 file2"
        exit(1)
    # Then compare
    compare_md5sum_files(argv[1], argv[2]) 
