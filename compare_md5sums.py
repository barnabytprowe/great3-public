#!/usr/bin/env python

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
        raise ValueError(
            "FAIL: file "+str(file1)+" has "+str(len(data1))+" records but file "+str(file2)+
            " has "+str(len(data2))+" records: FAIL")
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

    import argparse

    description = \
    """Compare two md5sums files output by get_md5sums.py, checking consistency.

    Reads in two full md5sum output files each with entries in the format

        <MD5 checksum> <filename>

    and compares the MD5 column for consistency (filenames, which may be different on different
    systems, are ignored for the purpose of this comparison).
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "file1", type=str,
        help="First ASCII file containing all md5sums calculated for one build of GREAT3")
    parser.add_argument(
        "file2", type=str,
        help="Second ASCII  file containing all md5sums calculated for another build of GREAT3")
    args = parser.parse_args()
    # Then compare
    compare_md5sum_files(args.file1, args.file2) 
