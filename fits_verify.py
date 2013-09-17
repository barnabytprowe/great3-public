#!/usr/bin/env python

import os

def fits_verify_file(filename, fitsverify_exec="fitsverify"):
    """Run fitsverify on the given `filename`.

    @param fitsverify_exec Executable file for fitsverify (default=`fitsverify`).

    @return True if fits file is verified, false if there was an error.
    """
    import subprocess
    retcode = subprocess.call([fitsverify_exec, "-q", "-e", filename])
    if retcode == 0:
        ret = True
    else:
        ret = False
    return ret

if __name__ is "__main__":

