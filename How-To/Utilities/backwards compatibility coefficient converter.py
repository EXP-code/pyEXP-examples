#!/usr/bin/env python3
# coding: utf-8

# Convert from the original EXP SphericalBasis coefficient
# normalization to the new EXP/pyEXP convention.  The new convention
# is pure orthonormal.

import os, sys, getopt
import pyEXP

def help(phrase: str) -> None:
    """Print some usage info"""
    print(phrase)

def index(k):
    """Return the values of l, m associated with a spherical index"""
    n = 0
    l = 0
    while (k>=n):
        m = k - n
        if m<=l: return (l, m)
        l += 1
        n += l

def main(prog, argv) -> int:
    """Convert from old EXP SphericalBasis coefficient normalization to
    new EXP/pyEXP orthonormal convention

    """

    ifile = ''
    ofile = ''

    phrase = prog + ': [-h] -i|--old=file -o|--new=file';

    try:
        opts, args = getopt.getopt(argv,"hi:o:",["old=","new="])
    except getopt.GetoptError:
        help(phrase)
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            help(phrase)
            sys.exit()
        elif opt in ("-i", "--old"):
            ifile = arg
        elif opt in ("-o", "--new"):
            ofile = arg

    # Check for input file
    if not os.path.exists(ifile):
        message = prog + ': can not open <' + ifile + '>'
        print(message)
        exit(-1)

    # Make sure output file name is available
    if os.path.exists(ofile):
        message = prog + ': <' + ofile + '> exists!'
        print(message)
        exit(-1)

    # Read the original coefficient db
    coefs = pyEXP.coefs.Coefs.factory(ifile)

    # m=0 factor
    fac0 = -1.0

    # m>0 factor
    fac1 = -2.0**0.5

    for t in coefs.Times():
        mat = coefs(t)
        for k in range(mat.shape[0]):
            l, m = index(k)
            if m==0: mat[k,:] *= fac0
            else:    mat[k,:] *= fac1
            coefs.setMatrix(t, mat)

    # Save the converted coefficient db
    coefs.WriteH5Coefs(ofile)

if __name__ == "__main__":
   main(sys.argv[0], sys.argv[1:])
