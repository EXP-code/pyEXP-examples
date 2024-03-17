#!/usr/bin/env python
# coding: utf-8

# ===============================================
# An example of creating coefficients using pyEXP
# ===============================================

# This example uses a prerun simulation whose output stored in the
# 'Tutorials/Data' directory. This can be recreated using 'EXP' and
# the example configuration and input files from the
# 'EXP-examples/Nbody' simulation.  But this notebook could be adapted
# for any simulation you like.
# 
# We begin by importing 'pyEXP' and friends and setting the working directory.

import os
import yaml
import pyEXP

# As described in the README, we assume that you have started Jupyter
# in the 'Tutorials' directory and have provided all of the necessary
# data in the 'Data' subdirectory.  We move to that directory as a
# first step:
#
os.chdir('../Data')


# First, we need to create the basis.  We'll only do the halo
# coefficients in this simple example.  The cylindrical coefficients
# would procede similarly.  See the 'viewing a basis' notebook for an
# example of creating the cylindrical basis.
#
# We start by importing the basis configuration.  The next bit of code
# loads a YAML stanza specifies all of the input parameters necessary
# to build the basis.
#
yaml_config = ""
with open('basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Alternatively, you could construct this on the fly in Python, e.g.
bconfig = """
---
id: sphereSL
parameters :
  numr: 2000
  rmin: 0.0001
  rmax: 1.95
  Lmax: 4
  nmax: 10
  rmapping: 0.0667
  modelname: SLGridSph.model
  cachename: sphereSL.cache
...
"""
print('-'*60)
print('Read from file')
print('-'*60)
print(yaml_config)
print('-'*60)
print('Constructed from string')
print('-'*60)
print(bconfig)
print('-'*60)

# Construct the basis instance
#
basis   = pyEXP.basis.Basis.factory(yaml_config)


# When an adaptive basis such as 'sphereSL' or 'Cylindrical' is
# constrructed, the orthogonality of the biorthogonal basis is
# automatically computed as a sanity check on the model file and input
# parameters.  The condition reported by 'orthoTest' is the absolute
# value of the inner product minus one or the inner product itself,
# for diagonal and off-diagonal matrix elements, respectively.  In
# this case, you can see that the worst inner product is ~0.0001 or on
# part in $10^4$.

# Now that we have a basis, we can use it to create coefficients from the particle snapshots.  'pyEXP' uses a 'ParticleReader' object for that.
# 
# The first step is to hand off the files that comprise a snapshot for
# every time slice.  The 'ParticleReader' provides a helper function
# for that.  There are two helper functions: 'parseFileList' and
# 'parseStringList'.  The first reads a list from a file and the
# second takes a list.  Otherwise they are the same.  The file names
# in the list are assumed to end with a snapshot index and an optional
# part index.  For example, if you have single files per snapshot, the
# list might look like: 'myrun.00000', 'myrun.00001', etc.  If you
# have multiple files per snapshot, they will look something like
# 'myrun.00000_0001', 'myrun.00000_0002', 'myrun.00001_0000',
# 'myrun.00001_0001', etc.
#
# One could use the parseStringList to create batches from a
# vector/list of files.  NB: a std::vector in C++ becomes a
# Python.list and vice versa
#
batches = pyEXP.read.ParticleReader.parseFileList('file.list', '')


# We now iterate the 'batches' created by the file parser to create
# the coefficients.  For each batch we create a new reader and pass
# the reader to the basis instance.  The 'basis.createFromReader'
# member creates and returns the coefficients.  The coefficients are
# added to a coefficient container called 'coefs'.  Note: on the first
# call 'coefs=None' so a new container is created on the first time
# through.

coefs   = None

for group in batches:

    print("file group is", group)

    # Make the reader for the desired type.  One could probably try to
    # do this by inspection but that's another project.
    #
    reader = pyEXP.read.ParticleReader.createReader('PSPout', group, 0, False);

    # Print the type list
    #
    print('The component names are:', reader.GetTypes())

    compname = 'dark halo'
    reader.SelectType(compname)
    print('Selected', compname)

    print('Call createFromReader at Time', reader.CurrentTime(), 'for', reader.CurrentNumber(), 'particles')

    coef = basis.createFromReader(reader)
    print("Created coef")

    # We need this idiom here because None is not mapping to a
    # null pointer in pybind11.  There is probably a way to do this.  
    # Suggestions anyone?
    #
    #                          This is optional---+
    #                                             |
    if coefs is None:           #                 v
        coefs = pyEXP.coefs.Coefs.makecoefs(coef, compname)
    else:
        coefs.add(coef)

    print('Added coef')
    print('-'*60)

print('\nCompleted the file group list\n')

print('The coefficient time list is', coefs.Times())


# The reader gives a verbose status update for each snapshot as it is
# read.  You can see the component names, simulation time, and
# particles numbers per snapshot as the files are read.  For large
# simulations, coefficient construction can take awhile.  Check out
# the following two scripts in the 'Utitilies' directory:
#
# 1. make coefficients MPI.py
#
# 2. make coefficients native MPI.py
#    
# The first one is an example for Gadget snapshots and the second is
# an example for EXP snapshots.
# 

# Now that we have our new coefficients, we can use the
# 'FieldGenerator' to view the BFE representation of the underlying
# fields.  Here is an example:

times = coefs.Times()[0:3]
pmin  = [-1.0, -1.0, 0.0]
pmax  = [ 1.0,  1.0, 0.0]
grid  = [  40,   40,   0]

fields = pyEXP.field.FieldGenerator(times, pmin, pmax, grid)

surfaces = fields.slices(basis, coefs)

print("We now have the following [time field] pairs")
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))

print("\nHere is the first one:")
for v in surfaces:
    for u in surfaces[v]:
        print('-'*40)
        print('----', u)
        print('-'*40)
        print(surfaces[v][u])
    break


# These could be make into images and so forth.  We'll do this in
# another example notebook.
# 
# At this point, it makes sense to save the coefficients that you have
# just created.  This is sone with the following call:


coefs.WriteH5Coefs('test_coefs')


# We now have a EXP HDF5 coefficient file called 'test_coefs'.  You
# can view the contents of 'test_coefs' directly using the 'h5dump'
# command supplied in your HDF5 installion:
#
# h5dump test_coefs | less
# 
# As an example of manipulating the newly made coefficients in pyEXP,
# let's try reading the newly created file into another coefficient
# container, 'coefs2'.  The container has a member function called
# 'CompareStanzas' which will check on the contents.  Let's do it.
#
coefs2 = pyEXP.coefs.Coefs.factory('test_coefs')
print("Type is", coefs2.getGeometry())

# Now compare with the original
#
coefs2.CompareStanzas(coefs)


# This member function will print differences.  No differenced should
# be printed, of course.
