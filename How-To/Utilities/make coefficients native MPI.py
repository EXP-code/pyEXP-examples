import os
import time                     # Used for timing the coefficient construction
import pyEXP
from mpi4py import MPI

from os.path import exists

#
# This script makes an HDF5 coefficient set from some EXP files in
# parallel using MPI.  This can easily be adapted for whatever
# snapshots you have and could be run on a cluster
#

# Usage:
#
# Run command is: "mpirun -np N python3 make_coefficients_MPI.py"
# where N is the number of processes.  For Slurm allocations, you can
# leave off "-np N" as usual.
#

if __name__ == "__main__":

    # Parameters
    #
    h5file  = 'Run1b_test_coefs'
    beg_seq = 0
    end_seq = 200
    nskip   = 1

    # Get basic information about the MPI communicator
    #
    world_comm = MPI.COMM_WORLD
    world_size = world_comm.Get_size()
    my_rank    = world_comm.Get_rank()

    # NB: we are not using standard MPI commands on the Python side,
    # but invoking mpi4py initializes MPI so it can be used correctly
    # by the C++ routines

    # Just for info and to convince yourself check that MPI is working
    #
    print("World size is {} and my rank is {}".format(world_size, my_rank))

    # Now switch the working directory where my Gadget simulation
    # lives.  Change this to your working directory.
    #
    os.chdir('/nas/astro-th/weinberg/Nbody/SimpleHalo')

    # Make the spherical basis config.
    #
    halo_config = """
id          : sphereSL
parameters  :
  numr      : 4000
  rmin      : 0.0001
  rmax      : 1.95
  Lmax      : 6
  nmax      : 20
  scale     : 0.0667
  modelname : SLGridSph.model
"""

    # Construct the basis instance
    #
    halo_basis = pyEXP.basis.Basis.factory(halo_config)

    # Make the file list for the snapshot sequence assumed have the
    # format: snapshot_XXXX.hdf5.  Change as necessary.
    #
    file_list = []
    for i in range(beg_seq, end_seq):
        file_list.append('OUT.run1b.{:05d}'.format(i))

    # Construct batches of files the particle reader.  One could use the
    # parseStringList to create batches from a vector/list of files.  NB:
    # a std::vector in C++ becomes a Python.list and vice versa
    #
    batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')

    # This will contain the coefficient container, need to start will a
    # null instance to trigger construction
    #
    coefs = None

    # One could replace the above with a read of an existent HDF5
    # coeffficient file and update various stanzas.  For this reason,
    # I am keeping a separate center time and center position list
    # below . . .

    for group in batches:
        okay = True
        for f in group:
            if not exists(f): okay = False
        
        # Skip a file that does not exist in the sequence without
        # going belly up.  I notice that Gadget-2 seems to not be
        # exactly sequential by 1 through a restart?
        #
        if not okay: continue

        if my_rank==0: print("file group is", group)

        # Make the reader for the desired type.  One could probably try to
        # do this by inspection but that's another project.
        #
        reader = pyEXP.read.ParticleReader.createReader('', group, 0, False);

        # Print the type list
        #
        if my_rank==0: print('The component names are:', reader.GetTypes())

        compname = 'dark'
        reader.SelectType(compname)
        if my_rank==0: print('We selected:', compname)

        startTime = time.time()

        # Now compute the coefficients with the default center
        #
        startTime = time.time()
        coef = halo_basis.createFromReader(reader)
        if my_rank==0:
            print('Created coefficients at Time {:5.3f} for {} particles '
                  'in {:4.2f} seconds'.
                  format(reader.CurrentTime(), reader.CurrentNumber(),
                         time.time() - startTime))

        # We need this stupid idiom here because None is not mapping to a
        # null pointer.  There is probably a way to do this.  Suggestions
        # anyone?
        #                          This is optional---+
        #                                             |
        if coefs is None:           #                 v
            coefs = pyEXP.coefs.Coefs.makecoefs(coef, compname)

        # Add the coefficients to the container
        #
        coefs.add(coef)

        if my_rank==0:
            print('Added coef to container')
            print('-'*60)

    if my_rank==0:
        print('\nCompleted the file group list\n')
        print('The coefficient time list is', coefs.Times())

        # You can call the file something convenient.  The suffix 'h5'
        # will be appended. You only want the root process to write
        # the file.
        #
        if exists(h5file):
            coefs.ExtendH5Coefs(h5file) # Update an existing HDF5
            print('Saved the coefficients to an existing HDF5 file')
        else:
            coefs.WriteH5Coefs(h5file) # Create a new HDF5
            print('Saved the coefficients to a new HDF5 file')

