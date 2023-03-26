import os
import time                     # Used for timing the coefficient construction
import pyEXP
from mpi4py import MPI

from os.path import exists

#
# This script makes an HDF5 coefficient set from some phase-space
# files in parallel using MPI.  This can easily be adapted for
# whatever snapshots you have and could be run on a cluster
#

# Usage:
#
# Run command is: "mpirun -np N python3 get_center_MPI.py"
# where N is the number of processes.  For Slurm allocations, you can
# leave off "-np N" as usual.
#

if __name__ == "__main__":

    # Parameters
    #
    ctrfile = 'new.centers'
    beg_seq = 0
    end_seq = 1000
    nskip   = 20

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

    # Make the file list for the snapshot sequence assumed have the
    # format: snapshot_XXXX.hdf5.  Change as necessary.
    #
    file_list = []
    for i in range(beg_seq, end_seq):
        file_list.append('snapshot_{:04d}.hdf5'.format(i))

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

    centime = []
    centers = []

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
        reader = pyEXP.read.ParticleReader.createReader('PSPspl', group, 0, False);

        # Print the type list
        #
        if my_rank==0: print('The component names are:', reader.GetTypes())

        compname = 'dark'
        reader.SelectType(compname)
        if my_rank==0: print('We selected:', compname)

        # This computes an expansion center from a mean density
        # weighted position.  You could compute and cache the center
        # array . . . or supply it in a different way
        #
        startTime = time.time()
        center = pyEXP.util.getDensityCenter(reader, stride=nskip, Nsort=1000)
        #                                            ^
        #                                            |
        # Choose every nskip particle for sample ----+
        # This is 10^6 samples for these snaps and nskip=20
        #
        if my_rank==0:
            print('Created center in {:4.2f} seconds'.
                  format(time.time() - startTime))
            print('Center is:', center)
            centime.append(reader.CurrentTime())
            centers.append(center)

    if my_rank==0:
        print('\nCompleted the file group list\n')

        # Save the center positions
        #
        with open(ctrfile, 'a') as f:
            for i in range(len(centime)):
                line = '{:13.6e} {:13.6e} {:13.6e} {:13.6e}\n'.format(centime[i], centers[i][0], centers[i][1], centers[i][2])
                f.write(line)
                
