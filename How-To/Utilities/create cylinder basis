#!/usr/bin/env python3
# coding: utf-8

from mpi4py import MPI
import numpy as np
import argparse
import pyEXP

# Helper for nice parallel printing
#
def pprint(str="", end="\n", comm=MPI.COMM_WORLD):
    """Print for MPI parallel programs: Only rank 0 prints *str*."""
    if comm.rank == 0:
        print(str+end, end='')
        
# A default and example cylindrical basis config
#
disk_config = """
id           : cylinder
parameters   :
  acyl       : 1.0
  hcyl       : 0.1
  lmaxfid    : 64
  nmaxfid    : 64
  mmax       : 10
  nmax       : 32
  ncylnx     : 256
  ncylny     : 128
  ncylr      : 2000
  rnum       : 200
  pnum       : 1
  tnum       : 80
  rcylmin    : 0.001
  rcylmax    : 20
  ashift     : 0
  logr       : true
  cachename  : eof.cache.file
"""

def main():
    global disk_config

    # Initialize MPI
    #
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    args = None

    # Only root process parses and sends arguments to all processes
    #
    if rank == 0:

        parser = argparse.ArgumentParser(
            prog="cylcache",
            description="Construction a cylindrical basis specified by a YAML config file",
            epilog="Running with no arguments will use the default config provided as a template. Typical usage: mpirun -np 8 cylcache -c my_config.yaml"
        )
        parser.add_argument("-c", "--config", dest="config_file", help="The YAML basis config used construct a basis")
        parser.add_argument("-t", "--template", dest="template_file", help="Write a sample default YAML basis config and exit")

        # Parse in input arguments
        #
        try:
            args = parser.parse_args()
        except:
            args = None
    
        # Pass the namespace to all processes
        #
        args = comm.bcast(args, root=0)
    
        # Help?
        #
        if args is None:
            MPI.Finalize()
            exit(0)
    else:
        # Receive the namespace from the root process
        #
        args = comm.bcast(None, root=0)

        # Help?
        #
        if args is None:
            MPI.Finalize()
            exit(0)

    # User asked for a template YAML file
    #
    if (args.template_file):
        if rank == 0:
            file = open(args.template_file, 'w')
            file.write(disk_config)
        MPI.Finalize()
        exit(0)

    # User passed a config file
    #
    if (args.config_file):
        file = open(args.config_file)
        disk_config = file.read()

    # Just for info and to convince yourself check that MPI is working
    #
    pprint("============================================================================")

    # Begin calculation and start stopwatch
    #
    t_start = MPI.Wtime()

    # Construct the basis instance
    #
    disk_basis = pyEXP.basis.Basis.factory(disk_config)

    comm.Barrier()

    pprint("---- Calculation finished")

    # Stop stopwatch
    #
    comm.Barrier()
    t_diff = MPI.Wtime() - t_start

    pprint("---- Computed basis in {}".format(t_diff))
    pprint("============================================================================")
    MPI.Finalize()

if __name__ == "__main__":
    main()
