# pyEXP-examples/How-To/Utilities

Python scripts for various pipeline processing tasks

## Organization

| Script name    | Task        |
| ---            | ---         |
| backwards compatibility coefficient converter.py | Convert from the original EXP SphericalBasis coefficient normalization to the new EXP/pyEXP convention.  This will only be relevant for old EXP coefficient sets |
| check Cylinder basis orthogonality.py | Read cache file for a Cylinder basis and test the biorthogonality using pyEXP |
| check Spherical basis orthogonality.py | Generate and test the biorthogonality of a spherical basis using pyEXP |
| convert ascii table to coefficients.py | Read a phase-space body table, makes a set of coefficients, and adds them to an HDF coefficient file |
| create Cylinder basis (parallel).py | Makes a biorthogonal basis using pyEXP.  This MPI version can be run on a cluster or multicore workstation... |
| find center (parallel).py | Reads a set of snapshots and estimates the density center for a particular component for later use in coefficient generation |
| visualize Spherical basis.py   | Read the SphericalBasis cache file and plot the basis functions using pyplot |
| visualize gravitational power.py | Plot gravitational power of a coefficient set |
| visualize Cylinder basis.py   | Read the Cylinder cache file and plot the basis functions using pyplot |
