# pyEXP-examples/How-To/Utilities

Python scripts for various pipeline processing tasks

## Organization

| Script name    | Task        |
| ---            | ---         |
| viewSphBasis   | Read the SphericalBasis cache file and plot the basis functions using pyplot |
| viewCylBasis   | Read the Cylinder cache file and plot the basis functions using pyplot |
| orthoCheckHalo.py | Generate and test the biorthogonality of a spherical basis using pyEXP |
| orthoCheckDisk.py | Read cache file for a Cylinder basis and test the biorthogonality using pyEXP |
| get_center_MPI.py | Reads a set of snapshots and estimates the density center for a particular component for later use in coefficient generation |
| ascii_to_coefficients.py | Read a phase-space body table, makes a set of coefficients, and adds them to an HDF coefficient file |
| make_cylinder_basis.py | Makes a biorthogonal basis using pyEXP.  This MPI version can be run on a cluster or multicore workstation... |
| coefficient converter | Convert from the original EXP SphericalBasis coefficient normalization to the new EXP/pyEXP convention.  This will only be relevant for old EXP coefficient sets |
