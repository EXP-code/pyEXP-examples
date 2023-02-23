"""
an example script for making some MW orbits from the Sgr simulation of Laporte+ 2018

this version of the script assumes you have already computed the bases elsewhere.

"""

# the first part of this script follows (a non parallel version of) EXP-code/pyEXP-examples/external-simulations/make_coefficients_MPI.py

# you will want to create arrays of coefficients for each of the halo, bulge, and disc components.

import os
import pyEXP
import matplotlib.pyplot as plt

# most functions have documentation, which may be accessed like this:
print(pyEXP.read.ParticleReader.createReader.__doc__)

# change this the directory where you have downloaded the snapshots AND cachefiles.
basedir = '/Volumes/External1/ChervinSgrSim/'
os.chdir(basedir)

# this can be a list of snapshots, or you can loop through one at a time. For now, I'm just practicing on the first snapshot.
group = ['snap_001']

# read in the actual data
reader = pyEXP.read.ParticleReader.createReader('GadgetNative', group, 0, True); # this will take a tens of seconds.
print('The component names are:', reader.GetTypes())
# The component names are: ['Bndry'=Sgr stars (Hernquist), 'Bulge'=MW bulge, 'Disk'=MW disk, 'Halo'=MW halo, 'Stars'=DM halo of Sgr]
# We currently only care about the MW halo, bulge, and disk. We'll expand them below.

# the first step is to find the center of the simulation. Let's assume that is the halo centre for now. We will use this centre for all components.
compname = 'Halo'
reader.SelectType(compname) # this will take a few seconds.
nskip = 50 # this needs to be a relatively small number to get a good centre. We might be able to get away with a little larger (if gadget sims are in kpc)
center = pyEXP.util.getDensityCenter(reader, nskip) # this will take tens of seconds.

# we do need to record the center at each timestep for later translation, so let's open a file and print
# only open the file and print the header at the first timestep!
centerfile = open('simulationcenters.txt','w')
print('{:13s} {:13s} {:13s} {:13s}'.format('time','xcentre','ycentre','zcentre'),file=centerfile)

# print to the file
print('{:13.6e} {:13.6e} {:13.6e} {:13.6e}'.format(reader.CurrentTime(), center[0], center[1], center[2]),file=centerfile)

# now, let's generate a basis for the halo. These are the settings, no need to touch these for now.
halo_config = """
id         : sphereSL
parameters :
  numr     : 4000
  rmin     : 0.06
  rmax     : 1000.0
  Lmax     : 6
  nmax     : 20
  rs       : 0.176
  modelname : SLGrid.empirical.mw.isolated.ieee
  cachename : SLGridSph.cache.mw
"""

halo_basis = pyEXP.basis.Basis.factory(halo_config)

# make the coefficients
compname = 'Halo'
reader.SelectType(compname) # this will take a few seconds.
halo_coef = halo_basis.createFromReader(reader, center) # this will take a few seconds

# make a makecoefs instance
# only do this at the first step: afterwords just add (see below).
halo_coefs = pyEXP.coefs.Coefs.makecoefs(halo_coef, compname)

# add the coefficients to the makecoefs instance
halo_coefs.add(halo_coef)

# write the first step (afterwards use ExtendH5Coefs)
halo_coefs.WriteH5Coefs('outcoef.mwhalo.L18')

# normally, we would go back and re-read a new file, make new coefficients, etc, but in this case we are just practicing with a single snapshot.
# so we'll adjust the time manually, add the coefficients again, and extend the written coefficient file
halo_coef.time = 1000.0
halo_coefs.add(halo_coef)
halo_coefs.ExtendH5Coefs('outcoef.mwhalo.L18')

# same procedure as above, but now for the bulge
# now, let's generate a basis for the bulge
bulge_config = """
id         : sphereSL
parameters :
  numr     : 4000
  rmin     : 0.01
  rmax     : 1000
  Lmax     : 0
  nmax     : 10
  rs       : 1
  modelname : SLGrid.empirical.bulge.isolated.ieee
  cachename : SLGridSph.cache.bulge
"""
bulge_basis = pyEXP.basis.Basis.factory(bulge_config)

# generate the coefficients
compname = 'Bulge'
reader.SelectType(compname) # this will take a few seconds.
bulge_coef = bulge_basis.createFromReader(reader, center) # this will take a few seconds

# make a makecoefs instance
# only do this at the first step: afterwords just add (see below).
bulge_coefs = pyEXP.coefs.Coefs.makecoefs(bulge_coef, compname)
bulge_coefs.add(bulge_coef)
bulge_coefs.WriteH5Coefs('outcoef.mwbulge.L18')

# do the false duplication for extra time
bulge_coef.time = 1000.0
bulge_coefs.add(bulge_coef)
bulge_coefs.ExtendH5Coefs('outcoef.mwbulge.L18')

# now, the disc. The configuration looks different, but the process is exactly the same.
disk_config = """
id         : cylinder
parameters :
  acyl       : 3.0
  hcyl       : 0.5
  mmax       : 2
  ncylorder  : 12
  ncylnx     : 128
  ncylny     : 64
  rnum       : 200
  pnum       : 80
  tnum       : 80
  lmax       : 48
  nmax       : 36
  logr       : true
  density    : true
  eof_file   : CylGrid.cache.mwdisk
"""
disk_basis = pyEXP.basis.Basis.factory(disk_config)

compname = 'Disk'
reader.SelectType(compname) # this will take a few seconds.
disk_coef = disk_basis.createFromReader(reader, center) # this will take a few seconds

# make a makecoefs instance
# only do this at the first step: afterwords just add (see below).
disk_coefs = pyEXP.coefs.Coefs.makecoefs(disk_coef, compname)
disk_coefs.add(disk_coef)
disk_coefs.WriteH5Coefs('outcoef.mwdisk.L18')

# do the false duplication for extra time
disk_coef.time = 1000.0
disk_coefs.add(disk_coef)
disk_coefs.ExtendH5Coefs('outcoef.mwdisk.L18')



# now, we can set up the whole instance to integrate orbits
# bring in all coeficients and bases fresh
cfiles = {'mwdisc': 'outcoef.mwdisk.L18.h5', 'mwhalo': 'outcoef.mwhalo.L18.h5', 'mwbulge': 'outcoef.mwbulge.L18.h5'}
bases  = {'mwdisc': disk_basis, 'mwhalo': halo_basis, 'mwbulge': bulge_basis}
coefs  = dict()
for v in cfiles:
    coefs[v] = pyEXP.coefs.Coefs.factory(cfiles[v])

# create the list of model tuples
model = [[bases['mwdisc'], coefs['mwdisc']], [bases['mwhalo'], coefs['mwhalo']], [bases['mwbulge'], coefs['mwbulge']]]

# define the scaling for G=1
astronomicalG = 0.0043009125 # pc/Msun/km/s/km/s

# define a solar-like orbit (x y z vx vy vz)
# don't forget that we have to scale velocities to account for G=1
phasespace = [[8.0,0.,0.,astronomicalG*0.,astronomicalG*220.,astronomicalG*0.]]

# we need to add the centre to the positions. so start from the sun, but add the offset computed above.
phasespace = [[8.0+center[0],0.+center[1],0.+center[2],astronomicalG*0.,astronomicalG*220.,astronomicalG*0.]]

# integrate the orbit: arguments are tinitial,tfinal,tstep,phasespace,model
times, orbits = pyEXP.basis.IntegrateOrbits(4., 1000.0, 0.1, phasespace, model, pyEXP.basis.AllTimeAccel())

# visualise the orbit
plt.plot(orbits[0][0],orbits[0][1],color='black')

# note that the time is also scaled because of G=1. From EXP, the time units can be multiplied by astronomicalG to get (approximately) Gyr.

# if we still have the centerfile open, close it.
centerfile.close()
