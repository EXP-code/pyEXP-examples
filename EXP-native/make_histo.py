import os, sys
import pyEXP
import numpy as np
import pickle

from os.path import exists

if (len(sys.argv)<2):
    print('Usage: {} runtag [rmax]'.format(sys.argv[0]))
    exit(1)

rmax = 0.03
if (len(sys.argv)>2):
    rmax = float(sys.argv[2])

nbin = 80
if (len(sys.argv)>3):
    nbin = int(sys.argv[3])


# Make the file list for the snapshot sequence
#
beg_seq = 0
end_seq = 10000
file_list = []
for i in range(beg_seq, end_seq):
    file_list.append('SPL.{}.{:05d}'.format(sys.argv[1], i))
#                     ^
#                     |
#   Change this depending on your phase-space type

# Construct batches of files the particle reader.  One could use the
# parseStringList to create batches from a vector/list of files.  NB:
# a std::vector in C++ becomes a Python.list and vice versa
#
batches = pyEXP.read.ParticleReader.parseStringList(file_list, '')
xy = {}
xz = {}
yz = {}

# For a unique map with fixed-point time as a key
#
def getTime(time):
    fixedD = 100000.0;
    return int(fixedD*time+0.5)/fixedD


times = []
lower = [-rmax, -rmax, -rmax]
upper = [ rmax,  rmax,  rmax]
ngrid = [ nbin,  nbin,  nbin]

fg = pyEXP.field.FieldGenerator(times, lower, upper, ngrid)
gd = {}
rd = {}

for group in batches:

    okay = True
    for f in group:
        if not exists(f):
            okay = False
        
    if not okay: continue

    # Make the reader for the desired type.  One could probably try to
    # do this by inspection but that's another project.
    #
    reader = pyEXP.read.ParticleReader.createReader('PSPspl', group, 0, False);

    compname = 'star'
    reader.SelectType(compname)
    
    tim = getTime(reader.CurrentTime())
    gd[tim] = fg.histo2d(reader)
    rd[tim] = fg.histo1d(reader, rsize, nbins, "xy")

keys =  list(gd.keys())
print("Time[0]={}  Time[{}]={}".format(keys[0], len(keys)-1, keys[-1]))

db = {'image': gd, 'histo': rd, 'lower': lower, 'upper': upper, 'ngrid': ngrid}
file = open('imagePickle', 'wb')
pickle.dump(db, file)
file.close()
