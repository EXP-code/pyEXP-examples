#!/usr/bin/env python
# coding: utf-8

# This script is an example of applying MSSA to the coefficients you
# created using pyEXP
# 
# This example assumes that you have run the `Part1-Coefficients`
# script to create the coefficients.
# 
# We begin by importing `pyEXP` and friends and setting the working directory.

import os
import yaml
import pyEXP
import numpy as np
import matplotlib.pyplot as plt

# In this test, I assume that you have already explored and run
# Part1-Coefficients. This script does some additional analysis and
# plotting
#
os.chdir('../Data')

# First, we create the basis This step is the same as in the previous
# notebook.  We are going to use the basis for field evaluation so we
# need the basis.

yaml_config = ""
with open('basis.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    yaml_config = yaml.dump(config)

# Construct the basis instance
#
basis = pyEXP.basis.Basis.factory(yaml_config)


# So far, this is the same as in `Part1` . . . now we read the halo
# coefficients We created the halo coefficients in the last notebook
# and stashed them in the EXP HDF5 format.

coefs = pyEXP.coefs.Coefs.factory('test_coefs')

print("Got coefs for name=", coefs.getName())


# Now, let's visualize the reconstructed fields.  We made a few test
# grids in Part 1.  This time, we'll make surface renderings:

# Now try some slices for rendering
#
times = coefs.Times()
pmin  = [-5.0, -5.0, 0.0]
pmax  = [ 5.0,  5.0, 0.0]
grid  = [  40,   40,   0]

fields = pyEXP.field.FieldGenerator(times, pmin, pmax, grid)

surfaces = fields.slices(basis, coefs)

print("We now have the following [time field] pairs")
final = 0.0
for v in surfaces:
    print('-'*40)
    for u in surfaces[v]:
        print("{:8.4f}  {}".format(v, u))
        final = v

# Print the potential image at the final time (I think there is a
# fencepost issue in this grid, no matter).
x = np.linspace(pmin[0], pmax[0], grid[0])
y = np.linspace(pmin[1], pmax[1], grid[1])
xv, yv = np.meshgrid(x, y)

# Visualize the final time slice in grid.  Obviously, we could use any 
# field we want.  Here is potential and radial force.

cont1 = plt.contour(xv, yv, surfaces[final]['potl'].transpose(), colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['potl'].transpose())
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Potential at T={}'.format(final))
plt.show()

cont1 = plt.contour(xv, yv, surfaces[final]['rad force'].transpose(), colors='k')
plt.clabel(cont1, fontsize=9, inline=True)
cont2 = plt.contourf(xv, yv, surfaces[final]['rad force'].transpose())
plt.colorbar(cont2)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Radial force at T={}'.format(final))


# Finally, let's analyze these coefficients with mSSA.  We begin by 
# generating the initial mSSA analysis and plot the run of eigenvalues

# The first step is to make a sequence of keys that identify the mSSA
# channels.
#
keylst = coefs.makeKeys([1])
print("Keys=", keylst)

# The mSSA config is a YAML structure that names the coefficients for
# mnemonic convenience and associates the name to the coefficient set
# and key list.

name = 'dark halo'
config = {name: (coefs, keylst, [])}

# This is the window length.  The largest meaningful length is half
# the size of the time series.  We'll choose that here but it *is*
# interesting to retry these analyses with different window lengths.

window = int(len(coefs.Times())/2)
npc = 10

print("Window={} PC number={}".format(window, npc))

ssa = pyEXP.mssa.expMSSA(config, window, npc)

ev = ssa.eigenvalues()

plt.plot(ev, 'o-')
plt.xlabel("Index")
plt.ylabel("EV")
plt.title("Eigenvalues by index")
plt.show()


# Next plot the principle components.

times = coefs.Times();
pc    = ssa.getPC();

rows, cols = pc.shape

for i in range(cols):
    plt.plot(times[0:rows], pc[:,i], '-', label="{:d}".format(i))

plt.xlabel('Time')
plt.ylabel('PC')
plt.legend()
plt.title("Principal components (left-singular vectors)")
plt.show()


# ### Now try a reconstruction

# In[ ]:


ssa.reconstruct([0, 1, 2, 3])

coefs.zerodata() # <---replace with reconstructed
newdata = ssa.getReconstructed()
print('newdata is a', type(newdata))


# Let's try to discover groups of channels that represent like
# dynamics using k-means analysis
#
# Each channel is the real or imaginary part (designated cosine [c] or
# sine [s]).  These could be combined and this will be done in a later
# update.  The k-means analysis illustrates the clear patterns in the
# PCs seen in the plots above: the individual features occur in pairs,
# dominated by the first two pairs as seen in the eigenvalue plot.
#

ssa.kmeans()


# Another useful diagnostic for channel correlation is the
# w-correlation matrix This is the sum of the w-correlation matrices
# for all channels.

mat = ssa.wCorrAll()
x = plt.pcolormesh(mat)
plt.title("All channels")
plt.colorbar(x)
plt.show()


# One can also look at individual channels.  Next, we retrieve and
# view the w-correlation matrices for each key individually
#
for k in keylst:
    mat = ssa.wCorr(name, k)
    # Plot it
    x = plt.pcolormesh(mat)
    plt.title('{}: {}'.format(name,k))
    plt.colorbar(x)
    plt.show()
