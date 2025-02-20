"""
an example script to port an externally-generated Hernquist basis coefficient set and fill an EXP coefficient file 
(e.g. from gala)

"""

# pyEXP
import pyEXP

# some basic imports
import os
import numpy as np
import matplotlib.pyplot as plt


def makemodel(func,M,funcargs,rvals = 10.**np.linspace(-2.,4.,2000),pfile='',plabel = '',verbose=True):
    """make an EXP-compatible spherical basis function table from a functional input

    inputs
    -------------
    func        : (function) the callable functional form of the density
    M           : (float) the total mass of the model, sets normalisations
    funcargs    : (list) a list of arguments for the density function.
    rvals       : (array of floats) radius values to evaluate the density function
    pfile       : (string) the name of the output file. If '', will not print file
    plabel      : (string) comment string, printed to the top of the file
    verbose     : (boolean)

    outputs
    -------------
    R           : (array of floats) the radius values
    D           : (array of floats) the density
    M           : (array of floats) the mass enclosed
    P           : (array of floats) the potential

    """
    R = np.nanmax(rvals)
    # query out the density values
    dvals = func(rvals,*funcargs)
    # make the mass and potential arrays
    mvals = np.zeros(dvals.size)
    pvals = np.zeros(dvals.size)
    pwvals = np.zeros(dvals.size)
    # initialise the mass enclosed an potential energy
    mvals[0] = 1.e-15
    pwvals[0] = 0.
    # evaluate mass enclosed and potential energy by recursion
    for indx in range(1,dvals.size):
        mvals[indx] = mvals[indx-1] +\
          2.0*np.pi*(rvals[indx-1]*rvals[indx-1]*dvals[indx-1] +\
                 rvals[indx]*rvals[indx]*dvals[indx])*(rvals[indx] - rvals[indx-1]);
        pwvals[indx] = pwvals[indx-1] + \
          2.0*np.pi*(rvals[indx-1]*dvals[indx-1] + rvals[indx]*dvals[indx])*(rvals[indx] - rvals[indx-1]);
    # evaluate potential (see theory document)
    pvals = -mvals/(rvals+1.e-10) - (pwvals[dvals.size-1] - pwvals)
    # get the maximum mass and maximum radius
    M0 = mvals[dvals.size-1]
    R0 = rvals[dvals.size-1]
    # compute scaling factors
    Beta = (M/M0) * (R0/R);
    Gamma = np.sqrt((M0*R0)/(M*R)) * (R0/R);
    if verbose:
        print("! Scaling:  R=",R,"  M=",M)
    rfac = np.power(Beta,-0.25) * np.power(Gamma,-0.5);
    dfac = np.power(Beta,1.5) * Gamma;
    mfac = np.power(Beta,0.75) * np.power(Gamma,-0.5);
    pfac = Beta;
    if verbose:
        print(rfac,dfac,mfac,pfac)
    # save file if desired
    if pfile != '':
        f = open(pfile,'w')
        print('! ',plabel,file=f)
        print('! R    D    M    P',file=f)
        print(rvals.size,file=f)
        for indx in range(0,rvals.size):
            print('{0} {1} {2} {3}'.format( rfac*rvals[indx],\
              dfac*dvals[indx],\
              mfac*mvals[indx],\
              pfac*pvals[indx]),file=f)
        f.close()
    return rvals*rfac,dfac*dvals,mfac*mvals,pfac*pvals


def twopower_density(r,a,alpha,beta):
    """a twopower density profile"""
    ra = r/a
    return (ra**-alpha)*(1+ra)**(-beta+alpha)


# realise an analytic basis. in this case, it's a hernquist (alpha=1,beta=4), see BT08 eq. 2.64
hernquist_rs = 15. # in kpc
total_mass   = 1.0 # this can be in any units; note that EXP enforces G=1, so the mass will be scaled.

R,D,M,P = makemodel(twopower_density,total_mass,[hernquist_rs,1,4],rvals = 10.**np.linspace(-1.,3.,2000),pfile='Hernquist_rs{}.dat'.format(hernquist_rs))


# now, let's generate a Hernquist basis inside EXP.
hernquist_config = """
id         : sphereSL
parameters :
  numr     : 2000
  rmin     : 0.1
  rmax     : 1000.0
  Lmax     : 6
  nmax     : 20
  rs       : 1.0
  modelname : Hernquist_rs15.dat
  cachename : Hernquist_rs15.cache
"""
hernquist_basis = pyEXP.basis.Basis.factory(hernquist_config)

# create a dummy array that you can modify.
hernquist_coef = hernquist_basis.createFromArray([1.],[[1.],[1.],[1.]], time=0.0)
print("Generated dummy array with dimensions {}.".format(hernquist_coef.getCoefs().shape))


# make an instance of coefs
hernquist_coefs = pyEXP.coefs.Coefs.makecoefs(hernquist_coef, 'halo')

# here is where you will loop through the times in your simulation. modify the array for each timestep

# change the time each loop through
hernquist_coef.time = 0.0

modcoefs = hernquist_coef.setCoefs()
# now you can change the array modcoefs and the changes will be reflected in hernquist_coef

# ...plug in your changes here

# the coefficients are written as imaginary now, so
# 0: Y_00, 1: Y_10, 2: Y_11, 3: Y_20, 4:Y_21, 5:Y_22, 6:Y_30, 7:Y_31, etc
# 0:0, 1:1 (1), 3:2 (2), 6:3 (3), 10:4 (4), 15:5 (5), 21:6 (6)

# following the Lowing et al. (2011) notation, the S coefficients are the real part of EXP coefficients, and the T coefficients are the imaginary part of EXP coefficients.

# add the coefficients to the makecoefs instance
hernquist_coefs.add(hernquist_coef)

# write the first time (afterwards use ExtendH5Coefs)
# note that the writer only ever writes one step, so it must be called each time
hernquist_coefs.WriteH5Coefs('outcoef.hernquist')
# now the file is outcoef.hernquist.h5
# for all subsequent steps, use hernquist_coefs.ExtendH5Coefs('outcoef.hernquist')

# end loop through times


# with these coefficients in hand, you can now create density, force, potential, etc.


# once done with looping...
# to bring the file back in,
hernquist_coefs_in = pyEXP.coefs.Coefs.factory('outcoef.hernquist.h5')

# now you can navigate the structure as
# hernquist_coefs_in.Times()
# hernquist_coefs_in.getAllCoefs()

# you can also use this generated set of coefficients to run mSSA
