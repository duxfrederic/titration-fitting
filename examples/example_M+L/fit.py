# -*- coding: utf-8 -*-
# use numpy to load the content of the textfiles to arrays
import numpy as np 
# import the needed classes
from titrationFitter import Component, System, Titration

# %% 
"""
LOADING AND SAMPLING THE DATA
"""
folder = './'
## ligand/metal ratios at the different points of the titration:
eqs = np.loadtxt(folder+'eqs.csv')
#
beg=0; end=-1; every=1
## span of the UV-VIS spectrum: 
span = np.loadtxt(folder+'span.csv')[::-1][beg:end:every]
# spanold = np.loadtxt(folder+'span.csv')[::-1][30:150:1]
spanold = np.loadtxt(folder+'span.csv')[::-1][30:120:2]
#
## values of the UV-VIS spectrum: (contained by columns in the file)
uvs  = np.loadtxt(folder+'uvs.csv')[beg:end:every,:]

volumes = np.loadtxt(folder+'volumes.csv')


#%% 
"""
CREATING THE SPECIES, ARRANGING THEM IN A SYSTEM, INITIALIAZING THE TITRATION
"""
# constructor of a component :
# Component( guess for equilibrium constant, initial concentration, {named arguments})
# (see docstring)
# if the Component is a building-block (M or L), the value of the equilibrium constant
# does not matter. (below, 0 is given as the initial guess.)

eqknown = 0
uvknown = 0

conc_ini = 0.1

M    = Component(0,    0.,  uv_known=True,  uvvis=np.zeros_like(span), eqconst_known=True, \
                 name='M', titrant=True) # the metal is the titrant.

L    = Component(0,    conc_ini,  uv_known=True,  uvvis=uvs[:,0]/conc_ini, eqconst_known=True,\
                 name='L') 

ML   = Component(10**2.97,  0.,  buildblocks=[M,L],  coeffs=[1,1], uv_known=uvknown,  
                 uvvis=uvs[:,-1]/conc_ini*volumes[-1],
                 eqconst_known=eqknown)

S    = System([M,L], [ML], conc_ini=conc_ini, span=span)

#      Titration(system, M/L ratios, experimental uv-vis spectra)
beg  = 0
end  = 11
T    = Titration(S, eqs[beg:end], uvs[:, beg:end], volumes=volumes[beg:end])
#%%
# OPTIMIZE 
T.optimize()
# plot and print the result
T.plotCurrentModel()
T.printCurrentModel()
