# -*- coding: utf-8 -*-
# use numpy to load the content of the textfiles to arrays
import numpy as np 
from scipy.interpolate import interp1d
# import the needed classes
from titrationFitter import Component, System, Titration

# %% 
"""
LOADING AND SAMPLING THE DATA
"""
## ligand/metal ratios at the different points of the titration:
eqs = np.loadtxt('eqs')

beg = 10; end = 230; every=5
## span of the UV-VIS spectrum: 
span = np.loadtxt('wavelengths')[::-1][beg:end:every]
## values of the UV-VIS spectrum: (contained by columns in the file)
uvs  = np.loadtxt('uvs')[beg:end:every,:]

volumes = np.loadtxt('volumes')

def interpolate(x, y, xnew):
    f = interp1d(x, y, kind='cubic', bounds_error=False, fill_value='extrapolate')
    return f(xnew)
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

conc_ini = 4.96e-6

M    = Component(0,    0.,  uv_known=True,  uvvis=np.zeros_like(span), eqconst_known=True, \
                 name='M', titrant=True) # the metal is the titrant.

P    = Component(0,  conc_ini,  uv_known=True,  uvvis=uvs[:,0]/conc_ini, eqconst_known=True,\
                 name='P') 
    
H    = Component(0,  1.5e-3,  uv_known=True,  uvvis=np.zeros_like(span), eqconst_known=True,\
                 name='H') 

MP   = Component(10**6.603,  0.,  buildblocks=[M,P],  coeffs=[1,1], uv_known=uvknown,  
                 uvvis=uvs[:,-1]/conc_ini,
                 eqconst_known=eqknown)

MH   = Component(228,  0, buildblocks=[M,H], coeffs=[1,1], uv_known=True,  
                 uvvis=np.zeros_like(span), eqconst_known=True)

MP2  = Component(390000000000,  0.,  buildblocks=[M,P],  coeffs=[1,2], uv_known=uvknown,  
                  uvvis=(uvs[:,0]+uvs[:,-1])/conc_ini,
                  eqconst_known=eqknown)

S    = System([M,P], [MP], conc_ini=conc_ini, span=span)

#      Titration(system, M/L ratios, experimental uv-vis spectra)
beg  = 0
end  = 12
T    = Titration(S, eqs[beg:end], uvs[:, beg:end], volumes=volumes[beg:end])
#%%
# OPTIMIZE 
T.optimize()
# plot and print the result
T.plotCurrentModel()
T.printCurrentModel()


# print(np.log10(MP.eqconst))

