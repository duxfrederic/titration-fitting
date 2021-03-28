This repository implements a simulation of a titration, for a system of two main components,
N species and N-1 equilibrium constants. Some species (composed of a combination of a metal (M) and of a ligand (L), with a given or unknown UV-VIS spectrum) are contained in a flask. More of a solution of a ligand or of the metal (the titrating agent, or titrant) is added to the flask, modifying the concentrations of the species depending on the metal/ligand concentrations ratio.  The UV-VIS spectra are recorded experimentally at each step of the titration, and the routines contained in this script try to obtain a best fit of the data by finding the equilibrium constants that will describe best the concentrations of the species. If the UV-VIS spectra of a given species is not known, its values will be added to the variables to be optimised.

### Dependencies 

- Python 3.7+
- Scipy
- Numpy
- Matplotlib

### Usage example

Examples are given in the repository of the package: https://gitlab.epfl.ch/dux/titration-fitting/.  This file is the most explicit: https://gitlab.epfl.ch/dux/titration-fitting/-/blob/master/examples/ML_ML2_ML3/tutorial.ipynb. It can be adapted to your situation, just use a [jupyter notebook](https://jupyter.org/).

A typical situation: a flask contains a solution of M.  A solution of L is added to the flask, generating the species ML, ML2 and ML3.

The spectra of L and M are known, but not those of ML and ML2 or ML3.
The initial guess for ML is given by the end result of the titration, when
L is in large excess compared to L. Initial guesses for the spectra of ML2 and ML3 are given from the knowledge of the chemistry at hand.

The required data files are listed below: 

- a list of the ligand/metal ratio (equivalents) at each step of the titration (file example/titration_data/eqs),
- the span of the UV-VIS spectrum (not mandatory: just for plotting, file example/titration_data/span),
- the UV-VIS spectra, one column per step of titration (file example/titration_data/uvs).
- The total volume at each step of the titration. 

For each step of the titration, the program computes the concentrations of the different species in solution given the initial guess of the equilibrium constants. Using the provided UV-VIS spectra (or guessed if the UV-VIS spectrum of a compound is not known), the spectrum is computed at each step of the titration. The program then performs a minimisation step on the squared differences between the modelled spectra and the experimental spectra. The iteration is repeated until a stationary point is reached. At this point the program declares convergence and the resulting estimated spectra and equilibrium constants are given.

## Remark

Say you need a complex with a metallic center and 1 to 3 ligands: ML, ML2, ML3. It is much better to define them like this:

```python
ML   =  Component(buildblocks=[M,L],  coeffs=[1,1], ......)
# same as above, except we have 1 M and 2 L.
ML2  = Component(buildblocks=[ML,L],  coeffs=[1,1], ......)

# same again, but with 1 M and 3 L.
ML3  = Component(buildblocks=[ML2,L],  coeffs=[1,1], ......)

```

The equilibrium constants are then given in terms of the intermediate species, and will all be of the same order of magnitude. This will greatly help the concentration solver and you will be a lot less frustrated with the fitting process. 
Try it in the example jupyter-notebook before moving on with your own systems. 
