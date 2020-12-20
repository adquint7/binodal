# binodal
An application of the gradient descent [Adam Optimization](https://arxiv.org/abs/1412.6980) to compute the equilibrium phase concentrations of a two component free-energy.

The binodal (equilibrium volume fractions), spinodal (region of phase instability), and equilibrium chemical potential is computed for any two component free-energy having exactly two inflection points within the region of volume fraction (0,1) such as the two component Flory-Huggins free-energy of polymer mixing with solvent as shown below.

![BinodalDiagram](FH2.pdf)

## Python Implementation
Starting an instance of python3 with binodal.py on the path:

`import binodal`

Currently, only the incompressible Flory-Huggins free-energy is implemented. Requests for other free-energy densities can be easily accomodated.
Compute the equilibrium volume fractions for a degree of polymerization 10 and interaction parameter 1:
`binodal.floryhugg(10,1)`

Output is stored as numpy arrays.

Binodal volume fractions:  
`binodal.vb`

Spinodal volume fractions:  
`binodal.vs`

Equilibrium chemical potential:  
`binodal.mub`

## MATLAB Implementation

`[vb,vs,mub] = binodal(f,mu,dmu)`

Input:  
`f = @(v)f(v)` : free-energy density  
`mu = @(v)mu(v)` : chemical potential (1st derivative of f)  
`dmu = @(v)dmu(v)` : 1st derivative of mu  

Output:  
vb = [dilute phase vol.frac., dense phase vol.frac.]  
vs = [dilute spinodal vol.frac., dense spinodal vol.frac]  
mub = binodal chemical potential
