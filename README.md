# MEP-Chemostat
This model is being develop as part of the NSF funded project: ***Investigating the connectivity of microbial food webs using thermodynamic models, stable isotope probing and genomics*** (Award #1655552). For more details see the project website: http://ecosystems.mbl.edu/MEP-FoodWeb/

ODE solver is BiM from here: http://web.math.unifi.it/users/brugnano/BiM/

Optimization is via hyperBOB found here: https://github.com/maxEntropyProd/hyperBOB, which uses BOBYQA, also here https://github.com/jacobwilliams/PowellOpt

Latin hypercube sampling is from here: https://people.sc.fsu.edu/~jburkardt/f_src/latin_random/latin_random.html

Much of the output from the program is stored in Tecplot (https://tecplot.com/) binary files using tecio for plotting with Tecplot. The library included here was compiled for Rocky 8. The Tecplot szplt binary files can be read or converted to ASCII with the PyTecplot scripts: https://github.com/Tecplot/handyscripts

Funding support from Simons Foundation (award 549941FY22) and NSF DEB (award 1655552) 
[![DOI](https://zenodo.org/badge/956167895.svg)](https://doi.org/10.5281/zenodo.15215254)

