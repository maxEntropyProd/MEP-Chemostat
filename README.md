# MEP-Chemostat
This model is being develop as part of the NSF funded project: Investigating the connectivity of microbial food webs using thermodynamic models, stable isotope probing and genomics (Award #1655552). The website for the project is here: http://ecosystems.mbl.edu/MEP-FoodWeb/

ODE solver is BiM from here: http://web.math.unifi.it/users/brugnano/BiM/

Optimization is via hyperBOB found here: https://github.com/maxEntropyProd/hyperBOB, which uses BOBYQA, also here https://github.com/jacobwilliams/PowellOpt

Latin hypercube sampling is from here: https://people.sc.fsu.edu/~jburkardt/f_src/latin_random/latin_random.html

Much of the output from the program is stored in Tecplot (https://tecplot.com/) binary files using tecio for plotting with Tecplot. The library included here was compiled for Rocky 8

