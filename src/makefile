# Simple makefile to compile MPI version of Darwin_Chemostat model
#
# To run with MPI from a subdirectory, use something like > mpiexec -machinefile ../mpi.hosts ../Darwin_Chemostat
#
# Intel Fortran MPI compilier
FC = mpiifort

# Note, -heap-arrays was needed for Windows compilation, not sure if needed here.
FFLAGS = -O3 -qmkl -xHost -qopenmp -ipo -qopt-matmul -fc=ifort -o Darwin_Chemostat

# BiM routines
BIMSRC  = ./BiM/bim__NoWrite_SOLOUTmod.for
BIMSRC += ./BiM/subbim_NoLinAlg.for
FFLAGS += -I./BiM

# Thermodynamic routines
THERMOSRC = ./ThermoData/Module_Thermo_V5.1.f90

# LINPACK to LAPACK
LIN2LAPSRC = ./Linpack2Lapack/Linpack2Lapack.f90

# using hyperBOB
BOBDIR = ./BOBYQA
BOBSRC  = $(BOBDIR)/kind_module.f90
BOBSRC += $(BOBDIR)/bobyqa.f90
BOBSRC += $(BOBDIR)/hyperBOB_V1.4.f90

# using latin hypercube
LATSRC = ./Latin_Random/latin_random.f90

# Tecplot library
# Note, libtecio.a must be built. First download TecIO Library from Tecplot
# You also have to download Boost, but don't install it, just extract into a temp directory
# then modify the CMakeLists.txt as described (I used "/home/jvallino/Downloads/boost_1_82_0/") and follow readme
# I did install and use cmake. Also, you still need to include the lstdc++ library as shown below.
TECIOLIB = ./tecio/libtecio.a
FFLAGS += -I./tecio/include/
EXTRALIBS = -lstdc++

# Main program source files
DRVSRC = ./Darwin_ChemoStat_V2.2c.f90

prog: $(BIMSRC) $(LIN2LAPSRC) $(THERMOSRC) $(BOBSRC) $(LATSRC) $(DRVSRC)
	$(FC) $(FFLAGS) $(BIMSRC) $(LIN2LAPSRC) $(THERMOSRC) $(BOBSRC) $(LATSRC) $(DRVSRC) $(TECIOLIB) $(EXTRALIBS)
