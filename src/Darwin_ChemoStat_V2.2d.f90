   ! This to try in next version:
   !  -  if bacteria or consumers are ever greater than t0, then treat them as present in the SD and TP calculations? Some
   !     solutions have bacteria and grazers in a high amplitude oscillation, but at tf, them may be < initial conditions
   !
   ! Current Versions
   !  22-May-2023, V2.2
   !     -  This version adds the ability to calculate effective richness (Hill number, https://en.wikipedia.org/wiki/Diversity_index)
   !        for the bacteria (associated with the B matrix) and effective trophic levels of the consumsers
   !        from the G matrix (see Levine1980). These are calculated for all the solutions and results
   !        are stored in basefilename_stats.dat. It is also stored in the *_opt.dat file for the best solution.
   !     -  Also calculate effective richness and trophic position based on rb and rg at tf and save them in
   !        *_rStats.dat. This was added on 26-May-2023, and called V2.2b in output
   !     -  Modified the stats so that if bi or gi is <= bi_ini or gi_ini, then set the SD or TP to zero
   !        This was done 28-May-2023, so call this V2.2c
   !     -  Modified SD and TP based on B and G so that values are set to zero if bi or gi never exceed the
   !        initial conditions during the simulation. This was added 2-Jun-2023 and called V2.2d
   !        This required that the solution output be invoked, but not saved.
   !  18-May-2023, V2.1
   !     -  This adds a time window over which EP is maximized instead of using the whole simulation time;
   !        however, the output is still simulated over the specified time. New parameters t0_mep, tday_mep,
   !        where t0_mep is the start time and EP is maximized over tdays_mep days.
   !     -  Normalized B and G so they are easiler to visualize (V2.1b)
   !  8-May-2023, V2.0
   !     -  Converting this to the MEP guided trait-based version. Here, Traits are set by MEP optimization
   !        instead of randomly choosing them. Like V1.5b, the jacobian is calculated numerically, but
   !        it is done using a separate subroutine instead of having BiM do it. It's not expected that openMP
   !        will speed up that calculation much given the small food web sizes that will be run here.
   !     -  Note, data was stored using Tecplot binary files because initial model was exected to run to 1000's
   !        of B's and G's; however, food web sizes now will be much smaller, so that's not really necessary
   !        now. All traits (i.e., epsB, B, epsG and G) will be stored in an ascii file that tecplot can read
   !     -  Found bug in rg calculation: had **F_Thermo(...) instead of *F_Thermo(...) 
   !     -  Also fixed bug in calculation of ws and rg. ws for rg had a bound of nSub, but it should have been nBac,
   !        and rg had bi(i) instead of bi(j).
   ! 3-Oct-2018
   !  -  V1.5b: There was a significant performance decrease using openMP, and this 
   !     is being used to determine why and how it might be returned to non-openmp
   !     speed.  Used Intel advisor to improve vector performance via SIMD directives
   ! 26-Sep-2018
   !  -  Ver 1.5: This is similar to 1.4, but tyring OpenMP instead of MPI to improve
   !     the performance of the numerical gradient.  Aslo, instead of modifying BiM code,
   !     put the OpenMP code in the subroutine called by BiM for the Jacobian.
   ! 25-Sep-2018
   !  -  Ver 1.4.  This will examine if any code could be sped up using OpenMP in places
   !     since version 1.3 did not show much improvement by passing array parts as larger
   !     blocks in the rxnProperties() routine. (see ver 1.3 on that)
   !  -  There is considerable cost calculating the numerical gradient in the BiM code, and
   !     First try (in this version) uses MPI, which provides some speadup.  Example
   !     running nSub, nBac, nGrz (5,400,400) on 1 core ran for 9.66 min, while on 5 cores
   !     took 5.52 min, so 1.75 times faster, not 5x.
   ! 19-Sep-2018
   !  -  Ver 1.2:Using new weighting on reactions based on concentration of either
   !     substrates or prey items.  B and G are now binary matricies.
   ! 12-Sep-2018
   !  -  Ver 1.1, starting from version 1.0, but implementing TecIO storage
   ! 16-Jul-2018
   !  -  Version 1.0 of DEB Darwin Chemostat Model
   !  -  Use MPI to distribute calculation of ODE right hand side calcultions
   !     This probably won't speed things up much, and could slow it down.  I also
   !     considered using Fortran Coarray, but parts of that are still not
   !     implemented in Intel Fortran 2018, such as co_broadcast, etc.
   !  -  This version, which was never operational, was intended to use HDF5 to save
   !     data; however, it turns out that tecplot is limited to only 50 variables in
   !     a group, which will not work for what I need.  Aslo, the HDF5 files got rather
   !     large, and were slow to write.  Instead, use TecIO
   module realPrec
      ! this sets the precision of reals
      implicit none
      save
      integer, parameter:: sp = kind(1.0E0)
      integer, parameter:: dp = kind(1.0D0)
      integer, parameter:: mp = dp ! In theory, just change defintion of mp to change precision, 
                                   ! but it will likely fail for numerous reasons.
   end module realPrec
           
   module globalVars
      ! This module contains global parameters
      use realPrec
      use thermoData
      use iso_c_binding   ! needed for declaring C variables for tecIO   
      implicit none
      include 'tecio.f90' ! interfaces for tecIO functions
      save
      ! Fixed parameters.  Changing these would require altering the program, so no need to include them in *.ini file
      integer, parameter:: nChm = 5 ! number of chemcial species (i.e., DO, H3PO4, H2CO3, etc)
      character(20), parameter:: version = '2.2d, 2-Jun-2023'      
      
      ! Parameters to be read from basefilename.ini
      character(80) basefilename      
      integer(8) nSub ! number of substrates
      integer(8) nbac ! number of bacteria
      integer(8) ngrz ! number of grazers/preditors
      namelist /dimens/ nSub, nBac, nGrz ! the /dimens/ list should occur first in the *.ini file to allow allocation of matricies
      real(mp), allocatable:: CHO(:,:) ! alpha, beta, gamma for each substrate (i.e., chemical formula)
      character(maxName), allocatable:: subNames(:) ! Names for substrates. See thermoData for maxName size
      
      real(mp) absZero ! Number to use as zero
      real(mp) T_K ! temperature (K)
      real(mp) pH
      real(mp) is ! ionic strength (M)
      real(mp) kappa ! half saturation constant (mmol/m3)
      real(mp) nuStar ! max growth rate (1/d)
      real(mp) delPsi ! Cell potential (V)
      real(mp) alf_Bio, bet_Bio, gam_Bio, del_Bio, zet_Bio ! C, H, O, N, P composition of Bac and Grz
      real(mp) V_L, F_L, V_G, F_G ! Bioreactor liquid volume (m3), flow rate (m3/d) and gas volume (m3) and flow rate (m3/d)
      real(mp) kL_o2, kL_co2 ! liquid-side mass transfer coefficient (piston velocity) for oxygen and co2, repsectively (m/d)
      real(mp) area_GL ! area of gas-liquid interface, including bubbles if present (m2)
      real(mp) nh3_hold ! Value to set nh3 concentration to (it's not currently a state variable) (mmol/m3)
      integer iSeed ! Used to initialize the random_number routine
      
      ! Intitial conditions
      real(mp) t0, tDays ! the start time and number of days to run for ODE integration (d)
      namelist /parameters/ t0, tDays
      real(mp), allocatable:: sj_f(:)  ! substrates. Note, allocatable arrays in namelist are now allowed in Fortran 2003.
      real(mp), target:: chm_f(nChm)   ! Chemical concentrations
      real(mp) bi_ini, bi_f            ! for now, all bacteria are set to the same initial condition and feed concentration
      real(mp) gi_ini, gi_f            ! likewise for grazeres
      real(mp) t0_mep ! the time (day) that entropy production maximization begins
      real(mp) tf_mep ! the end time of MEP optimization (this is calculated)
      real(mp) tDays_mep ! the number of days over which entropy production is maximized
      namelist /parameters/ t0_mep, tDays_mep
            
      ! Parameters associated with problem
      namelist /parameters/ absZero, T_K, pH, is
      namelist /parameters/ kappa, nuStar, delPsi
      namelist /parameters/ alf_Bio, bet_Bio, gam_Bio, del_Bio, zet_Bio
      namelist /parameters/ V_L, F_L, V_G, F_G
      namelist /parameters/ kL_o2, kL_co2 
      namelist /parameters/ area_GL
      namelist /parameters/ subNames
      namelist /parameters/ CHO, nh3_hold
      namelist /parameters/ sj_f, chm_f, bi_ini, bi_f, gi_ini, gi_f
      namelist /parameters/ iSeed
      real(mp):: sdtpMin = 1.1 ! This factor multiplies bi_ini or gi_ini to determined of bi or gi are important in a simulation
      namelist /parameters/ sdtpMin ! can just use the default value
      logical onlyCheckSoln
      
      ! ODE constrain
      real(mp) tf ! calculated end time.
      real(mp) maxODEtime ! maximum time allowed to solve the ODE = (tf-t0)/minComFac in days
      real(mp) ODE_t0 ! time at the start the ODE solution
      
      ! BiM parmeters
      real(mp) rtol, atol ! relative and absolute tolerances
      real(mp) hmax_BiM ! maximum step size.  = (tend-t0)/8 if set to 0
      integer maxstep_BiM ! maximum number of steps. Default 100000 if set to zero
      real(mp) maxattempts ! number of retrys    
      integer ompThreads ! Number of OMP threads to use for Jacobian calculation
      real(mp) minCompFac ! if a process takes longer than (tf-t0)/minCompFac, then it is terminated 
      namelist /BiMparams/ rtol, atol, maxattempts, maxstep_BiM, hmax_BiM, ompThreads, minCompFac

      ! hyperBOB parameters
      logical optimize ! if set to false then traits values are randomly set and a simulation is run. Mostly for testing
      real(mp) rhobeg   ! initial and final values of a trust region radius
      real(mp) rhoend 
      integer iprint   ! controls amount of printing (0, 1, 2 or 3)
      integer maxfun   ! maximum number of calls to CALFUN
      integer fcnUpdate ! After every fcnUpdate PDE integrations, infor is printed out.      
      integer seed_bobyqa ! used by hyperBOB for random hypercube. All MPI processes use same value
      namelist /hyperBOBparams/ optimize, rhobeg, rhoend, iprint, maxfun, fcnUpdate, seed_bobyqa
      integer fcnCalls, fcnCallsTot ! Not a parameter, but a counter
      integer ODEfailed, ODEfailedTot, ODEtimedOut, ODEtimedOutTot ! number of times the PDE integration failed or timedout
      real(mp):: fcnMax = 0.0_mp 
      real(mp), allocatable:: traitsL(:), traitsU(:)
            
      ! tecIO parameters
      ! nIOpts: number of points to store before saving in tecplot binary file.  These are saved as zones. Making nIOpts
      ! larger decreases the number of zones that will need to be plotted, and will reduce the number of duplicate
      ! points needed for zone overlaps, but the more internal memory that will be used.
      integer(8) nIOpts 
      integer:: debug_IO = 0 ! debuging on: 1, off: 0
      logical:: flushData = .false. ! is true, then data is flushed to temp files at each write.
      namelist /tecIOparams/ nIOpts, debug_IO, flushData 
      
      ! Global variables not set in *.ini file, but rather calculated
      integer myRank, noProc ! rank of MPI process and number of processes available
      ! Most output is via Tecplot, but a few are just ascii
      integer inpUnit ! units given to input and output streams  
      integer iunit_allOpt ! file to store the optimum traits from all processes
      integer iunit_opt    ! parameters found by optimization
      integer iunit_stats  ! Contains the stats (substrate diversity, trophic position) for all simulations
      integer iunit_rStats ! like above, but calculated from rb and rg at tf. 
      real(mp) WTime_it0 ! used as the global start time      
      real(mp):: aveODEtime = 0._mp ! Average PDE integration time 
      
      integer nState  ! number of state variables used by BiM
      integer nTraits ! number of traits or adjustable parameters, calculated in sysSetup()      
      real(mp), allocatable:: aCs(:), bCs(:), aAs(:), bAs(:) ! reaction stoichiometric coef. for bacterial substrate uptake
      real(mp) aCg, bCg ! reaction stoichiometric coef. for grazer predation.
      real(mp), allocatable:: epsB(:), epsG(:) ! theromodynamic efficiencies for Bac and Grz reactions
      real(mp), allocatable:: B(:,:) ! nbac x nSub matrix that determines which bacteria each which substrate
      real(mp), allocatable:: G(:,:) ! nGrz x (nBac+nGrz) matrix that determines connectivity between prey and predators
      type solInfo
         ! used to pass information about ODE integration
         integer fail ! = 1 if ODE integration failed; = 2 if time exceeded
         integer idid ! idid value from BACOLI
         real(mp) t0  ! Integration time where failure occured
      end type solInfo      
      
      ! Thermodynamic variables
      real(mp), allocatable:: drG0_rbC(:), drG0_rbA(:) ! Standard free energies of reaction for bacterial catabolic and anabolic reactions
      real(mp)                drG0_rgC,    drG0_rgA    ! Standard free energies of reaction for predator catabolic and anabolic reactions
      real(mp), allocatable:: drG_rb(:,:) ! Gibbs free energies of reaction for bacteria
      real(mp), allocatable:: drG_rg(:,:) ! Gibbs free energies of reaction for predators
      real(mp), allocatable:: rb(:,:), rg(:,:) ! reaction rates
      real(mp), allocatable:: sig_rb(:,:), sig_rg(:,:) ! entropy production associated with rb and rg reactions
      
      ! State variable names
      real(mp), allocatable:: sj(:)   ! concentration of added substrates
      real(mp), allocatable, target:: chm(:) ! vector that holds concentration of chemical species
      real(mp), pointer:: o2, po2, h2co3, pco2, h3po4 ! pointers to chm to allow easy use of names
      real(mp), pointer:: o2_f, po2_f, h2co3_f, pco2_f, h3po4_f ! pointers to chm in feed to allow easy use of names
      integer io2, ipo2, ih2co3, ipco2, ih3po4 ! Set in setStateNames to locate the variables on state vector used by BiM
      real(mp) nh3    ! Not currently being used, but is set by nh3_hold parameter
      real(mp), allocatable:: bi(:) ! Concentration of bacteria
      real(mp), allocatable:: biF(:) ! Flag set to 1. if bi exceeds bi_ini
      real(mp), allocatable:: gi(:) ! Concentration of predators
      real(mp), allocatable:: giF(:) ! Flag set to 1. if gi exceeds gi_ini
      real(mp) sumSigmaDot ! entropy production rate (J/K/d) summed over all reactions.
      real(mp) sumSigma ! Integrated entropy production summed over all reactions (J/K)
      real(mp), allocatable:: wS(:) ! this is used to determine weights for substate or prey consumption in reactions
      
      ! Variables associated with tecIO for output.  These are only allocated in the root process
      integer(8) :: nIO_cnt = 0 ! Used as a counter to fill tecplot time series zone.  Once nIOpts is reached, a zone is written for tecplot output.
      ! Output file handles
      type(c_ptr) :: sub_file = C_NULL_PTR, chm_file = C_NULL_PTR, Bac_file = C_NULL_PTR, Grz_file = C_NULL_PTR
      type(c_ptr) :: rb_file = C_NULL_PTR, rg_file = C_NULL_PTR, sigma_file = C_NULL_PTR 
      type(c_ptr) :: epsB_file = C_NULL_PTR, epsG_file = C_NULL_PTR, B_file = C_NULL_PTR, G_file = C_NULL_PTR
      type(c_ptr) :: gridFile = C_NULL_PTR ! only used for solution only writing.  Use null for fileType = 0 or 1
      ! Various variables assocated with tecIO
      integer(c_int32_t) :: fileFormat = 1 ! 0 for plt, 1 for szplt
      integer(c_int32_t) :: fileType = 0 ! 0=full, 1=grid, 2=solution
      integer(c_int32_t) :: debug ! turn debuging on for testing
      integer(c_int32_t) :: defDataType  = 1 ! 1: float, 2: double, 3: 32-bit integer, 4: 16-bit int, 5: 8-bit unsigned int
      integer(c_int32_t) :: strandID = 1 ! not really using this.
      integer(c_int32_t) :: faceNeighborMode = 0
      integer(c_int32_t) :: shareConnectivityFromZone = 0
      integer(c_int64_t) :: numFaceConnections = 0
      integer(C_int32_t) :: zonesToRetain(1) = 0
      integer(c_int32_t) :: dataType
      integer(c_int32_t) :: zone
      integer(c_int64_t) :: iMax, jMax, kMax
      real(c_double) :: solTime
      integer(c_int32_t), allocatable :: varTypes(:)
      integer(c_int32_t), allocatable :: shareVarFromZone(:)
      integer(c_int32_t), allocatable :: valueLocations(:)
      integer(c_int32_t), allocatable :: passiveVarList(:)
      character(:), allocatable:: varList
      integer ioErr
      
      ! Space to hold solution over the nIOpts points in a zone
      real(c_float), allocatable:: time_IO(:)      ! time
      real(c_float), allocatable:: sub_IO(:,:)     ! subsrates
      real(c_float), allocatable:: chm_IO(:,:)     ! Chemistry
      real(c_float), allocatable:: Bac_IO(:,:)     ! Bacteria
      real(c_float), allocatable:: Grz_IO(:,:)     ! Grazers
      real(c_float), allocatable:: sigma_IO(:,:)     ! summed entropy production rate and integrated entropy
      real(c_float), allocatable:: mat_IO(:,:)     ! this is dummy array used for writing tecIO output
      real(c_float), allocatable:: iIndx(:,:), jIndx(:,:)   ! Used to save indexs values for matrix saves   
      
      ! MPI parameters for writing integration errors
      character(len=2), parameter:: CRLF = char(13)//char(10) ! used by MPI writes for error output      
      integer mpiFHerr     ! Unit for MPI writing errors to for BiM, and is set by MPI_FILE_OPEN in the initialization routine.            
      
      ! For OpenMP, need to make the variables used by ODEs threadprivate
      !$OMP threadprivate(sj, chm, bi, gi, sumSigma)
      !$OMP threadprivate(wS, drG_rb, rb, sig_rb, drG_rg, rg, sig_rg, sumSigmaDot)
      !$OMP threadprivate(o2, po2, h2co3, pco2, h3po4)      
  contains
      subroutine initialize()
         ! gets base filename, openFiles for input and output, reads parameters and allocates some space. 
         use mpi
         implicit none
         integer narg, i, j, ioerr, ioE, mpierr
         character (len=80) chr80

         ! Initialize thermoData module
         call initThermo()
         
         ! First try getting basefilename from commandline
         ! if that fails, then have processes zero get name and broadcast it to other processes 
         narg = command_argument_count () ! see if a filename has been included in command line
         if (narg == 1) then
            call get_command_argument(narg,value=basefilename)
         else 
            if (myRank == 0) then
               write(*,'(a,$)') 'Enter parameter file base name, no extension: '
               read(*,*) basefilename
            end if
            call MPI_BCAST(basefilename,80,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
         end if
             
         ! open a readonly unit to the inp file. This has to be readonly so that each process can open
         ! it without conflict. This file must be available for all processes.
         open(newunit=inpUnit,file=trim(basefilename)//'.inp',action='read',status='old',iostat=ioerr)   
         if (ioerr /= 0) then 
            write(6,'(a,i5,a)') 'Process ',myRank,' could not open INP file' ! this may or maynot print.
         else
            ! get dimension of problem, and allocate space
            read(inpUnit,nml=dimens)
            ! set the number of traits needed
            nTraits = nBac & ! epsB
                    + nBac*nSub & ! B matrix
                    + nGrz & ! epsG
                    + nGrz*(nBac + nGrz) ! G matrix
            ! allocate space
            allocate( traitsL(nTraits), traitsU(nTraits) ) ! used for easy non-dimensionalization of trait vector            
            allocate( subNames(nSub) )
            allocate( CHO(nSub,3) )
            allocate( B(nBac,nSub), G(nGrz,nGrz+nBac) )
            allocate( epsB(nBac), epsG(nGrz) )
            allocate( aCs(nSub), bCs(nSub), aAs(nSub), bAs(nSub) )
            allocate( drG0_rbC(nSub), drG0_rbA(nSub) )
            allocate( sj_f(nSub) )
            
            ! read in model parameters
            read(inpUnit,nml=parameters)
            ! read in hyperBOB parameters
            read(inpUnit,nml=hyperBOBparams)            
            ! read tecIO parameters
            read(inpUnit,nml=tecIOparams)
            ! get BiM parameters
            read(inpUnit,nml=BiMparams)            
            close(inpUnit) 
            if (nIOpts < 2) then
               if (myRank == 0) write(6,'(/a/)') 'Error:: nIOpts must be >= 2. Aborting.'
               call cleanUp
            end if
            
            if (myRank == 0) then
               ! Allocate space to store data temporarily for tecplot zones, but only for the root process
               allocate( time_IO(nIOpts) )         ! time
               allocate( sub_IO(nSub,nIOpts) )     ! subsrates
               allocate( chm_IO(nChm,nIOpts) )     ! Chemistry
               allocate( Bac_IO(nBac,nIOpts) )     ! Bacteria
               allocate( Grz_IO(nGrz,nIOpts) )     ! Grazers
               allocate( sigma_IO(2,nIOpts) )      ! summed entropy production rate and summed sigma integrated over time
               allocate( mat_IO(max(nBac,nGrz), max(nSub,nBac+nGrz)) ) ! Allocate space as large as will ever be needed.
               allocate( iIndx(max(nBac,nGrz), max(nSub,nBac+nGrz)), jIndx(max(nBac,nGrz), max(nSub,nBac+nGrz)) ) ! Make these as big as will ever be needed
               ! fill out index matrices, they never change
               ! Note, some memory could be saved here by using just two vectors, then using a loop for tecplot write.
               do i=1,max(nBac,nGrz)
                  do j=1,max(nSub,nBac+nGrz)
                     iIndx(i,j) = real(i)
                     jIndx(i,j) = real(j)
                  end do
               end do               
            end if
         end if
         
         ioerr = abs(ioerr) ! remove negative sign, if present
         ! sum up ioerr across processes and place in ioE in ALL processes.  ioE should be 0 if no errors occured.
         call MPI_ALLREDUCE( ioerr, ioE, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)   
         ! check for error in opening the PRM file by any of the processes.
         if (ioE /= 0) then
            ! an error occured reading the prm file by one or more processes
            if (myRank == 0) then
               write(6,'(/a/)') 'Error opening or reading '//trim(basefilename)//'.inp file. Aborting.'
            end if
            call cleanUp
         end if         
            
         ! BiM file for error output for all processes
         ! open up file using MPI for error output from each process.  This is where mpiFHerr gets set.
         ! use MPI_MODE_WRONLY for write only, creat it if necessary, and for sequential access
         ! However, there is no easy way to clear out a file if it already exists, so open with close_on_delete, to get rid of trash, then reopen
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL+MPI_MODE_DELETE_ON_CLOSE, &
                  MPI_INFO_NULL, mpiFHerr, mpiErr)
         call MPI_FILE_CLOSE(mpiFHerr, mpiErr) ! This will delete the file.    
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basefilename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL, &
                  MPI_INFO_NULL, mpiFHerr, mpiErr)
         ! write a header to the file
         if (myRank == 0) then
            chr80 = 'Errors associated with BiM integration:' 
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(chr80)//CRLF, len_trim(chr80)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
         end if   
         
         ! Initialize pointer assignments and where the chm variables are located in the BiM solution vector.
         io2    = nSub+1
         ipo2   = nSub+2
         ih2co3 = nSub+3
         ipco2  = nSub+4
         ih3po4 = nSub+5
         nh3   = nh3_hold ! this is not a state variable yet, but it's value is needed for kinetics and thermodynamics.
         ! pointers for feed concentrations
         o2_f    => chm_f(1)
         po2_f   => chm_f(2)
         h2co3_f => chm_f(3)
         pco2_f  => chm_f(4)
         h3po4_f => chm_f(5)         
                   
         ! Calculated variables that do not change over time
         nState = nChm + nSub + nBac + nGrz + 1 ! Number of state variables.  Last variable is entropy production
         ! Reaction stoichiometric coefficients
         ! Bacteria catabolic and anabolic reactions
         aCs(:) = -(-4._mp*CHO(:,1) - CHO(:,2) + 2._mp*CHO(:,3))/(4._mp*CHO(:,1))
         bCs(:) = -1._mp + CHO(:,2)/(2._mp*CHO(:,1))
         ! opps, typed in ones for reduced substrate, which is not being used, but keep in case it is later used.
         !aAs(:) = ( -(bet_Bio - 2._mp*gam_Bio - 3._mp*gam_Bio + 5._mp*zet_Bio)/alf_Bio + (CHO(:,2) - 2._mp*CHO(:,3))/CHO(:,1) )/4._mp
         !bAs(:) = ( (-bet_Bio + 3._mp*(del_Bio + zet_Bio))/alf_Bio + CHO(:,2)/CHO(:,1) )/2._mp
         aAs(:) = ( CHO(:,1)*(bet_Bio - 2._mp*gam_Bio - 3._mp*del_Bio + 5._mp*zet_Bio) - alf_Bio*(CHO(:,2) - 2._mp*CHO(:,3)) ) &
                / (CHO(:,1)*(4._mp*alf_Bio + bet_Bio - 2._mp*gam_Bio - 3._mp*del_Bio + 5._mp*zet_Bio))
         bAs(:) = ( CHO(:,1)*(-3._mp*bet_Bio + 2._mp*gam_Bio + 9._mp*del_Bio + zet_Bio) + CHO(:,2)*(3._mp*alf_Bio - gam_Bio + 4._mp*zet_Bio) &
                + CHO(:,3)*(-2._mp*alf_Bio + bet_Bio - 3._mp*(del_Bio + zet_bio)) ) &
                / (CHO(:,1)*(4._mp*alf_Bio + bet_Bio - 2._mp*gam_Bio - 3._mp*del_Bio + 5._mp*zet_Bio))
         ! Grazer catabolic reaction
         aCg = -(-4._mp*alf_Bio - bet_Bio + 2._mp*gam_Bio + 3._mp*del_Bio - 5._mp*zet_Bio)/(4._mp*alf_Bio)
         bCg = -(2._mp*alf_Bio - bet_Bio + 3._mp*del_bio + 3._mp*zet_Bio)/(2._mp*alf_Bio) 
         
         ! Standard free energy calculcations, as T, is and pH are assume to be more or less constant in the chemostat.
         ! For the bacterial catabolic and anaboic reactions
         do i=1,nSub
            drG0_rbC(i) = dGf0n('h2co3',T_k, is, pH) + bCs(i)*dGf0n('h2o',T_k, is, pH) &
                        - ( dGf0n(subNames(i),T_k, is, pH)/CHO(i,1) + aCs(i)*dGf0n('o2aq',T_k, is, pH) )
            drG0_rbA(i) = (1._mp - aAs(i))*dGf0n('yeast',T_k, is, pH)/alf_Bio + aAs(i)*dGf0n('h2co3',T_k, is, pH) + bAs(i)*dGf0n('h2o',T_k, is, pH) &
                        - ( dGf0n(subNames(i),T_k, is, pH)/CHO(i,1) + (1._mp - aAs(i))*(del_Bio*dGf0n('ammonia',T_k, is, pH) + zet_Bio*dGf0n('pi',T_k, is, pH))/alf_Bio )
         end do  
         ! Standard free energies of reaction for predator reactions
         drG0_rgC = dGf0n('h2co3',T_k, is, pH) + del_Bio*dGf0n('ammonia',T_k, is, pH)/alf_Bio + zet_Bio*dGf0n('pi',T_k, is, pH)/alf_Bio + bCg*dGf0n('h2o',T_k, is, pH) &
                  - ( dGf0n('yeast',T_k, is, pH)/alf_Bio + aCg*dGf0n('o2aq',T_k, is, pH) )
         drG0_rgA = 0._mp ! Assume composition is the same, so only concentrations differences alter reaction free eneries calculated below.
         
         if (myRank /=0 ) return ! only the root process will write to files.
         debug = debug_IO ! turn debugging on or off (1: on, 0: off)  In *.ini file
         call tecIOsetup() ! setup the tecplot output files
         ! Optimal parameter values found
         open(newunit=iunit_opt,   file=trim(basefilename)//'_opt.dat'  ,status='unknown')         
         return
      end subroutine initialize
      
      subroutine initThreads()
         ! This routine is needed for OpenMP, because these are decleared threadprivate and must
         ! be allocated within a parallel region; otherwise, threads other than the master will not have them
         allocate( drG_rb(nBac,nSub), drG_rg(nGrz,nBac+nGrz) )
         allocate( rb(nBac,nSub), rg(nGrz,nBac+nGrz) )
         allocate( sig_rb(nBac,nSub), sig_rg(nGrz,nBac+nGrz) )
         allocate( sj(nSub), bi(nBac), gi(nGrz) )
         allocate( biF(nBac), giF(nGrz) )
         allocate( chm(nChm) )
         allocate( wS(max(nBac,nGrz)) )
      
         ! Initialize pointer assignments and where the chm variables are located in the BiM solution vector.
         o2    => chm(1)
         po2   => chm(2)
         h2co3 => chm(3)
         pco2  => chm(4)
         h3po4 => chm(5)
         return
      end subroutine initThreads
      
      subroutine initializeState(state)
         ! This initializes the state vector
         real(mp), intent(out):: state(:)
         ! local declarations
         integer i
         
         ! Substrates: initialize with what is in the feed
         do i=1,nSub
            state(i) = sj_f(i)
         end do
         ! Chemistry: initialize with what is in the feed
         do i=1,nChm
            state(nSub+i) = chm_f(i)
         end do
         ! Bacteria: Set them all to what the bi_ini paramters is set to
         state(nSub+nChm+1:nSub+nChm+nBac) = bi_ini! 
         ! Grazers
         state(nSub+nChm+nBac+1:nSub+nChm+nBac+nGrz) = bi_ini! 
         ! Entropy production rate
         state(nState) = 0.           
         return
      end subroutine initializeState   
      
      subroutine setupEpsBG()
         ! This routine is used to populate epsB, epsG, B and G
         ! In this version (1.2), B and G are binary (either 0 or 1)
         ! For now, just generate random bits.
         integer i
         integer, allocatable:: seed(:) ! used to set the random_seed.  
      
         ! First setup the random number generator
         call RANDOM_SEED (size=i) ! get number of integers random_seed uses
         allocate ( seed(i) )
         seed = iseed ! set a value based in iseed in the parameters file.  This allows a run to be repeated
         call random_seed (put=seed(1:i))
         deallocate ( seed )
         
         ! Generate random values for epsB and epsG
         call RANDOM_NUMBER (epsB)
         call RANDOM_NUMBER (epsG)
         ! random_number can produce zero values, so fix those up to be >= absZero
         where (epsB < absZero) epsB = absZero
         where (epsG < absZero) epsG = absZero
         
         ! For now, just generate random connection matrices
      
         call RANDOM_NUMBER(B)
         call RANDOM_NUMBER(G)
         ! Now round number either to 0 or 1
         B = anint(B)
         G = anint(G)
         return
      end subroutine setupEpsBG      
      
      subroutine tecIOsetup()
         ! Sets up tecIO output files using tecFileWriterOpen, where the last argument is the file handle  
         integer i, j
         character (len=2000) string
         character (len=20) iStr, jStr         

         ! Setup up tecplot binary files.  Each file must have the same variables, so it is easier just to setup separate files
         dataType = defDataType
         ! Substrates
         varList = 'Time (d)'
         do i=1,nSub
            varList = trim(varList)//', '//trim(subNames(i))//' (uM)'
         end do         
         ioErr = tecFileWriterOpen(trim(basefilename)//'_sub'//C_NULL_CHAR, & ! data set filename
            'Substrates data'//C_NULL_CHAR, & ! data set title
            trim(varList)//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, sub_file)
         ioErr = tecFileSetDiagnosticsLevel(sub_file, debug)
         
         ! open a file for chemistry data
         ioErr = tecFileWriterOpen(trim(basefilename)//'_chm'//C_NULL_CHAR, & ! data set filename
            'Chemistry data'//C_NULL_CHAR, & ! data set title
            'Time (d), O2 (uM), pO2 (mbar), H2CO3 (uM), pCO2 (mbar), H3PO4 (uM)'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, chm_file)
         ioErr = tecFileSetDiagnosticsLevel(chm_file, debug)
         
         ! Bacteria
         varList = 'Time (d)'
         do i=1,nBac
            write(iStr,*) i
            varList = trim(varList)//', b('//trim(adjustL(iStr))//') (uM)'
         end do         
         ioErr = tecFileWriterOpen(trim(basefilename)//'_Bac'//C_NULL_CHAR, & ! data set filename
            'Bac data'//C_NULL_CHAR, & ! data set title
            trim(varList)//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, Bac_file)
         ioErr = tecFileSetDiagnosticsLevel(Bac_file, debug)
         
         ! Grazers
         varList = 'Time (d)'
         do i=1,nGrz
            write(iStr,*) i
            varList = trim(varList)//', g('//trim(adjustL(iStr))//') (uM)'
         end do         
         ioErr = tecFileWriterOpen(trim(basefilename)//'_Grz'//C_NULL_CHAR, & ! data set filename
            'Grz data'//C_NULL_CHAR, & ! data set title
            trim(varList)//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, Grz_file)
         ioErr = tecFileSetDiagnosticsLevel(Grz_file, debug)
         
         ! Bacterial reaction rates and entropy production
         ioErr = tecFileWriterOpen(trim(basefilename)//'_rb'//C_NULL_CHAR, & ! data set filename
            'rb and sigmaDot data'//C_NULL_CHAR, & ! data set title
            'nBac, nSub, rb (uM/d), sigma rb (J/K/d)'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, rb_file)
         ioErr = tecFileSetDiagnosticsLevel(rb_file, debug)
         
         ! Grazer reaction rates and associate entropy production
         ioErr = tecFileWriterOpen(trim(basefilename)//'_rg'//C_NULL_CHAR, & ! data set filename
            'rg and sigmaDot data'//C_NULL_CHAR, & ! data set title
            'nGrz, nBac+nGrz, rg (uM/d), sigma rg (J/K/d)'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, rg_file)
         ioErr = tecFileSetDiagnosticsLevel(rg_file, debug)
         
         ! Summed sigma dot for all reactions
         ioErr = tecFileWriterOpen(trim(basefilename)//'_sigma'//C_NULL_CHAR, & ! data set filename
            'summed sigma dot'//C_NULL_CHAR, & ! data set title
            'Time (d), sigmaDot (J/K/d), sumSigma (J/K)'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, sigma_file)
         ioErr = tecFileSetDiagnosticsLevel(sigma_file, debug)
         
         ! eps vector for bacteria, this is only saved once (does not change with time)
         ioErr = tecFileWriterOpen(trim(basefilename)//'_epsB'//C_NULL_CHAR, & ! data set filename
            'eps for bacteria'//C_NULL_CHAR, & ! data set title
            'nBac, epsB'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, epsB_file)
         ioErr = tecFileSetDiagnosticsLevel(epsB_file, debug)
         
         ! eps vector for grazers, this is only saved once (does not change with time)
         ioErr = tecFileWriterOpen(trim(basefilename)//'_epsG'//C_NULL_CHAR, & ! data set filename
            'eps for grazers'//C_NULL_CHAR, & ! data set title
            'nGrz, epsG'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, epsG_file)
         ioErr = tecFileSetDiagnosticsLevel(epsG_file, debug)
         
         ! B matrix for bacteria, this is only saved once (does not change with time)
         ioErr = tecFileWriterOpen(trim(basefilename)//'_B'//C_NULL_CHAR, & ! data set filename
            'B matrix for bacteria'//C_NULL_CHAR, & ! data set title
            'nBac, nSub, B'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, B_file)
         ioErr = tecFileSetDiagnosticsLevel(B_file, debug)
         
         ! G matrix for grazers, this is only saved once (does not change with time)
         ioErr = tecFileWriterOpen(trim(basefilename)//'_G'//C_NULL_CHAR, & ! data set filename
            'G matrix for grazers'//C_NULL_CHAR, & ! data set title
            'nGrz, nBac+nGrz, G'//C_NULL_CHAR, & ! Names of variables
            fileFormat, fileType, dataType, gridFile, G_file)
         ioErr = tecFileSetDiagnosticsLevel(G_file, debug)         
         return
      end subroutine tecIOsetup
   
      subroutine saveRxnRates(t)
         ! This routine saves the reaction rate and sigma dot matricies for bacteria
         ! and grazers.  This is done at each solution output, where each timepoint
         ! is saved in its own zone.
         ! NOTE, the reaction rates and entropy production must be calculated before calling this routine
         ! that is, the folling should have been called
         ! call setStateNames(x) ! copy state vector to useful names, as well as prevent numbers < absZero
         ! call rxnProperties(t) ! calculated as the reaction rates, etc.
      
         real(mp), intent(in):: t  ! Time solution is to be saved for
         ! Local declarations
         character(25) string
         
         solTime = t ! sets tecplot solution time variable
         write(string,*) solTime
         
         ! Save bacteria reaction rates and sigma dot
         iMax = nBac
         jMax = nSub
         kMax = 1
         ioErr = tecZoneCreateIJK(rb_file, 'time = '//trim(adjustL(string))//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
            
         ! write data to zone, in block format
         ioErr = tecZoneSetUnsteadyOptions(rb_file, & ! Here solTime is used.
            zone, solTime, strandID)
         ioErr = tecZoneVarWriteFloatValues(rb_file, &
            zone, 1, 0, nBac*nSub, iIndx(1:nBac,1:nSub))
         ioErr = tecZoneVarWriteFloatValues(rb_file, &
            zone, 2, 0, nBac*nSub, jIndx(1:nBac,1:nSub))
         mat_IO(1:nBac,1:nSub) = rb     ! copy to c-float number (i.e., single precision)
         ioErr = tecZoneVarWriteFloatValues(rb_file, &
            zone, 3, 0, nBac*nSub, mat_IO(1:nBac,1:nSub))     
         mat_IO(1:nBac,1:nSub) = sig_rb     ! copy to c-float number (i.e., single precision)         
         ioErr = tecZoneVarWriteFloatValues(rb_file, &
            zone, 4, 0, nBac*nSub, mat_IO(1:nBac,1:nSub))            
         if (flushData) ioErr = tecFileWriterFlush(rb_file, 0, zonesToRetain) ! flush data in memory to disk         
         
         ! Save grazer reaction rates and sigma dot
         iMax = nGrz
         jMax = nBac+nGrz
         kMax = 1
         ioErr = tecZoneCreateIJK(rg_file, 'time = '//trim(adjustL(string))//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
            
         ! write data to zone, in block format
         ioErr = tecZoneSetUnsteadyOptions(rg_file, & ! Here solTime is used.
            zone, solTime, strandID)
         ioErr = tecZoneVarWriteFloatValues(rg_file, &
            zone, 1, 0, nGrz*(nBac+nGrz), iIndx(1:nGrz,1:nBac+nGrz))
         ioErr = tecZoneVarWriteFloatValues(rg_file, &
            zone, 2, 0, nGrz*(nBac+nGrz), jIndx(1:nGrz,1:nBac+nGrz))
         mat_IO(1:nGrz,1:nBac+nGrz) = rg     ! copy to c-float number (i.e., single precision)
         ioErr = tecZoneVarWriteFloatValues(rg_file, &
            zone, 3, 0, nGrz*(nBac+nGrz), mat_IO(1:nGrz,1:nBac+nGrz))     
         mat_IO(1:nGrz,1:nBac+nGrz) = sig_rg     ! copy to c-float number (i.e., single precision)         
         ioErr = tecZoneVarWriteFloatValues(rg_file, &
            zone, 4, 0, nGrz*(nBac+nGrz), mat_IO(1:nGrz,1:nBac+nGrz))            
         if (flushData) ioErr = tecFileWriterFlush(rg_file, 0, zonesToRetain) ! flush data in memory to disk
         return
      end subroutine saveRxnRates      

      subroutine saveEpsBG()
         ! This routine saves epsB, epsG, B and G.  This is only done once.
      
         ! Write the epsB vector
         iMax = nBac
         jMax = 1
         kMax = 1
         ioErr = tecZoneCreateIJK(epsB_file, 'epsB'//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format
         ioErr = tecZoneVarWriteFloatValues(epsB_file, &
            zone, 1, 0, nBac, iIndx(1:nBac,1))
         mat_IO(1:nBac,1) = epsB
         ioErr = tecZoneVarWriteFloatValues(epsB_file, &
            zone, 2, 0, nBac, mat_IO(1:nBac,1))            
         if (flushData) ioErr = tecFileWriterFlush(epsB_file, 0, zonesToRetain) ! flush data in memory to disk
  
         ! Write the epsG vector
         iMax = nGrz
         jMax = 1
         kMax = 1
         ioErr = tecZoneCreateIJK(epsG_file, 'epsG'//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format
         ioErr = tecZoneVarWriteFloatValues(epsG_file, &
            zone, 1, 0, nGrz, iIndx(1:nGrz,1))
         mat_IO(1:nGrz,1) = epsG
         ioErr = tecZoneVarWriteFloatValues(epsG_file, &
            zone, 2, 0, nGrz, mat_IO(1:nGrz,1))            
         if (flushData) ioErr = tecFileWriterFlush(epsG_file, 0, zonesToRetain) ! flush data in memory to disk
         
         ! Write the B matrix
         iMax = nBac
         jMax = nSub
         kMax = 1
         ioErr = tecZoneCreateIJK(B_file, 'B'//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format
         ioErr = tecZoneVarWriteFloatValues(B_file, &
            zone, 1, 0, nBac*nSub, iIndx(1:nBac,1:nSub))
         ioErr = tecZoneVarWriteFloatValues(B_file, &
            zone, 2, 0, nBac*nSub, jIndx(1:nBac,1:nSub))
         mat_IO(1:nBac,1:nSub) = B
         ioErr = tecZoneVarWriteFloatValues(B_file, &
            zone, 3, 0, nBac*nSub, mat_IO(1:nBac,1:nSub))            
         if (flushData) ioErr = tecFileWriterFlush(B_file, 0, zonesToRetain) ! flush data in memory to disk
         
         ! Write the G matrix
         iMax = nGrz
         jMax = nBac+nGrz
         kMax = 1
         ioErr = tecZoneCreateIJK(G_file, 'G'//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format
         ioErr = tecZoneVarWriteFloatValues(G_file, &
            zone, 1, 0, nGrz*(nBac+nGrz), iIndx(1:nGrz,1:nBac+nGrz))
         ioErr = tecZoneVarWriteFloatValues(G_file, &
            zone, 2, 0, nGrz*(nBac+nGrz), jIndx(1:nGrz,1:nBac+nGrz))
         mat_IO(1:nGrz,1:nBac+nGrz) = G
         ioErr = tecZoneVarWriteFloatValues(G_file, &
            zone, 3, 0, nGrz*(nBac+nGrz), mat_IO(1:nGrz,1:nBac+nGrz))            
         if (flushData) ioErr = tecFileWriterFlush(G_file, 0, zonesToRetain) ! flush data in memory to disk               
         return
      end subroutine saveEpsBG
   
      subroutine saveTimeSeries(t)
         ! This routine saves the time series variables.  However, it saves them in
         ! in tecplot zones that are nIOpts long.
         ! saveInitialSoluton must be called first, and just once, prior to a call to this routine
         ! Also, before every call, setNames must be called.
         real(mp), intent(in):: t ! time of simulation
         
         ! increament the counter
         nIO_cnt = nIO_cnt + 1
         ! Add the current solution to the IO variables that save the data ofr nIOpts points
         time_IO(nIO_cnt)        = t     ! time
         sub_IO(1:nSub,nIO_cnt)  = sj    ! subsrates
         chm_IO(1:nChm,nIO_cnt)  = chm   ! Chemistry
         Bac_IO(1:nBac,nIO_cnt)  = bi    ! Bacteria
         Grz_IO(1:nGrz,nIO_cnt)  = gi    ! Grazers
         sigma_IO(1,nIO_cnt)     = sumSigmaDot    ! summed entropy production rate
         sigma_IO(2,nIO_cnt)     = sumSigma       ! Integrated summed entropy production rate
         
         if (nIO_cnt < nIOpts) return ! arrays have not been filled yet.
         call saveIOarrays(nIOpts) ! save a full zone of data
         nIO_cnt = 1 ! reset the counter
         ! Copy the last point in the arrays to the first point  This insures when the 
         ! zones are plotted, there is overlap with the previous zone, so you get a continuous line in tecplot
         time_IO(1)        = time_IO(nIOpts) 
         sub_IO(1:nSub,1)  = sub_IO(1:nSub,nIOpts)
         chm_IO(1:nChm,1)  = chm_IO(1:nChm,nIOpts)
         Bac_IO(1:nBac,1)  = Bac_IO(1:nBac,nIOpts)
         Grz_IO(1:nGrz,1)  = Grz_IO(1:nGrz,nIOpts)
         sigma_IO(1:2,1)   = sigma_IO(1:2,nIOpts)   
         return
      end subroutine saveTimeSeries
   
      subroutine saveIOarrays(nIO) 
         ! This routine save the time serie array to a tecplot zone.  
         integer(8), intent(in):: nIO ! the number of time points in the array to save.  Most of the time it will be nIOpts
         ! local declarations
         integer i
         character(25) numStr
         
         write(numStr,*) time_IO(1) ! mark zone names with the start of the time in that zone
         iMax = nIO ! These are the same for all zones here
         jMax = 1
         kMax = 1
         
         ! Substrates
         ioErr = tecZoneCreateIJK(sub_file, 'Substrates '//trim(numStr)//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format, first time
         ioErr = tecZoneVarWriteFloatValues(sub_file, &
            zone, 1, 0, nIO, time_IO(1:nIO))
         do i=1,nSub
            ioErr = tecZoneVarWriteFloatValues(sub_file, &
               zone, i+1, 0, nIO, sub_IO(i,1:nIO))
         end do         
         if (flushData) ioErr = tecFileWriterFlush(sub_file, 0, zonesToRetain) ! flush data in memory to disk

         ! Chemistry
         ioErr = tecZoneCreateIJK(chm_file, 'Chemistry '//trim(numStr)//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format, first time
         ioErr = tecZoneVarWriteFloatValues(chm_file, &
            zone, 1, 0, nIO, time_IO(1:nIO))
         do i=1,nChm
            ioErr = tecZoneVarWriteFloatValues(chm_file, &
               zone, i+1, 0, nIO, chm_IO(i,1:nIO))
         end do         
         if (flushData) ioErr = tecFileWriterFlush(chm_file, 0, zonesToRetain) ! flush data in memory to disk
         
         ! Bacteria
         ioErr = tecZoneCreateIJK(bac_file, 'Bacteria '//trim(numStr)//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format, first time
         ioErr = tecZoneVarWriteFloatValues(bac_file, &
            zone, 1, 0, nIO, time_IO(1:nIO))
         do i=1,nBac
            ioErr = tecZoneVarWriteFloatValues(bac_file, &
               zone, i+1, 0, nIO, bac_IO(i,1:nIO))
         end do         
         if (flushData) ioErr = tecFileWriterFlush(bac_file, 0, zonesToRetain) ! flush data in memory to disk
         
         ! Grazers
         ioErr = tecZoneCreateIJK(grz_file, 'Grazers '//trim(numStr)//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format, first time
         ioErr = tecZoneVarWriteFloatValues(grz_file, &
            zone, 1, 0, nIO, time_IO(1:nIO))
         do i=1,nGrz
            ioErr = tecZoneVarWriteFloatValues(grz_file, &
               zone, i+1, 0, nIO, grz_IO(i,1:nIO))
         end do         
         if (flushData) ioErr = tecFileWriterFlush(grz_file, 0, zonesToRetain) ! flush data in memory to disk
         
         ! sumSigmaDot and sumSigma
         ioErr = tecZoneCreateIJK(sigma_file, 'sum SigmaDot '//trim(numStr)//C_NULL_CHAR, &
            iMax, jMax, kMax, varTypes, shareVarFromZone, &
            valueLocations, passiveVarList, &
            shareConnectivityFromZone, numFaceConnections, &
            faceNeighborMode, zone)
         ! write data to zone, in block format, first time
         ioErr = tecZoneVarWriteFloatValues(sigma_file, &
            zone, 1, 0, nIO, time_IO(1:nIO))
         ioErr = tecZoneVarWriteFloatValues(sigma_file, &
            zone, 2, 0, nIO, sigma_IO(1,1:nIO))
         ioErr = tecZoneVarWriteFloatValues(sigma_file, &
            zone, 3, 0, nIO, sigma_IO(2,1:nIO))
         if (flushData) ioErr = tecFileWriterFlush(sigma_file, 0, zonesToRetain) ! flush data in memory to disk
         return
      end subroutine saveIOarrays
   
      subroutine saveInitialSolution (t0)
         ! This routine saves the inital solution as well as the values of epsB, epsG, B and G
         ! Note, this routine should only be called once, and must be called before saveSolution
         ! It is assumed that setStateNames and rxnProperties have already been called to set global 
         ! variable values, and that epsB, epsG, B and G have already been generated.
         real(mp), intent(in):: t0  ! value of time for start of simulation (d)
         
         if (myRank /= 0) return
         nIO_cnt = 0 ! initialize the counter for time series data that are save as blocks of nIOpts in a zone
         ! The initial solution for these variables do not get formally save until  nIOpts is reached, or
         ! if the simulation ends before nIOpts is reached.         
         call saveTimeSeries(t0) ! the initial solution must be written to state variables before call to saveInitialSolution
         call saveRxnRates(t0) ! save rb, rg, sigB and sigG.  call to rxnProperties must have already occured.
         call saveEpsBG() ! Save epsB, epsG, B and G, this is only done once.
         ! close the EpsBG files
         ioErr = tecFileWriterClose(epsB_file)
         ioErr = tecFileWriterClose(epsG_file)
         ioErr = tecFileWriterClose(B_file)
         ioErr = tecFileWriterClose(G_file)         
         return 
      end subroutine saveInitialSolution
      
      subroutine saveSolution(t, nx, x, rpar, ipar)
         ! This routine is called by the BiM dummy routine solout
         real(mp), intent(in)   ::  t        ! time (d)
         integer,  intent(in)   ::  nx       ! size of state vector
         real(mp), intent(in)   ::  x(nx)    ! state vector solution at time t
         real(mp), intent(inout)::  rpar(*)  ! user past real(mp) vector
         integer,  intent(inout)::  ipar(*)  ! user past integer vector
      
         ! Set state name
         call setStateNames(x) ! Set all the state variables and insure they are >= absZero
         ! Calculate free energy of reactions, reaction rates, and thermodynamics. MPI is used here
         call rxnProperties(t)          
         ! save the solution
         call saveTimeSeries(t) ! Save the times series data.  This is done in chunks
         call saveRxnRates(t) ! save rb, rg, sigB and sigG.  
         return
      end subroutine saveSolution      
   
      subroutine cleanUp()
         ! close units where necessary and deallocate
         integer mpierr
         deallocate( subNames )
         deallocate( CHO )
         deallocate( B, G )
         deallocate( epsB, epsG )
         deallocate( aCs, bCs, aAs, bAs )
         deallocate( drG0_rbC, drG0_rbA )
         deallocate( drG_rb, drG_rg )
         deallocate( rb, rg )
         deallocate( sig_rb, sig_rg )         
         deallocate( sj, bi, gi )
         deallocate( biF, giF )
         deallocate( chm )         
         deallocate( sj_f )
         deallocate( wS )

         ! deallocate storage for tecIO
         if (myRank == 0) then
            ! Allocate space to store data temporarily for tecplot zones, but only for the root process
            deallocate( time_IO )   ! time
            deallocate( sub_IO )    ! subsrates
            deallocate( chm_IO )    ! Chemistry
            deallocate( Bac_IO )    ! Bacteria
            deallocate( Grz_IO )    ! Grazers
            deallocate( sigma_IO )  ! summed entropy production rate
            deallocate( mat_IO )
            deallocate( iIndx, jIndx )
            
            ! Close the tecIO files
            ioErr = tecFileWriterClose(sub_file)
            ioErr = tecFileWriterClose(chm_file)
            ioErr = tecFileWriterClose(Bac_file)
            ioErr = tecFileWriterClose(Grz_file)
            ioErr = tecFileWriterClose(rb_file)
            ioErr = tecFileWriterClose(rg_file)
            ioErr = tecFileWriterClose(sigma_file)
            
            ! close fortran units
            close(unit=iunit_opt)
         end if
         
         call MPI_FINALIZE(mpierr)         
         stop         
      end subroutine cleanUp            
   end module globalVars      
   
   module functions
      ! Note, functions are placed in this module so that real(mp) can be used instead of real(8), etc.
      ! It's not really necessary...
      ! However, but using a module, calling routines do not need an interface, which they would with the SIMD directive.
      use realPrec
      implicit none      
      
   contains 
      elemental real(mp) function F_Thermo (delG, ne, Tk)
         ! This function calculates the thermodynamic driver using LaRowe2012
         ! Function returns a unitless value between 0 and 1
         ! This directive allows the function to be processes as a vector function.  Since Tk is the same, uniform can be use
         ! but if Tk varies over a loop, then uniform should be removed.
         !$OMP DECLARE SIMD(F_Thermo) uniform(Tk)
         use globalVars, only: delPsi  ! The electric potential accross the cell wall (V).  Note, this is sigmodial 
                                       ! function, so F_Thermo is not zero if delG is less than it.  See LaRowe2012
         real(mp), intent(in):: delG  ! free energy of reaction (kJ/mol reaction)
         real(mp), intent(in):: ne    ! mole of electrons exchanged per mole of reaction extent
         real(mp), intent(in):: Tk    ! Temperature (K)
         ! Local declarations
         real(mp), parameter:: RkJ  = 8.3144598d-3   ! gas constant (kJ/(g-mol K) or J/(g-mmol K))
         real(mp), parameter:: F = 96485.3329 ! Faraday constant (C/mol-e)
         
         F_Thermo = 1.0_mp/( 1.0_mp + Exp( (delG/ne + F*delPsi/1000._mp)/(RkJ*Tk) ) )  ! Divide F*delPsi to get kJ/mol instead of J/mol
      end function  F_Thermo   
      
      real(mp) function stepUp(x, xm, sig)
         ! This routine generates a continuous sigmodial step up function from 0 to 1 around xm, with a gradient described by sig
         implicit none
         real(mp) x   ! value of x where function is evaluated
         real(mp) xm  ! value of x where setUp is 0.5
         real(mp) sig ! how steap the exponential function is around xm
         
         stepUp = 1._mp/( 1._mp + exp(-sig*(x - xm)) )
         return
      end function stepUp
      
      elemental real(mp) function fZero(absZero, x)
         ! This routine is used to prevent state variables from going negative or to zero, but in a continuous manner
         ! That is, while x >= 2 absZero, the function returns x, but as x decreases below 2 absZero, fZero also returns an number: absZero < fZero <= 2 absZero
         !$OMP DECLARE SIMD(fZero) uniform(absZero)
         implicit none
         real(mp), intent(in):: absZero  ! Value that x must stay above or at
         real(mp), intent(in):: x        ! value of x that could be zero or less than zero that needs to be prevented from that

         if (x > 2._mp*absZero) then
            fZero = x
            return
         end if
         if (x <= absZero) then
            fZero = absZero
            return
         end if
         ! x is between absZero and 2 absZero, so use quadratic for transition
         fZero = x**2/(4._mp*absZero) + absZero
         return
      end function fZero               
   end module functions
         
   Program Darwin_Chemostat
      use realPrec
      use globalVars
      use mpi
      use omp_lib
      use hyperBOB_module
      implicit none
      integer nameLen
      character (len=MPI_MAX_PROCESSOR_NAME) nodeName
      character(len=8) timeStr
      character(len=9) dateStr
      integer i, j, k, ifail, ifailtot
      integer mpierr, maxThreads, noThreads
      real(mp) tEnd0
      real(mp) WTime_0, WTime_f, cpuTime
      character(len=5) jc, kc, ic      
      character(len=80) flag
      character(len=1000) str1000
      type(solInfo) ODEinfo ! info on ODE success, etc.  
      real(mp), allocatable:: workVec(:) ! contains the trophic position of the each consumer (bacteria are all 1 in this model)
      real(mp):: qHill = 2.0_mp ! order of Hill number
      
      ! hyperBOB related variabes
      external entropy ! function hyperBOB calls to get objective function value
      integer npt_bobyqa ! number of interpolation conditions on each quadratic model used in BOBYQA 
      real(mp) maxTime ! not used by hyperBOB yet
      real(mp) MEP ! maximum entropy production found integrated over time (J/K)
      real(mp), allocatable:: traits(:) ! non-dimensional traits vector between 0 and 1
      real(mp), allocatable:: zeroVec(:), oneVec(:) ! vector of 0's and 1's
      integer ierr ! used by hyperBOB, but it's trivial error, so not cheched here.      
      real(mp), allocatable:: fxmat(:,:)

      ! Initialize MPI
      call MPI_INIT( mpierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpierr) ! get rank of this process in world

      if (myRank == 0) write(*,'(/a/)') 'Darwin Chemostat Version: '//version
      call initialize() ! read in parameters, allocate space, calculate variables.
      allocate( workVec(max(nBac,nGrz)) ) ! vector to store effective trophic position of the consumers
      tf = t0 + tDays ! set the end of the simulation    
      tf_mep = t0_mep + tDays_mep ! the end of the mep optimization
      if (t0_mep < t0) then
         t0_mep = t0
         if(myRank == 0) then
            write(*,'(a)') 'WARNING:: t0_mep < t0, optimization will start at t0, not t0_mep'
         end if
      end if      
      if (tf_mep > tf) then
         tf_mep = tf
         if(myRank == 0) then
            write(*,'(a)') 'WARNING:: tf_mep > tf, optimization will end at tf, not tf_mep'
         end if
      end if      
      maxODEtime = (tf-t0)/minCompFac ! maximum allowed CPU wall time (days) allowed to solve the PDE over t0 to tf.            
      if (myRank == 0) then 
         write(*,'(/a)') 'Parameters being used for simulation:'
         write(*,nml=dimens)
         write(*,nml=parameters)
         write(*,nml=hyperBOBparams)
         write(*,nml=tecIOparams)
         write(*,nml=BiMparams)
         write(*,'(/,a,f0.2,a)') 'Maximum ODE solution time set to:  ', maxODEtime*24.0*60., ' (min)'
         write(*,'(a,f0.2,a,/)') 'Maximum possible optmization time: ', maxODEtime*real(maxFun), ' (days)'
      end if

      ! Get the number of processes running and the maxium number of threads
      call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)
      maxThreads = OMP_GET_MAX_THREADS()
      call omp_set_num_threads(min(maxThreads, ompThreads)) ! Don't use more threads than available
      !$OMP parallel
      noThreads = OMP_GET_NUM_THREADS() ! Outside of omp parallel regions, this returns 1
      call initThreads() ! Allocate space for threadprivate variables while in parallel region
      !$OMP end parallel
      if (myRank == 0) write(*,'(a,i4,a)') 'Running program with ',noProc, ' processes'     
      call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      call MPI_GET_PROCESSOR_NAME(nodeName, nameLen, mpiErr)
      do i=0,noProc-1
         if (i == myRank) then
            write(*,'(a,i4,2a,2(a,i4))') ' Process ', myRank, ' on node "', trim(nodeName), '" has ', noThreads, ' OpenMP threads out of: ', maxThreads
            !write(*,'(a,i12,/)')  ' The thread stack size is set to: ', KMP_GET_STACKSIZE_S()   
         end if
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      end do     
      WTime_0 = MPI_Wtime() ! Start the clock
      
      call time(timeStr); call date(dateStr) 
      WTime_it0 = MPI_Wtime() ! Get the clock start time
      fcnCalls = 0
      ODEfailed = 0; ODEtimedOut = 0
      ! set npt_bobyqa, which is Number of interpolation points between [N+2,(N+1)(N+2)/2]. best if < 2*N+1    
      ! This could be made a parameter, but the suggested value below seems to work pretty well in general
      npt_bobyqa = 2*nTraits + 1      
      ! allocate space
      allocate( traits(nTraits) )
      allocate( zeroVec(nTraits), oneVec(ntraits) )
      zeroVec = 0._mp; oneVec = 1._mp ! Used for upper and lower bounds on normalized trait vector
      call setBounds() ! Sets lower and upper bounds on traits
      call hyperBOB_initialize(nTraits) ! setup space for hyperBOB      
      
      if (optimize) then
         onlyCheckSoln = .false.
         if (myRank == 0) then
            write(*,'(/,a)') 'Begining hyperBOB optimization at '//timeStr//', '//dateStr//' ...'
         end if                  
         call hyperBOB(nTraits, npt_bobyqa, traits, zeroVec, oneVec, rhoBeg, rhoEnd, iprint, maxFun, maxTime, entropy, seed_bobyqa, MEP, ierr)
         MEP = -MEP ! hyperBOB finds the minimum, so convert to max here
         WTime_f = MPI_Wtime() ! get the clock time
         cpuTime = real(WTime_f-WTime_it0)/60.0/60. ! hrs              
         do i=0,noProc-1
            if (i == myRank) then
               flag = ''
               if (fcnCalls > maxfun) flag = ' **** maxfun calls exceeded ****'
               write(*,'(a,i0,a,2(i0,'', ''),i0,a,g0.7,a)') &
                  '-> Process ', myRank, ' finished (fCalls, F, TO: ', fcnCalls, ODEfailed, ODEtimedOut,'). fVal = ', -fVal_hyperBOB(myRank+1), trim(flag) 
            end if
            call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
         end do  
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr) 
                 
         call MPI_ALLREDUCE(fcnCalls   , fcnCallsTot   , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum fcnCalls for all processes.  
         call MPI_ALLREDUCE(ODEfailed  , ODEfailedTot  , 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum PDEfailed for all processes.           
         call MPI_ALLREDUCE(ODEtimedOut, ODEtimedOutTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum PDEtimedOut for all processes.           
                  
         call time(timeStr); call date(dateStr)       
         if (myRank == 0) then
            ! Search fVal_hyperBOB for the minimum location in MEP   
            write(*,'(/,a,f0.2,a)') 'hyperBOB optimization finished at at '//timeStr//', '//dateStr//' taking ', cpuTime, ' hrs'
            write(*,'(3(a,i0))') '   Total number of ODE integrations: ', fcnCallsTot, ' ODEfailed: ', ODEfailedTot, ' ODEtimedOut: ', ODEtimedOutTot
            write(*,'(a,i0,a,f0.0)') '   Process ', minloc(fVal_hyperBOB)-1, ' found best solution: ', MEP
            write(*,'(a)') ' '
            
            ! Output all the solutions but only the obj value and traits.
            ! sort fVal_hyperBOB and xMat_hyperBOB and store all optimum solutions in *_allOpt.dat
            allocate( fxmat(noProc,nTraits+1) )
            fxmat(1:noProc,1) = fVal_hyperBOB(1:noProc)
            fxmat(1:noProc,2:nTraits+1) = transpose(xMat_hyperBOB(1:nTraits,1:noProc))
            call dsortArray (noProc, nTraits+1, fxmat, 1)
            open(newunit=iunit_allOpt, file=trim(basefilename)//'_allOpt.dat' ,status='unknown')         
            ! tecplot header 
            ! Note, it's not worth trying to make this file easy to view, so just wrap the data as Tecplot can read that.
            ! Also, the best solution is stored in the tecplot binaries, but that's not needed so much with the optimized version now.
            write(iunit_allOpt,'(a)') 'Variables = "Index" "obj"'
            ! Bacteria
            do i=1,nBac
               write(ic,'(i0)') i
               str1000 = '"epsB('//trim(ic)//')"'               
               do j=1,nSub
                  write(jc,'(i0)') j
                  str1000 = trim(str1000)//' "B('//trim(ic)//','//trim(jc)//')"'
               end do
               write(iunit_allOpt,'(a)') trim(str1000)
            end do
            ! Consumers
            do i=1,nGrz
               write(ic,'(i0)') i
               str1000 = '"epsG('//trim(ic)//')"'               
               do j=1,nBac+nGrz
                  write(jc,'(i0)') j
                  str1000 = trim(str1000)//' "G('//trim(ic)//','//trim(jc)//')"'
               end do
               write(iunit_allOpt,'(a)') trim(str1000)
            end do        
            write(iunit_allOpt,'(a)') 'Zone T="Traits: all solutions"'
            
            ! open file for output of stats iunit_stats
            open(newunit=iunit_stats, file=trim(basefilename)//'_stats.dat' ,status='unknown')
            open(newunit=iunit_rStats, file=trim(basefilename)//'_rStats.dat' ,status='unknown')
            write(iunit_stats, '(a)') 'Variables = "Index" "obj-Int" "obj-full"'
            write(iunit_rStats,'(a)') 'Variables = "Index" "obj-Int" "obj-full"'
            ! bacterial substratre use diversity index
            str1000 = ' "SD-B(1)"'
            do i=2,nBac
               write(ic,'(i0)') i
               str1000 = trim(str1000)//' "SD-B('//trim(ic)//')"'
            end do
            write(iunit_stats,'(a)') trim(str1000)
            write(iunit_rStats,'(a)') trim(str1000)
            ! consumer equivelent trophic position
            str1000 = ' "TP-C(1)"'
            do i=2,nGrz
               write(ic,'(i0)') i
               str1000 = trim(str1000)//' "TP-C('//trim(ic)//')"'
            end do
            write(iunit_stats,'(a)') trim(str1000)
            write(iunit_stats,'(a)') 'Zone T="Stats: all solutions"'
            write(iunit_rStats,'(a)') trim(str1000)
            write(iunit_rStats,'(a)') 'Zone T="Rxn-based Stats at tf: all solutions"'

            do i=1,noProc
               ! save all the trait values
               call mapTraits(nTraits, fxmat(i,2:nTraits+1)) ! This does all the mapping and conversion to dimensional form
               write(iunit_allOpt,'(i0,1x,g0.7)') i, -fxmat(i,1) ! remove the negative sign from the objective
               do j=1,nBac
                  write(iunit_allOpt,'(g0.7,*(1x,g0.7))') epsB(j), B(j,1:nSub)
               end do               
               do j=1,nGrz
                   write(iunit_allOpt,'(g0.7,*(1x,g0.7))') epsG(j), G(j,1:nBac+nGrz)
               end do  
               ! run the model over the whole time domain to get objective fcn for that
               biF = 0._mp; giF = 0._mp ! set flags to 0
               onlyCheckSoln = .true.
               call integrateState(t0, tf, .false., ODEinfo) ! don't check for errors. sumSigma, rb and rg are set here
               write(iunit_stats,'(i0,2(1x,g0.7))') i, -fxmat(i,1), sumSigma ! save process and objective fuctions
               write(iunit_rStats,'(i0,2(1x,g0.7))') i, -fxmat(i,1), sumSigma ! save process and objective fuctions
               ! calculate the effective diversity of substrates used by each bacteria
               ! based on a Hill number with q = 2 (see https://en.wikipedia.org/wiki/Diversity_index)         
               ! do the same, but use rb
               workVec = 0.0_mp
               do j=1,nBac
                  do k=1,nSub
                     workVec(j) = workVec(j) + B(j,k)**qHill
                  end do
                  workVec(j) = workVec(j)**(1._mp/(1._mp - qHill))
               end do
               workVec(1:nBac) = biF(1:nBac)*workVec(1:nBac) ! remove values if bi is never significant
               write(iunit_stats,'(*(g0.7,1x))') workVec(1:nBac)
               ! do the same thing, but based on the normalized version of rb at tf
               workVec = 0.0_mp
               do j=1,nBac
                  rb(j,1:nSub) = rb(j,1:nSub)/sum(rb(j,1:nSub)) ! normalize the row of rb
                  do k=1,nSub
                     workVec(j) = workVec(j) + rb(j,k)**qHill
                  end do
                  workVec(j) = workVec(j)**(1._mp/(1._mp - qHill))
               end do                
               where(bi(1:nBac) <= bi_ini) workVec(1:nBac) = 0._mp ! if bi at tf is <= bi_ini, set their SD to zero
               write(iunit_rStats,'(*(g0.7,1x))') workVec(1:nBac)
               
               ! get trophic position of the consumers
               call trophicPos(G, workVec)         
               workVec(1:nGrz) = giF(1:nGrz)*workVec(1:nGrz) ! remove values if gi is never significant
               write(iunit_stats,'(*(g0.7,1x))') workVec(1:nGrz)    
               ! so the same, but for normalized rg
               do j=1,nGrz
                  rg(j,1:nBac+nGrz) = rg(j,1:nBac+nGrz)/sum(rg(j,1:nBac+nGrz))
               end do               
               call trophicPos(rg, workVec)
               where(gi(1:nGrz) <= gi_ini) workVec(1:nGrz) = 0._mp ! if consumer is <= gi_ini, set TP to zero
               write(iunit_rStats,'(*(g0.7,1x))') workVec(1:nGrz)               
            end do
            close(iunit_allOpt)
            close(iunit_stats)
            close(iunit_rStats)
            deallocate( fxmat )
         end if
      else
         ! No optimization, just read the *_opt.dat file, only master process does this
         if (myRank == 0) then 
            write(*,'(/,a)') 'Solve ODE once with trait values from '//trim(basefilename)//'_opt.dat'//', No optimization'
            read(iunit_opt,'(3/)',iostat=ioerr) ! This skips the header in the *_opt.dat file
            if (ioerr /= 0) then
               write(*,'(a)') 'File appears to be empty, Stopping run.'
               call cleanUp()
               deallocate( traits )               
               stop
            end if
            do i=1,nBac
               read(iunit_opt,*) epsB(i), B(i,1:nSub)
            end do 
            do i=1,nGrz
               read(iunit_opt,*) epsG(i), G(i,1:nBac+nGrz)
            end do                     
            rewind(iunit_opt) ! This file is written to below, so rewind it to the top.
            ! noramlize the B and G trait matracies. These get normalized again by substrate or prey densities, but this makes them easiler to visualize.
            ! and they need to be normalized for the stats calculations (substrate diversity and trophic position)
            do i=1,nBac
               B(i,1:nSub) = B(i,1:nSub)/sum(B(i,1:nSub))
            end do
            do i=1,nGrz
               G(i,1:nBac+nGrz) = G(i,1:nBac+nGrz)/sum(G(i,1:nBac+nGrz))
            end do                  
         end if               
      end if
      
      ! Optimium traits have been found, or read from the opt file if not using optimization
      if (myRank == 0) then
         if (optimize) call mapTraits(nTraits, traits) ! map traits bact to dimensional values and associated variables
         write(*,'(/,a)') 'Integrating ODE with best (or read) trait values'
         onlyCheckSoln = .false.
         call integrateState(t0, tf_mep, .false., ODEinfo) ! this gets MEP over t0_mep to tf_mep
         MEP = sumSigma
         cpuTime = MPI_Wtime() ! get the clock time   
         biF = 0._mp; giF = 0._mp ! set flags to 0   
         onlyCheckSoln = .false.
         call integrateState(t0, tf, .true., ODEinfo) ! this stores solutoin AND gets MEP over t0 to tf
         if (nIO_cnt > 1) call saveIOarrays(nIO_cnt) ! This saves the timeseries data to the last zone
         if (ODEinfo%fail /= 0) then
            write(*,'(a)') 'ERROR: integration of ODE with best (or read) solution FAILED!!'
            write(*,'(a,i0)') '   ODEinfo%fail = ', ODEinfo%fail
            write(*,'(a,i0)') '   ODEinfo%t0   = ', ODEinfo%t0 
            write(*,'(a,i0)') '   ODEinfo%idid = ', ODEinfo%idid            
         end if         
         ! Save best parameters only (this also re-saves the trait values if no optimization was done)
         write(iunit_opt,'(a,g0.7,1x,g0.7)') 'Maximum integrated entropy production over short and long intervals (J/K): ', MEP, sumSigma
         write(iunit_opt,'(a)') 'Optimum parameters save as follows:' 
         write(iunit_opt,'(a)') '  epsB(i), B(i,1:nSub): one row for each nBac'
         write(iunit_opt,'(a)') '  epsG(i), G(i,1:nBac+nGrz): one row for each nGrz'
         do i=1,nBac
            write(iunit_opt,'(g0.7,*(1x,g0.7))') epsB(i), B(i,1:nSub)
         end do 
         do i=1,nGrz
            write(iunit_opt,'(g0.7,*(1x,g0.7))') epsG(i), G(i,1:nBac+nGrz)
         end do 
         ! calculate the effective diversity of substrates used by each bacteria
         ! based on a Hill number with q = 1 (see https://en.wikipedia.org/wiki/Diversity_index)         
         workVec = 0.0_mp
         do i=1,nBac
            do j=1,nSub
               workVec(i) = workVec(i) + B(i,j)**qHill
            end do
            workVec(i) = workVec(i)**(1._mp/(1._mp - qHill))
         end do
         workVec(1:nBac) = biF(1:nBac)*workVec(1:nBac) ! remove values if bi is never significant
         write(iunit_opt,'(/,a)') 'Substrate diversity for each bacteria'
         write(iunit_opt,'(*(g0.7,1x))') workVec(1:nBac)
         ! do the same thing, but base it on rb at tf
         workVec = 0.0_mp
         do i=1,nBac
            rb(i,1:nSub) = rb(i,1:nSub)/sum(rb(i,1:nSub)) ! normalize the row of rb            
            do j=1,nSub
               workVec(i) = workVec(i) + rb(i,j)**qHill
            end do
            workVec(i) = workVec(i)**(1._mp/(1._mp - qHill))
         end do
         where(bi(1:nBac) <= bi_ini) workVec(1:nBac) = 0._mp ! if bi at tf is <= bi_ini, set their SD to zero
         write(iunit_opt,'(/,a)') 'Substrate diversity for each bacteria based on rb at tf'
         write(iunit_opt,'(*(g0.7,1x))') workVec(1:nBac)         
         ! get trophic position of the consumers
         call trophicPos(G, workVec)
         workVec(1:nGrz) = giF(1:nGrz)*workVec(1:nGrz) ! remove values if gi is never significant
         write(iunit_opt,'(/,a)') 'Trophic position of the consumers'
         write(iunit_opt,'(*(g0.7,1x))') workVec(1:nGrz)
         ! get trophic position of the consumers based on rg at tf
         do i=1,nGrz
            rg(i,1:nBac+nGrz) = rg(i,1:nBac+nGrz)/sum(rg(i,1:nBac+nGrz))
         end do               
         call trophicPos(rg, workVec)         
         where(gi(1:nGrz) <= gi_ini) workVec(1:nGrz) = 0._mp ! If gi at tf is <= gi_ini, set their TP to zero
         write(iunit_opt,'(/,a)') 'Trophic position of the consumers based on rg at tf'
         write(iunit_opt,'(*(g0.7,1x))') workVec(1:nGrz)
         cpuTime = MPI_Wtime() - cpuTime
         write(*,'(a,g0.7)') '   Finished: CPU time (min): ', cpuTime/60._mp         
      end if
            
      ! Deallocate global vars, close tecIO files, and stop
      deallocate( workVec )
      call hyperBOB_cleanUp()      
      call cleanUp
      deallocate( traits, zeroVec, oneVec )
   end program Darwin_Chemostat
   
   subroutine setBounds()
      ! This routine sets the lower and upper bounds on the traits vector, which is used normalize the triat values
      ! between 0 and 1 for call to hyperBOB. Most of them are just 0 and 1 anyway, but this allows some changes to that.
      use realPrec
      use globalVars
      implicit none
   
      traitsL = 0._mp; traitsU = 1._mp ! most will already be in the range [0,1]
      traitsL(1:nBac) = absZero ! keep lower bound on epsB above 0.
      traitsL(nBac+nBac*nSub+1:nBac+nBac*nSub+nGrz) = absZero ! keep lower bound on epsG above 0.
      return
   end subroutine setBounds
   
   subroutine mapTraits(nT, traitsND)
      ! This routine maps dimensionless traits back to the dimensionalized values held in globalVars module
      use realPrec
      use globalVars
      implicit none
      integer,  intent(in)               :: nT     ! number of traits
      real(mp), intent(in), dimension(nT):: traitsND ! non-dimensional trait values between 0 and 1
      ! local declarations
      integer i
      real(mp), dimension(nT):: traits ! demensionalized trait values
      real(mp) demom
      
      traits = traitsND*(traitsU - traitsL) + traitsL ! dimensionalized trait values
      epsB(1:nBac) = traits(1:nBac)
      B(1:nBac,1:nSub) = reshape(traits(nBac+1:nBac + nBac*nSub), [nBac,nSub])
      epsG(1:nGrz) = traits(nBac + nBac*nSub + 1:nBac + nBac*nSub + nGrz)
      G(1:nGrz,1:nBac+nGrz) = reshape(traits(nBac + nBac*nSub + nGrz + 1:nBac + nBac*nSub + nGrz + nGrz*(nBac+nGrz)), [nGrz,nBac+nGrz])
      ! noramlize the B and G trait matracies. These get normalized again by substrate or prey densities, but this makes them easiler to visualize.
      do i=1,nBac
         B(i,1:nSub) = B(i,1:nSub)/sum(B(i,1:nSub))
      end do
      do i=1,nGrz
         G(i,1:nBac+nGrz) = G(i,1:nBac+nGrz)/sum(G(i,1:nBac+nGrz))
      end do      
      return
   end subroutine mapTraits
   
   subroutine entropy (n, x, f)
      ! this routine is called by hyperBOB to find the function value, in this case entropy production over simulation time t0 to tf
      use realPrec
      use globalVars
      use mpi
      use kind_module
      implicit none
      ! Dummy variables
      ! ****** NOTE, the function must be defined EXACTLY as shown by the interface for bobyqa in hypterBOB .  *******
      ! that is, don't instead declare as "real(wp) x(n)", etc as that will cause problems      
      integer,intent(in)               :: n ! number of control variables
      real(wp),dimension(:),intent(in) :: x ! control variables (i.e., traits)
      real(wp),intent(out)             :: f ! objective function, integrated entropy production
      ! Local declarations
      integer i, j    
      type(solInfo) ODEinfo ! info on ODE success, etc.
      real(mp) WTime_f, cpuTime, WTime_s, intTime, eProd(3)
      character(len=8) timeStr
      character(len=9) dateStr
      
      fcnCalls = fcnCalls + 1
      call mapTraits(n, x)  ! set epsB, B, epsG and G  
      WTime_s = MPI_Wtime() ! get the clock time
      call integrateState(t0, tf_mep, .false., ODEinfo) ! only integrate up to tf_mep.
      WTime_f = MPI_Wtime() ! get the clock time
      f = -sumSigma ! sumSigma was set when integrateState exits. hyperBOB finds minimum, so change sign.
      if (ODEinfo%fail /= 0) then
         f = 0._mp ! BOBYQA does not handle failures, so just set to 0.
         if (ODEinfo%fail == 1) ODEfailed   = ODEfailed   + 1
         if (ODEinfo%fail == 2) ODEtimedOut = ODEtimedOut + 1
      end if
      if (-f > fcnMax) fcnMax = -f
      cpuTime = real(WTime_f-WTime_it0)/60.0/60. ! hrs
      intTime = real(WTime_f-WTime_s)/60.0 ! min
      aveODEtime = aveODEtime + intTime ! determine the average PDE integration time over fcnUpdate steps
      if (mod(fcnCalls, fcnUpdate) == 0._mp) then
         call time(timeStr); call date(dateStr)         
         write(*,'(a,i0,a,f0.1,3(a,i0),a,f0.0,a,f0.2,a)')  'Process ', myRank, ' after ', cpuTime, ' hrs, ODE integtrated ', &
            fcnCalls, ' times (',ODEfailed, ' fail, ', ODEtimedOut, ' TO), fcnMax = ', fcnMax, &
            ' (ODE time: ', aveODEtime/real(fcnUpdate),' min; '//timeStr//', '//dateStr//')'
         aveODEtime = 0._mp         
      end if      
      return
   end subroutine entropy     
   
   subroutine integrateState(t0_ode, tf_ode, strSoln, ODEinfo)
      ! This routine integrates the state equations
      ! Most parameters are passed via module parameters.
      ! Since this may be later developed to run with OpenMP, all variables should be local 
      ! or allocated; that is, keep thread safe with respect to cntl and any variable that is assigned.
      use mpi
      use realPrec
      use globalVars
      implicit none
      real(mp),      intent(in) :: t0_ode ! Start time of integration (d)
      real(mp),      intent(in) :: tf_ode ! Stop time of integration (d)
      logical,       intent(in) :: strSoln ! if set to true, the solution is stored.
      type(solInfo), intent(out):: ODEinfo 
      ! The solInfo derived type has three variables
      ! If the integration fails or times out, ifail is set to non zero value
      !  ifail = 0 ODE integrated without error and within the max time allowed
      !  ifail = 1 ODE generated an error
      !  ifail = 2 ODE took longer than it's maximum time allowed.
      !  idid: value of idid from BACOLI
      !  t0: value of t0 (where the solution got to)

      ! local declarations
      integer i, j, mpierr
      real(mp) ODE_ti ! current wall time at iteration
      logical repeat
      integer attempts, iout
      character(len=200) string
      character(len=8) timeStr
      character(len=9) dateStr
      
      ! Local declarations needed for BiM
      real(mp) x(nState), tEnd          
      integer lenw, leniw
      integer  ierflg, ijac, mljac, mujac, idid, ierr, ipar(1)
      real(mp) rtolL, atolL, h_BiM, rpar(1), t0L
      real(mp), allocatable:: wk(:)
      integer, allocatable:: iwk(:)
      External ODEs, ompJAC

      ODEinfo%fail = 0 
      ODEinfo%t0 = t0_ode
      ODEinfo%idid = 0
      if(tf_ode <= t0_ode) return ! nothing to integrate.
      
      ! Space needed by BIM
      lenw = 14+10+8*nstate+4*10*nstate+2*nstate**2
      leniw = nstate+37
      allocate (  wk(lenw) )
      allocate ( iwk(leniw) )
      
      repeat = .True.
      attempts = 1
      rtolL = rtol ! relative and absolute tolerances in parameter module
      atolL = atol
      
      !************* Begin ODE Solution from t0_ode to tf_ode ***********************
      ODE_t0 = MPI_Wtime() ! start clock for this solution
      t0L = t0_ode ! on output from BiM, t0_ode is where the solution got to, so make this local variable   
      ipar(1) = 1 ! If a solution is not being saved, then assume entropy production is integrated over t0_mep to tf_mep (see ODE routine)
      if (strSoln) ipar(1) = 0 ! if the solution is being saved, then calculated entropy production starting at t0
      call initializeState(x) ! initialize state variables
      if (strSoln) then ! save initial conditions         
         call setStateNames(x) ! Set all the state variables and insure they are >= absZero         
         call rxnProperties(t0_ode) ! Calculate free energy of reactions, reaction rates, and thermodynamics.                    
         call saveInitialSolution(t0_ode) ! Save iniatial solution in tecplot files.
      end if
      Do While (repeat)
         repeat = .False. 
         h_BiM = 0. ! use default initial BiM stepsize
         ijac = 1 ! BiM: 0=calculate jacobian numerically, otherwise use analytical. Note, the supplied "analytical" jacobian uses numerical approach here.
         iout = 0
         if (strSoln .or. onlyCheckSoln) iout = 1 ! BiM: set to 1 to call solout.  BiM determines how frequently to save solution points
         mljac = nstate; mujac = nstate !BiM: banding aspect of jacobian? just set to nstate.
         iwk(1:8) = 0
         wk(1:14) = 0.0d0
         iwk(1) = maxstep_BiM ! default iterations 100000 used if set to zero
         wk(1) = epsilon(wk) ! set true precision
         ! note, after email exchanges with Cecilia Magherini, hmax actual refers to internal mesh spacing, that depends on the order of the method
         ! At maximum orde (12), there can be 10 interval, so you need to divide hmax by 10 if an output point needs to be < hmax
         wk(2) = hmax_BiM/10._dp ! largest integration step size (use to handle high frequency drivers)
         idid = 0
         CALL BiM (nstate, ODEs, t0L, tf_ode, x, h_BiM, rtolL, atolL, ompJAC, ijac, mljac, mujac, &
                     wk, lenw, iwk, leniw, rpar, ipar, iout, idid) ! ncntl declared in parameters module
         ODE_ti = MPI_Wtime() ! get the clock time
         ODEinfo%t0 = t0L
         ODEinfo%idid = idid
         call setStateNames(x)  ! this is mostly used to set sumSigma at t0L (=tf_ode if no failures)
         call rxnProperties(t0L) ! This sets r
         ! First check if maximum ODE time has been exceeded.
         if ((ODE_ti-ODE_t0)/86400._mp > maxODEtime) then
            ! Maximum solution time exceeded. Note, BiM does not handle ierr in the ODE subrouinte correctly, so it can return idid = 0 when it should have been -6.
            ODEinfo%fail = 2
            call time(timeStr); call date(dateStr)         
            write(string,'(a,i0,a,f0.2,a)') 'Process: ', myRank, ' ODE TIME EXCEEDED (t0 = ',  t0L, ') ('//timeStr//', '//dateStr//')'    
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)            
            return ! exits this ODE solution loop
         end if
         ! Check for integration errors
         If (idid < 0) Then             
            ! Save the error info
            call time(timeStr); call date(dateStr)         
            write(string, '(a,i4,a,f0.2,2(a,i0),a)') 'Error in BiM::idid = ', idid, ' at time: ', t0L, ' Attempts: ', attempts, &
               ' for Process: ', myRank,' ('//timeStr//', '//dateStr//')'
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)            
            rtolL = 10.*rtolL
            atolL = 10.*atolL
            repeat = .True.
            if (attempts > maxattempts) then
               write(string, '(a,i4,a,f0.2,a,i0,a)') '**MAX** Attempts in BiM::idid = ', idid, ' at time: ', t0L, ' for Process: ', &
                  myRank,' ('//timeStr//', '//dateStr//')'
               call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)               
               ODEinfo%fail = 1
               return
            end if
            attempts = attempts + 1
         End If
      End Do      
      deallocate (  wk )
      deallocate ( iwk )
      return
   end subroutine integrateState

   subroutine ODEs(nx, t, x ,xdot, ierr, rpar, ipar)
      ! This routine contains the state equations
      ! for which a solution is being sought.
      use realPrec
      use globalVars
      use thermoData
      use mpi
      implicit none
      integer,  intent(in) :: nx       ! Number of ODEs
      real(mp), intent(in) :: t        ! Time for evalution or RHS of ODEs (d) 
      real(mp), intent(in) :: x(nx)    ! Value of state variables at time t
      real(mp), intent(out):: xdot(nx) ! LHS of ODE equations (time derivative) 
      integer,  intent(out):: ierr     ! Return 0 if no errors occured, otherwise set to /= 0 (I assume BiM just terminates then)
      real(mp), intent(in) :: rpar(*)  ! user real vector that is passed from call to BiM
      integer,  intent(in) :: ipar(*)  ! user integer vector that is passed from call to BiM

      ! local declarations
      integer i, j
      real(mp) h_o2, h_co2 ! Henrys constant for O2 and CO2 in mmol/m3/mbar
      real(mp) ODE_ti ! current wall time at iteration
      
      ODE_ti = MPI_Wtime() ! start clock for this solution      
      if ((ODE_ti-ODE_t0)/86400._mp > maxODEtime) then      
         ierr = 1 ! any value other than 0 works.
         return
      end if      
      
      call setStateNames(x) ! Set all the state variables and insure they are >= absZero
      ! Calculate free energy of reactions, reaction rates, and thermodynamics. MPI is used here
      call rxnProperties(t)    
      h_o2  = solO2 (T_K, is, pH)/1.01325_mp/1000._mp ! solO2  is mmol/m3/atm, so convert to mmol/m3/mbar  See thermoData module
      h_co2 = co2g2l(T_K, is, pH)/1.01325_mp/1000._mp ! co2g2l is mmol/m3/atm, so convert to mmol/m3/mbar  See thermoData module
      ! Define ODE equations
      ! Substrates
      do j=1,nSub
         xdot(j) = F_L/V_L*(sj_f(j) - sj(j)) - sum(rb(1:nBac,j))/CHO(j,1)
      end do
      
      ! The nChm chemical species
      ! o2
      xdot(io2) = F_L/V_L*(o2_f - o2) + kL_o2*area_GL/V_L*(po2*h_o2 - o2) &
                - dot_product( matmul(rb(1:nBac,1:nSub),aCs(1:nSub)), (1._mp - epsB(1:nBac)) ) &
                - aCg*sum( matmul( (1._mp - epsG(1:nGrz)), rg(1:nGrz,1:nBac+nGrz) ) )
      
      ! po2
      xdot(ipo2) = F_G/V_G*(po2_f - po2) - kL_o2*area_GL/V_L*Rbar*T_K*(po2*h_o2 - o2)
      
      ! h2co3
      xdot(ih2co3) = F_L/V_L*(h2co3_f - h2co3) + kL_co2*area_GL/V_L*(pco2*h_co2 - freeCO2(T_K, is, pH)*h2co3) &
                   + dot_product( matmul(rb(1:nBac,1:nSub),aAs(1:nSub)), epsB(1:nBac) ) &
                   + sum( matmul( (1._mp - epsB(1:nBac)), rb(1:nBac,1:nSub) ) )         &
                   + sum( matmul( (1._mp - epsG(1:nGrz)), rg(1:nGrz,1:nBac+nGrz) ) )
      
      ! pco2
      xdot(ipco2) =  F_G/V_G*(pco2_f - pco2) - kL_co2*area_GL/V_G*Rbar*T_K*(pco2*h_co2 - freeCO2(T_K, is, pH)*h2co3) 
      
      ! h3po4
      xdot(ih3po4) = F_L/V_L*(h3po4_f - h3po4) &
                   - zet_Bio/alf_Bio*dot_product( matmul(rb(1:nBac,1:nSub),(1._mp - aAs(1:nSub))), epsB(1:nBac) ) &
                   + zet_Bio/alf_Bio*sum( matmul( (1._mp - epsG(1:nGrz)), rg(1:nGrz,1:nBac+nGrz) ) )
      
      ! bac
      do i=1,nBac
         xdot(nSub+nChm+i) = F_L/V_L*(bi_f - bi(i)) &
                           + epsB(i)/alf_Bio*dot_product( (1._mp - aAs(1:nSub)),rb(i,1:nSub) ) &
                           - sum( rg(1:nGrz,i) )/alf_Bio
      end do
      
      ! grz
      do i=1,nGrz
         xdot(nsub+nChm+nBac+i) = F_L/V_L*(gi_f - gi(i)) &
                                + epsG(i)/alf_Bio*sum( rg(i,1:nBac+nGrz) ) &
                                - sum( rg(1:nGrz,i+nBac) )/alf_Bio
      end do
      
      ! integrated entropy
      if (t < t0_mep .and. ipar(1) == 1) sumSigmaDot = 0._mp ! don't add EP if before t0_mep
      xdot(nsub+nChm+nBac+nGrz+1) = sumSigmaDot
      ierr = 0   
      return      
   end subroutine ODEs   
   
   subroutine dummyJAC(n,t,x,jacmat,ldjac,ierr,rpar,ipar)
      ! This a dummy jacobian routine
      use realPrec
      implicit none
      integer n, ldjac, ierr, ipar(*)
      real(mp) t,x(n), jacmat(ldjac,n), rpar(*)   
      return
   end subroutine dummyJAC

   subroutine ompJAC(n,t,x,jacmat,ldjac,ierr,rpar,ipar)
      ! This calculates a numberical jacobian using the same method in BiM code
      ! but uses OpenMP to speed things up.
      use realPrec
      implicit none
      integer,  intent(in)    :: n              ! Size of system
      real(mp), intent(in)    :: t              ! time
      real(mp), intent(inout) :: x(n)           ! state variables. Not changed, but do need to temp alter.
      real(mp), intent(out)   :: jacmat(ldjac,n)! Jacobian output
      integer,  intent(in)    :: ldjac          ! Size of jacmat
      integer,  intent(out)   :: ierr           ! if errors occur
      real(mp), intent(inout) :: rpar(*)        ! User defined vector
      integer,  intent(inout) :: ipar(*)        ! User defined vector
      ! Local declarations
      integer i, j
      real(mp) xsave, delt, uround 
      real(mp) xdot(n), xdot0(n)
  
      ! Values +1 is still 1. BiM uses 1.0d-16, but could use epsilon here.
      uround = 1.0d-16 
      
      ! get the value of xdot at current x
      call ODEs(n, t, x ,xdot0, ierr, rpar, ipar)
      continue
      ! Begin parallel calculation of jacobian matrix
      !$OMP PARALLEL DO firstprivate(x) private(i,j,xsave,delt,xdot) schedule(static)      
      do i=1,n
         xsave = x(i)
         delt = sqrt(uround*max(1.D-5,abs(xsave)))
         x(i) = xsave+delt
         call ODEs(n, t, x ,xdot, ierr, rpar, ipar)
         do j=1,n
            jacmat(j,i)=(xdot(j)-xdot0(j))/delt
         end do
         x(i)=xsave
      end do   
      !$OMP END PARALLEL DO      
      ierr = 0
      return
   end subroutine ompJAC
   
   subroutine rxnProperties(t)
      ! This version (1.2) using weighting based on B or G and relative abundance of substrate or prey
      ! This routine calculate Gibbs free energy of reactions, reaction rates and entropy production at time t
      ! The calculations are stored in global arrays defined in globalVars module
      ! **** NOTE, routine setStateNames should be called before calling this routine to define the named
      ! **** state variables.
      use mpi
      use realPrec
      use globalVars
      use functions
      use thermoData
      implicit none
      real(mp), intent(in) :: t           ! current time (days)
      ! local declarations
      integer i, j, mpierr
      real(mp) k_eps4, drG_A, drG_C
      real(mp), parameter:: ln106 = -6._dp*log(10._dp) ! this is the natural logorithm of 10^-6 to convert uM to M    
      
      ! Calculations associated with bacterial reactions.  It is assumed there will typically be
      ! far more bacteria than number of substrates
      do i=1,nBac
         ! Determine weighting vector based on B and concentrations of substrates. Calculate denominator, but insure >= absZero
         wS(i) = fZero(absZero, dot_product(B(i,1:nSub),sj(1:nSub))) 
      end do
      do j=1,nSub
         !$OMP SIMD
         do i=1, nBac
            ! Calculation reaction free energy (J/mmol)
            ! Anabolic reaction, 
            drG_A = drG0_rbA(j) + RkJ*T_K*( (1._mp - aAs(j))/alf_Bio*(log(bi(i))+ln106) + aAs(j)*(log(h2co3)+ln106) &
                  - (log(sj(j))+ln106)/CHO(j,1) - (1._mp - aAs(j))/alf_Bio*(del_Bio*(log(nh3)+ln106) + zet_Bio*(log(h3po4)+ln106)) )
            ! Catabolic Reaction
            drG_C = drG0_rbC(j) + RkJ*T_K*( (log(h2co3)+ln106) - (log(sj(j))+ln106)/CHO(j,1) - aCs(j)*(log(o2)+ln106) )
            ! Combined reaction
            drG_rb(i,j) = epsB(i)*drG_A + (1._mp - epsB(i))*drG_C
            ! Reaction rate
            k_eps4 = kappa*epsB(i)**4
            rb(i,j) = nuStar*epsB(i)**2*(B(i,j)*sj(j)/wS(i))*bi(i) &
                    *(sj(j)/(sj(j)+k_eps4/CHO(j,1)))*(o2/(o2+k_eps4*aCs(j)))*(nh3/(nh3+k_eps4*del_Bio*(1._mp - aAs(j))/alf_Bio)) &
                    *(h3po4/(h3po4+k_eps4*zet_Bio*(1._mp - aAs(j))/alf_Bio))*F_Thermo(drG_rb(i,j), 4._mp*aCs(j), T_k)
            ! Entropy production (J/d/K)
            sig_rb(i,j) = -V_L*drG_rb(i,j)*rb(i,j)/T_K 
         end do
      end do
            
      ! Calculations for predators, but do MPI column wise this time
      do i=1,nGrz
         ! First calculate the denominator for the weighting, so all processes have it
         wS(i) = dot_product(G(i,1:nBac),bi(1:nBac)) + dot_product(G(i,nBac+1:nBac+nGrz),gi(1:nGrz))
         wS(i) = fZero(absZero,wS(i))
      end do   
      ! First cover the bacteria
      do j=1, nBac
         !$OMP SIMD
         do i=1,nGrz
            drG_A = drG0_rgA + RkJ*T_K/alf_Bio*log(gi(i)/bi(j))
            ! Catabolic
            drG_C = drG0_rgC + RkJ*T_K*( (log(h2co3)+ln106) + ( del_bio*(log(nh3)+ln106) &
                  + zet_Bio*(log(h3po4)+ln106) - (log(bi(j))+ln106) )/alf_Bio - aCg*(log(o2)+ln106) )
            ! Combined reaction
            drG_rg(i,j) = epsG(i)*drG_A + (1._mp - epsG(i))*drG_C
            ! Reaction rate
            k_eps4 = kappa*epsG(i)**4
            rg(i,j) = nuStar*epsG(i)**2*(G(i,j)*bi(j)/wS(i))*gi(i) &
                     *(bi(j)/(bi(j)+k_eps4/alf_Bio))*(o2/(o2+k_eps4*aCg))*F_Thermo (drG_rg(i,j), 4._mp*aCg, T_k)
            ! Entropy production (J/d/K)
            sig_rg(i,j) = -V_L*drG_rg(i,j)*rg(i,j)/T_K 
         end do
      end do      
      ! now the grazers
      do j=nBac+1, nBac+nGrz
         !$OMP SIMD
         do i=1,nGrz
            drG_A = drG0_rgA + RkJ*T_K/alf_Bio*log(gi(i)/gi(j-nBac))
            ! Catabolic
            drG_C = drG0_rgC + RkJ*T_K*( (log(h2co3)+ln106) + ( del_bio*(log(nh3)+ln106) &
                  + zet_Bio*(log(h3po4)+ln106) - (log(gi(j-nBac))+ln106) )/alf_Bio - aCg*(log(o2)+ln106) )
            ! Combined reaction
            drG_rg(i,j) = epsG(i)*drG_A + (1._mp - epsG(i))*drG_C
            ! Reaction rate
            k_eps4 = kappa*epsG(i)**4
            rg(i,j) = nuStar*epsG(i)**2*(G(i,j)*gi(j-nBac)/wS(i))*gi(i) &
                     *(gi(j-nBac)/(gi(j-nBac)+k_eps4/alf_Bio))*(o2/(o2+k_eps4*aCg))*F_Thermo (drG_rg(i,j), 4._mp*aCg, T_k)
            ! Entropy production (J/d/K)
            sig_rg(i,j) = -V_L*drG_rg(i,j)*rg(i,j)/T_K 
         end do
      end do      
       
      ! Total entropy production
      sumSigmaDot = sum(sig_rb) + sum(sig_rg)
      return
   end subroutine rxnProperties   

   subroutine solout(m, t, x, f, k, ord, irtrn, rpar, ipar) 
      ! This routine is called by BiM.  It turns out that BiM runs much faster
      ! if you allow it to determine when to output points.  It will store points
      ! often when the state is changing fast, and slowley when not.  Makes for the
      ! smoothest output with the least amount of points (and it's faster)
      ! NOTE must use bim_SOLOUTmod.f with this SOLOUT routine.
      ! This routine is just a wrapper for saveSolution
      use realPrec
      use mpi
      use globalVars
      implicit none
      integer,  intent(in)   :: m       ! size of the problem
      real(8),  intent(in)   :: t(k)    ! current time
      real(8),  intent(in)   :: x(m,k)  ! current solution
      real(8),  intent(in)   :: f(m,k)  ! derivative, dx/dt
      integer,  intent(in)   :: k       ! block-size of method
      integer,  intent(in)   :: ord     ! order of method
      integer,  intent(out)  :: irtrn   ! return code (0 means everything OK)
      real(mp), intent(inout):: rpar(*) ! real vector set by user on call to BiM
      integer,  intent(inout):: ipar(*) ! integer vector set by user on call to BiM
      ! local declarations
      integer mpiErr

      ! The solution is not stored, but state variables bi and gi are checked
      ! to see if they are significantly greater than their initial values
      where(bi > sdtpMin*bi_ini) biF = 1._mp ! bi has had some growth
      where(gi > sdtpMin*gi_ini) giF = 1._mp ! gi has had some growth
      if (onlyCheckSoln) return

      if (myRank == 0) call saveSolution(t(k), m, x(1:m,k), rpar, ipar)
      irtrn = 0
      return
   end subroutine solout 

   subroutine setStateNames(x)
      ! This routine assigns useful names to state variables, and more importanly insures they are >= absZero
      use realPrec
      use globalVars
      use functions, only: fZero ! This is used to prevent state variable from being < absZero.  The transition to absZero occurs at 2*absZero
      implicit none
      real(mp), intent(in):: x(nState) ! vector of state variables defined for BiM
      ! Local declarations
      integer i
      ! Begin assignments
      do i=1,nSub
         sj(i)   = fZero(absZero, x(i))
      end do
      ! Now set the named concentration variables, there should be nChm of them (see parameter in params module)
      do i=1,nChm
         chm(i) = fZero(absZero,  x(nSub+i))
      end do
      ! Set the bacteria vector
      do i=1,nBac
         bi(i)   = fZero(absZero, x(nSub+nChm+i))
      end do
      ! Set the predator vector
      do i=1,nGrz
         gi(i)   = fZero(absZero, x(nSub+nChm+nBac+i))
      end do
      ! Integrated entropy produced over all reactions
      sumSigma = x(nSub+nChm+nBac+nGrz+1) 
      return
   end subroutine setStateNames   
      
   subroutine dsortArray (n, m, y, icol)
      ! This sorts an real(mp) nxm array, y(n,m) in ascending order based on
      ! icol column. All rows of y are moved correspondingly.
      use realPrec
      implicit none
      integer , intent(in)   :: n ! Row dimension of y
      integer , intent(in)   :: m ! column dimention of y
      real(mp), intent(inout):: y(n,m) ! The icol column is sorted and the corresponding row moves accordingly.
      integer , intent(in)   :: icol ! which column is used to sort the rows.
      ! local declarations
      integer i,j
      real(mp) rowMin(m)
    
      j = 1
      do while (j < n)
         rowMin = y(j,1:m)
         do i=j+1,n
            if (y(i,icol) < rowMin(icol)) then
               y(j,1:m) = y(i,1:m)
               y(i,1:m) = rowMin
               rowMin = y(j,1:m)
            end if
         end do
         j = j+1
      end do
      return   
   end subroutine dsortArray  
      
   subroutine trophicPos(M, trophicP)
      ! This routine calculates the trophic position for each consumer given the G matrix. The G matrix
      ! is a global variable, so it must be set before calling this routine
      ! This is based on Lavine1980 "Several measures of trophic structure applicable to complex food webs"
      use realPrec
      use globalVars
      implicit none
      real(mp), intent(in) :: M(nGrz,nBac+nGrz)
      real(mp), intent(out):: trophicP(nGrz) ! ! the trophic position of the nGrz consumers. Bacteria are all level 1, substrates 0.
      ! local declarations
      real(mp) Qmat(nBac+nGrz,nBac+nGrz), work(nBac+nGrz)
      integer i, ipiv(nBac+nGrz), info
      
      Qmat = 0.0_mp ! this matrix has all zeros in the first nBac  rows
      Qmat(nBac+1:nBac+nGrz,1:nBac+nGrz) = -M ! set the bottom of Q to -M.
      do i=1,nBac+nGrz
         Qmat(i,i) = Qmat(i,i) + 1.0_mp ! (this give I-Q found in Lavine1980)
      end do
      ! Use LAPACK (in intel MKL) to invert the (I-Q) matrix stored in Qmat. First generage the LU factorization
      call dgetrf( nBac+nGrz, nBac+nGrz, Qmat, nBac+nGrz, ipiv, info )
      if (info /= 0) then
         write(*,'()') 'Warning from trophicPos:: Qmat could not be LU factored!'
         trophicP = 0.0_mp ! just set to zero and return
         return
      end if
      call dgetri( nBac+nGrz, Qmat, nBac+nGrz, ipiv, work, nBac+nGrz, info )
      if (info /= 0) then
         write(*,'()') 'Warning from trophicPos:: Qmat could not be inverted!'
         trophicP = 0.0_mp ! just set to zero and return
         return
      end if
      ! a row in Qmat now sums to the trophic positoin of that organism. Note, the first nBack rows will all sum to one
      ! so that information is not passed back, just the position of the consumers
      do i=1,nGrz
         trophicP(i) = sum(Qmat(i+nBac,1:nBac+nGrz))
      end do      
      return
   end subroutine trophicPos