   module thermoData 
      ! This thermodynamic data is from Alberty2003 and his Mathematica 
      ! BiochemThermo`BasicBiochemData3.m data.
      ! For bacterial structure, have used data on yeast from Battley1998.
      ! Energy units are in kJ/mole
      ! V3: This version is identical to files labeled with V2.  Not sure why this was given that different name
      !    This file was copied from: C:\Users\jvallino\Projects\PIE-LTER\Modeling\PIE-1D-MEP\MEP-phototroph-1D-OC
      ! V4, 17-Aug-2016: Giving this a new version, just to distinguish it from V3.  A few changes were made
      !    Made R a global module parameter, and changed name to RJ to distinguish it from other R values.
      !    Used NIST values for R (also used Wikipedia's calculations for others.)
      !    Replaced data statements with parameters statements
      ! V4.1, 28-Aug-2016
      !     Added hydrogen (aq) species to parameter list
      ! V4.2, 20-Apr-2017: Nothing changed, but did add some new comments and want to keep those as versions go forward.
      !     No, now it has changed.  There is a slight descrepancy in Alberty2003 tables 8.1 and 8.2.  The values for co2aqsp have 
      !     been updated here.  Also, species for H2CO3, HCO3- and CO3-- have been added so that the fraction of h2co3, hco3 and co3 can be extracted
      !     from co2tot.  See freeHCO3 and freeCO3 functions.  It has been varified that freeCO2 + freeHCO3 + freeCO3 = 1 (see CarbonateChemistry2.nb in Mma directory) 
      !     Also, added free energy of formation for aveAA, aveNA, SiO2 and Si(OH)4
      ! V4.3, 21-Jul-2017: Added function to calculate H2S in equilibrium between air and water.
      ! V4.4, 17-May-2017: Added
      !     glucosesp, avePLsp, aveChlsp
      ! V5.0, 24-Jul-2018:
      !     Replaced parameter attribute of GHdata with target attribute
      !     Introduced derived type GHdb and array GHznData that points string names to sp data
      !     Introduced new function dGf0n(name,T,is,pH) that can use a string name for compound
      !     Introduced new routine initThermo(), that must be called to populate GHznData needed by dGf0n
      !     Introduced new parameters maxSpec, maxName and maxDim (see below).
      !     added compounds: acetatesp, ethanolsp, xylosesp
      !     Added gas constant for bar, Rbar
      ! V5.1, 22-Jan-2019:
      !     Added Mg(OH)2 (aq) and Mg2+ (aq) to the database for Chl a synthesis reaction.
      implicit none
      save ! Probably not necessary
      integer, parameter:: maxSpec = 5 ! This is the maximum number of species allowed for compound
      integer, parameter:: maxName = 10 ! This is the maximum length allowed for species names
      integer, parameter:: maxDim = 38  ! maximum allowed lenght of GHznData that holds the character names and pointer assignment    
      type GHdata
         real(8) dGf ! Standar free energy of formation, kJ/mol (25 C, 1bar)
         real(8) dHf ! Standard enthalpy of formation kJ/mol (25 C, 1 bar)
         real(8) zi  ! Charge on molecule
         real(8) nH  ! Number of hydrogen atoms at the given charge
      end type GHdata
      type GHdb
         ! This is introduced so that character names can be used to refer to species
         character(maxName) name
         type(GHdata), pointer:: sp(:)
      end type GHdb
      
      type(GHdb) GHznData(maxDim) ! This is needed to use string names in function dGf0n
      
      real(8), parameter:: Rp   = 0.082057338d0  ! gas constant (atm L/(g-mol K))
      real(8), parameter:: Rpm3 = 0.082057338d-6 ! gas constant (atm m^3/(g-mmol K)) (note, mmol not mol)
      real(8), parameter:: RkJ  = 8.3144598d-3   ! gas constant (kJ/(g-mol K) or J/(g-mmol K))
      real(8), parameter:: RJ   = 8.3144598d0    ! (J/g-mol/K) or (Pa m^3/(K g-mol))
      real(8), parameter:: Rbar = 8.3144598d-5   ! (bar m^3/K/g-mol) or (mbar m^3/K/g-mmol)
   
      ! Note, in order for all quantities to be treated as arrays, it is necessary to dimension species
      ! even if they only have one speces, such as h2osp(1).
      type(GHdata), target:: atpsp(3)		   = [GHdata(-2768.1,-3619.21,-4,12),GHdata(-2811.48,-3612.91,-3,13),GHdata(-2838.18,-3627.91,-2,14)]
      type(GHdata), target:: methaneaqsp(1)  = [GHdata(-34.33,-89.04,0,4)]
      type(GHdata), target:: methanegsp(1)	= [GHdata(-50.72,-74.81,0,4)]
      type(GHdata), target:: methanolsp(1)	= [GHdata(-175.31,-245.93,0,4)]
      type(GHdata), target:: methanolgsp(1)  = [GHdata(-162.3,-201.0,0,4)] ! gas phase data from CRC Handbook of Chemistry and Physics, Vol 91.
      type(GHdata), target:: adpsp(3)		   = [GHdata(-1906.13,-2626.54,-3,12),GHdata(-1947.1,-2620.94,-2,13),GHdata(-1971.98,-2638.54,-1,14)]
      type(GHdata), target:: pisp(2)		   = [GHdata(-1096.1,-1299.,-2,1),GHdata(-1137.3,-1302.6,-1,2)] ! phosphate h3po4
      type(GHdata), target:: h2osp(1)		   = [GHdata(-237.19,-285.83,0,2)]  
      type(GHdata), target:: o2gsp(1)		   = [GHdata(0.,0.,0.,0.)]
      type(GHdata), target:: o2aqsp(1)		   = [GHdata(16.4,-11.7,0,0)]
      type(GHdata), target:: yeastsp(1)      = [GHdata(-87.96,-133.09,0,1.613)] !Battley1998: CH1.613 O0.557 N0.158 P0.012 S0.003 K0.022 Mg0.003 Ca0.001
      type(GHdata), target:: ch2osp(1)		   = [GHdata(-152.65,-210.365,0,2)] ! this is glucose, but on unit C basis (it will be 6 times less than glucose)
      type(GHdata), target:: nadoxsp(1)		= [GHdata(0,0,-1,26)]
      type(GHdata), target:: nadredsp(1)	   = [GHdata(22.65,-31.94,-2,27)]
      type(GHdata), target:: sulfatesp(2)	   = [GHdata(-744.53,-909.27,-2,0),GHdata(-755.91,-887.34,-1,1)]
      type(GHdata), target:: ammoniasp(2)    = [GHdata(-26.5,-80.29,0,3),GHdata(-79.31,-132.51,1,4)]
      type(GHdata), target:: co2gsp(1)		   = [GHdata(-394.36,-393.5,0,0)]
      type(GHdata), target:: co2aqsp(1)      = [GHdata(-385.92,-413.80,0,0)] ! [GHdata(-385.97,-413.81,0,0)] ! this is not in the alberty database, but is given in Table 8.2, pg 151, alberty2003
                                                          ! This values are slighty different so that h2co3(sp) equals co2(sp) + h2o
      type(GHdata), target:: h2co3sp(1)      = [GHdata(-623.11,-699.63,0,2)]      
      type(GHdata), target:: hco3sp(1)	      = [GHdata(-586.77,-691.99,-1,1)] ! This is only the bicarbonate species (HCO3-)
      type(GHdata), target:: co3sp(1)	      = [GHdata(-527.81,-677.14,-2,0)] ! this is only the carbonate species (CO3--)
      type(GHdata), target:: co2totsp(3)	   = [GHdata(-527.81,-677.14,-2,0),GHdata(-586.77,-691.99,-1,1),GHdata(-623.11,-699.63,0,2)]
      type(GHdata), target:: nitratesp(2)    = [GHdata(-108.74,-205.,-1,0),GHdata(-111.25,-207.36,0,1)]
      type(GHdata), target:: nitritesp(2)	   = [GHdata(-32.2,-104.6,-1,0), GHdata(-50.6,-119.2,0,1)]
      type(GHdata), target:: n2aqsp(1)		   = [GHdata(18.7,-10.54,0,0)]
      type(GHdata), target:: h2saqsp(3)      = [GHdata(85.8,33.1,-2,0),GHdata(12.08,-17.6,-1,1),GHdata(-27.83,-39.7,0,2)]
      ! These data for H2S in g are from Thauer1977 (for the dGf), and also
      ! from the book "Chemical Principles" by P. Atkins and L. Jones, gotten from google books for dHf.
      ! These data should also be in "The NBS tables of chemical thermodynamic properties : 
      ! selected values for inorganic and C1 and C2 organic substances in SI units" by Donald D Wagman; et al.
      type(GHdata), target:: h2sgsp(1)		   = [GHdata(-33.56,-20.63,0.,2.)]
      ! New additions with V4.1
      type(GHdata), target:: h2aqsp(1)		   = [GHdata(17.6,-4.2,0.,2.)]
      ! New with version 4.2
      ! See excel file BIOMASS-2017.xlsx in C:\Users\jvallino\models\BIOENERG for calculations to derived aveAA and aveNA
      type(GHdata), target:: aveAAsp(1)      = [GHdata(-297.2, 0., 0.031, 7.06)] ! Aver AA in E. coli, with stoichiometry: C3.5 H7.06 O1.7 N Note, don't know delHf, but does not matter very much
      type(GHdata), target:: aveNAsp(1)      = [GHdata(-1332.7, 0., 0., 13.8)] ! Aver NA in E. coli, with stoichiometry: C9.6 H13.8 O7.8 N3.9 P Note, don't know delHf, but does not matter very muc
      type(GHdata), target:: sio2sp(1)       = [GHdata(-856.281, -910.700, 0., 0.)] ! From Gunnarsson2000
      type(GHdata), target:: h4sio4aqsp(1)   = [GHdata(-1309.181, -1461.722, 0., 4)] ! From Gunnarsson2000
      
      ! New with 4.4
      type(GHdata), target:: glucosesp(1)    = [GHdata(-915.9, -1262.19, 0, 12)]
      ! Here use Dihydrogeranylgeranyl diphosphate C20H38O7P2 (Kegg C17439) for phospholipids. 
      ! These data are from the Equilibrator site (http://equilibrator.weizmann.ac.il/compound?compoundId=C17439) 
      ! and they determined these values using Noor2013 group contribution methods.  No enthalphy of formation was available, so just set to zero.
      type(GHdata), target:: avePLsp(4)      = [GHdata(-1595.0, 0., 0., 38.), GHdata(-1584.9, 0., -1., 37.), GHdata(-1566.6, 0., -2.,36.), GHdata(-1524.3, 0., -3., 35.)]
      ! Chlorophyll a (C55H57O5N4Mg) is not in Equilibrator, but strangely enough, Chl b is.  
      ! Consequently, just using the thermo for Chl b here, as it can't be too different.
      type(GHdata), target:: aveChlsp(3)     = [GHdata(-289.1, 0, -1, 71), GHdata(-340.2, 0, 0, 72), GHdata(-341.9, 0, 1,73)]
      
      ! New with 5.0.  These are from {alberty2003}
      type(GHdata), target:: acetatesp(2)    = [GHdata(-369.31, -486.01, -1, 3), GHdata(-396.45, -485.76, 0, 4)]
      type(GHdata), target:: ethanolsp(1)    = [GHdata(-181.64, -288.3, 0, 6)]
      type(GHdata), target:: xylosesp(1)     = [GHdata(-750.49, -1045.94, 0, 10)]
      
      ! New with 5.1.  
      ! Mg(OH)2 (aq) from wikapedia, which  provided the delG, and the delH determined by addition from OH- and Mg2+
      type(GHdata), target:: mgoh2aqsp(1)    = [GHdata(-769.4, -926.8, 0, 2)]
      ! Mg2+ (aq) from Sawyer, McCarty, Parkin, Chemistry for Environmental Engineering and Science, 5th Edition, McGraw-Hill, 2003. 
      type(GHdata), target:: mg2aqsp(1)      = [GHdata(-454.8, -466.85, 2, 0)]

contains
      subroutine initThermo()
         ! This must be called before a call to the function dGf0n, as this populates GHznData array
         GHznData( 1)%name = 'atp';         GHznData( 1)%sp => atpsp
         GHznData( 2)%name = 'methaneaq';   GHznData( 2)%sp => methaneaqsp
         GHznData( 3)%name = 'methaneg';    GHznData( 3)%sp => methanegsp
         GHznData( 4)%name = 'methanol';    GHznData( 4)%sp => methanolsp
         GHznData( 5)%name = 'methanolg';   GHznData( 5)%sp => methanolgsp
         GHznData( 6)%name = 'adp';         GHznData( 6)%sp => adpsp
         GHznData( 7)%name = 'pi';          GHznData( 7)%sp => pisp
         GHznData( 8)%name = 'h2o';         GHznData( 8)%sp => h2osp
         GHznData( 9)%name = 'o2g';         GHznData( 9)%sp => o2gsp
         GHznData(10)%name = 'o2aq';        GHznData(10)%sp => o2aqsp
         GHznData(11)%name = 'yeast';       GHznData(11)%sp => yeastsp
         GHznData(12)%name = 'ch2o';        GHznData(12)%sp => ch2osp
         GHznData(13)%name = 'nadox';       GHznData(13)%sp => nadoxsp
         GHznData(14)%name = 'nadred';      GHznData(14)%sp => nadredsp
         GHznData(15)%name = 'sulfate';     GHznData(15)%sp => sulfatesp
         GHznData(16)%name = 'ammonia';     GHznData(16)%sp => ammoniasp
         GHznData(17)%name = 'co2g';        GHznData(17)%sp => co2gsp
         GHznData(18)%name = 'co2aq';       GHznData(18)%sp => co2aqsp
         GHznData(19)%name = 'h2co3';       GHznData(19)%sp => h2co3sp
         GHznData(20)%name = 'hco3';        GHznData(20)%sp => hco3sp
         GHznData(21)%name = 'co3';         GHznData(21)%sp => co3sp
         GHznData(22)%name = 'co2tot';      GHznData(22)%sp => co2totsp
         GHznData(23)%name = 'nitrate';     GHznData(23)%sp => nitratesp
         GHznData(24)%name = 'nitrite';     GHznData(24)%sp => nitritesp
         GHznData(25)%name = 'n2aq';        GHznData(25)%sp => n2aqsp
         GHznData(26)%name = 'h2saq';       GHznData(26)%sp => h2saqsp
         GHznData(27)%name = 'h2sg';        GHznData(27)%sp => h2sgsp
         GHznData(28)%name = 'aveAA';       GHznData(28)%sp => aveAAsp
         GHznData(29)%name = 'aveNA';       GHznData(29)%sp => aveNAsp
         GHznData(30)%name = 'h4sio4aq';    GHznData(30)%sp => h4sio4aqsp
         GHznData(31)%name = 'glucose';     GHznData(31)%sp => glucosesp
         GHznData(32)%name = 'avePL';       GHznData(32)%sp => avePLsp
         GHznData(33)%name = 'aveChl';      GHznData(33)%sp => aveChlsp
         GHznData(34)%name = 'acetate';     GHznData(34)%sp => acetatesp
         GHznData(35)%name = 'ethanol';     GHznData(35)%sp => ethanolsp
         GHznData(36)%name = 'xylose';      GHznData(36)%sp => xylosesp    
         GHznData(37)%name = 'mg(oh)2aq';   GHznData(37)%sp => mgoh2aqsp
         GHznData(38)%name = 'mg2+aq';      GHznData(38)%sp => mg2aqsp
         return
      end subroutine initThermo

      real(8) function dGf0n(name,T,is,pH)     
         character*(*), intent(in):: name        ! Name of compound, lowercase
         real(8),   intent(in):: T, is, pH   ! temperature (K), ionic strenght (M) and pH
         
         ! local declarations
         integer lloc
         
         ! Get location of GHzn data in GHznData based on string name
         lloc = findloc(GHznData(:)%name, value=name, dim=1)
         if (lloc == 0) then
            write(*,'(a)') 'Error, compound "'//trim(name)//'" does not exisit in thermodynamic GHznData data'
            dGf0n = 0.d0
            return
         end if         
         !dGf0n = dGf0(GHznData(lloc)%GHzn(1:GHznData(lloc)%nSpecies),T,is,pH)
         dGf0n = dGf0(GHznData(lloc)%sp,T,is,pH)
         return
      end function dGf0n

      real(8) function dGf0(sp,T,is,pH)
         ! This function returns the Gibbs energy of formation for the compound sp at temperature T (K),
         ! ionic strength, is (M), and pH in kJ/mol.  This is based on Alberty2003.
         !   sp:  contains, dGf, dHf, zi and nH at 298.15 and zero is and pH.
         ! 
         ! 30 Mar 2006: First written; Joe Vallino
         ! 21-Apr-2017: Allocating and deallocating eterms seems to chew up a bunch of CPU time. Instead, just use fix dimentsion
         !     But, if new items are added to GHDATA that are larger than 5, then eterms will need to be made larger
         implicit none
         type(GHdata) sp(:)
         real(8) T, is, pH
    
         ! Local declarations
         integer i, nsp
         real(8) gibbscoeff, dGzeroT, pHterm, istermG, gpfnsp, Tref
         !real(8), allocatable:: eterms(:) 
         real(8) eterms(5) 
         real(8) maxet
         real(8), parameter:: ln10 = log(10.d0)
        
         !**** Begin Program *****
         nsp = size(sp)
         !allocate ( eterms(nsp) )
         Tref = 298.15d0
         do i=1,nsp  ! sum for each species of the compound sp
            gibbscoeff=(9.20483d0*T)/1.D3-(1.284668d0*T**2)/1.D5+(4.95199d0*T**3)/1.D8 
            dGzeroT=(sp(i)%dGf*T)/Tref+sp(i)%dHf*(1-T/Tref)
            pHterm=-sp(i)%nH*pH*RJ*T*ln10/1000.d0
            istermG=gibbscoeff*(sp(i)%zi**2-sp(i)%nH)*sqrt(is)/(1.d0+1.6d0*sqrt(is));
            eterms(i)=-1.d3*(dGzeroT-pHterm-istermG)/(RJ*T)
         end do
         ! Need to bring the exponents down to values that can be calculated
         ! without overflow (underflow ok), so subtract the largest one.
         maxet = maxval(eterms(1:nsp))
         dGf0 = 0.d0
         do i=1,nsp
            dGf0 = dGf0 + exp(eterms(i) - maxet)
         end do
         dGf0 = -RJ*T*(log(dGf0) + maxet)/1.d3
         !deallocate ( eterms )
         return
      end function dGf0      
          
      real(8) function solCH4(T, is, pH)
        ! This function returns the solubiility of methane in water (uM/atm) at temperature T (K),
        ! pH and ionic strength is (M).  It is based on Alberty's values for CH4 in gas and aquious phase.
        ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
        ! but this function provides consistency with Alberty's data and approach, so will keep.
        implicit none
        real(8) T, pH, is
        ! Local declarations
        real(8) dGrxn
        ! Use the following reaction:  CH4(g) -> CH4(aq)
        dGrxn = dGf0(methaneaqsp,T,is,pH) - dGf0(methanegsp,T,is,pH)
     
        solCH4 = exp(-1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
     
        return
      end function solCH4

      real(8) function solCH3OH(T, is, pH)
         ! This function returns the solubiility of methanol in water (uM/atm) at temperature T (K),
         ! pH and ionic strength is (M).  It is based on Alberty's values for CH4 in aquious phase
         ! and CH4 in the gas phase from CRC data.
         ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
         ! but this function provides consistency with Alberty's data and approach, so will keep.
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dGrxn
         ! Use the following reaction:  CH3OH(g) -> CH3OH(aq)
         dGrxn = dGf0(methanolsp,T,is,pH) - dGf0(methanolgsp,T,is,pH)
        
         solCH3OH = exp(-1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
        
         return
      end function solCH3OH

      real(8) function solH2CO3(T, is, pH)
         ! This function returns the solubiility of carbon dioxide in water (uM/atm) at temperature T (K)
         ! pH, and ionic strength is (M).  That is, given the CO2 partial pressure in ATM, this returns the 
         ! total concentration of all carbonate species (ie., DIC = solH2CO3(T,is,pH)*pCO2 ).  This is from Alberty2003, pg 151.
         ! Rxn is:  H2CO3(aq) -> CO2(g) + H2O
         ! *** Note, since mass transfer coef are based on free CO2, use co2g2l to get mass transfer driver, that is, use
         ! flux = (pCO2*co2g2l(T, is, pH) - freeCO2(T, is, pH)*[DIC])
         ! NOT (pCO2*solH2CO3(T, is, pH) - [DIC])
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dGrxn
         dGrxn = dGf0(co2gsp,T,is,pH) + dGf0(h2osp,T,is,pH) - dGf0(co2totsp,T,is,pH)

         solH2CO3 = exp(1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM

         return
      end function solH2CO3    

      real(8) function freeCO2(T, is, pH)
         ! This function returns the ratio of CO2 to totalCO2 (i.e. DIC) at temperature T (K)
         ! pH, and ionic strength is (M).  freeCO2 is unitless.  That is, given the DIC concentration, this returns the 
         ! free CO2 in water ([CO2] = freeCO2*[DIC]). 
         ! Rxn is:  co2totsp -> h2co3
         ! Note, this use to be written as H2CO3(aq) -> CO2(aq) + H2O and co2aqsp was used, but instead this has been replaced with that above.
         ! see below for freeHCO3 and freeCO3 for similar calculations
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dG0rxn
         dG0rxn = dGf0(h2co3sp,T,is,pH) - dGf0(co2totsp,T,is,pH)
         ! now dGrxn = dG0rxn + RT ln ([h2co3]/[co2tot]), but a equilibrium, dGrxn = 0, so dG0rxn = - RT ln ([h2co3]/[co2tot])
         ! so freeCO2  = [h2co3]/[co2tot] = Exp( -dG0rxn/(RT) ), but dG0rxn is kJ/mol, so must use RkJ = RJ/1000 to get units to cancel
         freeCO2 = exp( -dG0rxn/(T*RkJ) ) ! This is unitless ([CO2]/[DIC])  

         return
      end function freeCO2    
    
      real(8) function freeHCO3(T, is, pH)
         ! This function returns the ratio of HCO3 to totalCO2 (i.e. DIC) at temperature T (K)
         ! pH, and ionic strength is (M).  freeHCO3 is unitless.  That is, given the DIC concentration, this returns the 
         ! free HCO3 in water ([HCO3] = freeHCO3*[DIC]). 
         ! Rxn is:  co2tot(aq) -> HCO3(-1) + H(+); however, the H(+) does not have to be added to the free energy calculations as described by Alberty
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dG0rxn
         dG0rxn = dGf0(hco3sp,T,is,pH) - dGf0(co2totsp,T,is,pH)
         freeHCO3 = exp( -dG0rxn/(T*RkJ) ) ! This is unitless ([HCO3]/[DIC])  
         return
      end function freeHCO3   
       
      real(8) function freeCO3(T, is, pH)
         ! This function returns the ratio of CO3 to totalCO2 (i.e. DIC) at temperature T (K)
         ! pH, and ionic strength is (M).  freeCO3 is unitless.  That is, given the DIC concentration, this returns the 
         ! free CO3 in water ([CO3] = freeHO3*[DIC]). 
         ! Rxn is:  co2tot(aq) -> CO3(2-) + 2H(+); however, the H(+) does not have to be added to the free energy calculations as described by Alberty
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dG0rxn
         dG0rxn = dGf0(co3sp,T,is,pH) - dGf0(co2totsp,T,is,pH)
         freeCO3 = exp( -dG0rxn/(T*RkJ) ) ! This is unitless ([HCO3]/[DIC])  
         return
      end function freeCO3   
       
      real(8) function co2g2l(T, is, pH)
         ! This function returns the ratio of free CO2 (uM) to pCO2 (atm) at temperature T (K)
         ! pH, and ionic strength is (M).  units are (uM/atm).  That is, given the pCO2 (atm) , this returns the 
         ! free CO2 in water ([CO2] = co2gas2l*pCO2).  
         ! Rxn is:  CO2(aq) -> CO2(g) 
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dGrxn
         dGrxn = dGf0(co2gsp,T,is,pH) - dGf0(co2aqsp,T,is,pH)
         co2g2l = exp(1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
         return
      end function co2g2l  
       
      real(8) function solO2(T, is, pH)
         ! This function returns the solubiility of oxygen in water (uM/atm) at temperature T (K),
         ! pH and ionic strength is (M).  It is based on Alberty's values for O2 in gas and aquious phase.
         ! However, these reactions are not affected by pH nor is, so could still have used NIST data,
         ! but this function provides consistency with Alberty's data and approach, so will keep.
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dGrxn
         ! Use the following reaction:  O2(g) -> O2(aq)
         dGrxn = dGf0(o2aqsp,T,is,pH) - dGf0(o2gsp,T,is,pH)
     
         solO2 = exp(-1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
     
         return
      end function solO2
   
      real(8) function solH2S(T, is, pH)
         ! This function returns the solubiility of hydrogen sulfide in water (uM/atm) at temperature T (K),
         ! pH and ionic strength is (M).  It is based on Alberty's values for H2S aquious phase and H2S in the gas phase mentioned above.
         implicit none
         real(8) T, pH, is
         ! Local declarations
         real(8) dGrxn
         ! Use the following reaction:  H2S(g) -> H2S(aq)
         dGrxn = dGf0(h2saqsp,T,is,pH) - dGf0(h2sgsp,T,is,pH)
     
         solH2S = exp(-1.d3*dGrxn/(T*RJ))*1.0d6*1.01325d0 ! Last number converts Bar to ATM
     
         return
      end function solH2S
      
      !    real(8) function solCH4_old(T)
      !      ! This function returns the solubiility of methane in water (uM/atm) at temperature T (K).
      !      ! Data from http://webbook.nist.gov/chemistry/
      !        implicit none
      !        real(8) T
      !        
      !        solCH4_old = 1.4d-3*dexp( 1700.*(1.0/T-1.0/298.0) )*1.0d6
      !        
      !        return
      !    end function solCH4_old            
      !    real(8) function solO2(T)
      !      ! This function returns the solubiility of oxygen in water (uM/atm) at temperature T (K).
      !      ! Data from http://webbook.nist.gov/chemistry/
      !        implicit none
      !        real(8) T
      !        
      !        solO2 = 1.3d-3*dexp( 1700.*(1.0/T-1.0/298.0) )*1.0d6
      !        
      !        return
      !    end function solO2
      !
      !    real(8) function solCO2(T)
      !        ! This function returns the solubiility of carbon dioxide in water (uM/atm) at temperature T (K).
      !        ! Data from http://webbook.nist.gov/chemistry/
      !        implicit none
      !        real(8) T
      !
      !        solCO2 = 3.4d-2*dexp( 2400.*(1.0/T-1.0/298.0) )*1.0d6
      !
      !        return
      !    end function solCO2    
    
   end module thermoData