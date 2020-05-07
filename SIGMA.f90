module Tools
IMPLICIT NONE
! "kind PARAMETER"  for DOuble precision
INTEGER, PUBLIC, PARAMETER        :: dp      = SELECTED_REAL_KIND(P=15)
REAL (KIND=dp), PUBLIC, PARAMETER :: pi      = 3.1415926535897932384_dp
REAL (KIND=dp), PUBLIC, PARAMETER :: clight  = 2.99792458e14_dp         ! Celerity of light c (micron s-1)
REAL (KIND=dp), PUBLIC, PARAMETER :: N_Avo      = 6.0221e23_dp          ! Avogadro constant

LOGICAL, PUBLIC                   :: verbose
LOGICAL, PUBLIC                   :: geom
REAL (KIND=dp)                    :: a0                                 ! monomer size
REAL (KIND=dp)                    :: Df                                 ! monomer size
LOGICAL, PUBLIC                   :: crt
LOGICAL, PUBLIC                   :: hyperion
REAL (KIND=dp), ALLOCATABLE       :: lam(:)                 ! wavelength
INTEGER                           :: nlam                   ! nr of wavelengths
CHARACTER (len=500)               :: data_path='/Users/lefevre/DUST/SIGMA/DATA/'!!!!! TO_BE_CHANGED
! ------------------------------------------------------------------------
! Mueller matrix to be filled, others are equal to zero
! ------------------------------------------------------------------------
TYPE MUELLER
	REAL (KIND=dp)                  :: F11(180)               ! Mueller coefficients related to phase function
	REAL (KIND=dp)                  :: F12(180)               ! Mueller coefficients
	REAL (KIND=dp)                  :: F22(180)               ! Mueller coefficients
	REAL (KIND=dp)                  :: F33(180)               ! Mueller coefficients
	REAL (KIND=dp)                  :: F44(180)               ! Mueller coefficients
	REAL (KIND=dp)                  :: F34(180)               ! Mueller coefficients
end type MUELLER
! ------------------------------------------------------------------------
! Output variables
! ------------------------------------------------------------------------
TYPE PARTICLE
	REAL (KIND=dp)                  :: rv                     ! radius
	REAL (KIND=dp)                  :: rvmin                  ! minimum radius
	REAL (KIND=dp)                  :: rvmax                  ! maximal radius
	REAL (KIND=dp)                  :: rho                    ! mass density in g.cm^-3
	REAL (KIND=dp)                  :: mass                    ! mass density in g.cm^-3
	REAL (KIND=dp), ALLOCATABLE     :: Kabs(:)                ! absorption coefficient
	REAL (KIND=dp), ALLOCATABLE     :: Ksca(:)                ! scattering coefficient
	REAL (KIND=dp), ALLOCATABLE     :: Kext(:)                ! extinction coefficient
	REAL (KIND=dp), ALLOCATABLE     :: qabs_size(:,:)         ! individual qabs by size
	REAL (KIND=dp), ALLOCATABLE     :: qsca_size(:,:)         ! individual qsca by size
	REAL (KIND=dp), ALLOCATABLE     :: r_size(:)              ! individual sizes
	REAL (KIND=dp), ALLOCATABLE     :: g(:)                   ! asymmetry parameter Henyey Greenstein
	REAL (KIND=dp), ALLOCATABLE     :: pol(:)                 ! average polarisation
	TYPE(MUELLER),  ALLOCATABLE     :: F(:)
END TYPE PARTICLE
end module Tools


program SIGMA
use Tools
IMPLICIT NONE
REAL (KIND=dp)                    :: lam1                   ! minimum wavelength
REAL (KIND=dp)                    :: lam2                   ! maximum wavelength
CHARACTER (len=500)               :: particlefile           ! name of output

INTEGER                           :: na                     ! nr of sizes for size distribution
REAL (KIND=dp)                    :: apow                   ! power law index
REAL (KIND=dp)                    :: fmax                   ! maximum fraction of vaccum for DHS
REAL (KIND=dp)                    :: porosity               ! porosity = 1 - filling_factor

INTEGER                           :: nm                     ! nr of grain material in a composite grain
INTEGER                           :: n_add                  ! nr of grain material to be added
INTEGER                           :: n_mix                  ! nr of grain material to be mixed
REAL (KIND=dp)                    :: V_ices                  ! volume fraction of ices

TYPE(PARTICLE)                    :: p
INTEGER                           :: i, j
INTEGER                           :: ilam
INTEGER                           :: nlines
INTEGER                           :: count_mix
INTEGER                           :: count_add
CHARACTER (len=100)               :: tmp
CHARACTER (len=100)               :: value

INTEGER, ALLOCATABLE              :: number(:)                      ! number associated to a component
CHARACTER*500, ALLOCATABLE        :: dust_type(:)										! type of component
CHARACTER*500, ALLOCATABLE        :: ref_index(:)
CHARACTER*3, ALLOCATABLE          :: rule(:)
CHARACTER*500, ALLOCATABLE		    :: sizedis(:)
REAL (KIND=dp), ALLOCATABLE       :: vfrac(:)                   ! volume fraction of each component
REAL (KIND=dp), ALLOCATABLE	      :: m_dust(:)
REAL (KIND=dp), ALLOCATABLE       :: a_min(:)                   ! minimum size
REAL (KIND=dp), ALLOCATABLE       :: a_max(:)                   ! maximum size

CHARACTER (len=500), ALLOCATABLE  :: dust_type_mix(:)
CHARACTER (len=500), ALLOCATABLE  :: ref_index_mix(:)
CHARACTER (len=500), ALLOCATABLE  :: sizedistrib_mix(:)
REAL (KIND=dp), ALLOCATABLE       :: vfrac_mix(:)
REAL (KIND=dp), ALLOCATABLE       :: a_min_mix(:)                   ! minimum size
REAL (KIND=dp), ALLOCATABLE       :: a_max_mix(:)                   ! maximum size

CHARACTER (len=500), ALLOCATABLE  :: dust_type_add(:)
CHARACTER (len=500), ALLOCATABLE  :: ref_index_add(:)
CHARACTER (len=500), ALLOCATABLE  :: sizedistrib_add(:)
REAL (KIND=dp), ALLOCATABLE       :: vfrac_add(:)
REAL (KIND=dp), ALLOCATABLE       :: a_min_add(:)                   ! minimum size
REAL (KIND=dp), ALLOCATABLE       :: a_max_add(:)                   ! maximum size

REAL (KIND=dp), ALLOCATABLE     :: Kabs_tot(:)                ! absorption coefficient
REAL (KIND=dp), ALLOCATABLE     :: Ksca_tot(:)                ! scattering coefficient
REAL (KIND=dp), ALLOCATABLE     :: Kext_tot(:)                ! extinction coefficient
REAL (KIND=dp), ALLOCATABLE     :: g_tot(:)
REAL (KIND=dp)                  :: Mtot

REAL (KIND=dp)                  :: lambda_ref                      ! reference wavelength for normalization

! ------------------------------------------------------------------------
! Default parameters - may be updated by users - to be moved to an external file?
! ------------------------------------------------------------------------
verbose        = .false.
geom           = .false.
a0             =  0.2_dp !default value for Min et al.2016
Df             =  3.0_dp !default value for Min et al.2016

crt            = .true.
hyperion       = .true.
apow           =  3.50_dp

fmax           =  0.0_dp
porosity       =  0.0_dp

V_ices          =  0.0_dp !by default no ice mantles

na             =  50
nlam           =  1000
lam1           =  0.05_dp
lam2           =  1000.0_dp
nm             =  3 !by default mixture of silicates, carbonaceous and iron sulfide
n_add          =  0
n_mix          =  0
particlefile   = "opacities"
lambda_ref     =  2.20 !by default K band is used

! ------------------------------------------------------------------------
! Input parameters
! ------------------------------------------------------------------------
CALL getarg(1,tmp)
i = 1
DO while(tmp.ne.' ')
	select CASE(tmp)
		CASE('-apow')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) apow
		CASE('-nm')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) nm
		CASE('-na')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) na
		CASE('-lmin')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) lam1
		CASE('-lmax')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) lam2
		CASE('-nlam')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) nlam
		CASE('-porosity')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) porosity
		CASE('-fmax')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) fmax
		CASE('-lambda_ref')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) lambda_ref
		CASE('-file','-filename')
			i = i+1
			CALL getarg(i,particlefile)
		CASE('-lref')
			i = i+1
			CALL getarg(i,value)
			READ(value,*) lambda_ref
                CASE('-df')
                        i = i+1
                        CALL getarg(i,value)
                        READ(value,*) Df
                CASE('-a0')
                        i = i+1
                        CALL getarg(i,value)
                        READ(value,*) a0
		CASE('-v','-verbose')
			verbose = .true.
		CASE('-geom')
			geom = .true.
		CASE('--help','-help','help')
			WRITE(*,'(" ")')
			WRITE(*,'("			SIGMA FEB-19 HELP			")')
			WRITE(*,'("		AVAILABLE OPTIONS ARE THE FOLLOWING:			")')
			WRITE(*,'(" ")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Dust Composition:")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	-nm = number of dust components apart from vacuum")')
			WRITE(*,'("")')
			WRITE(*,'("	Set according to DATA/input/COMPO.dat")')
			WRITE(*,'("	Examples of composition are provided in DATA/input folder")')
			WRITE(*,'("")')
			WRITE(*,'("	The different components could be mixed or added")')
			WRITE(*,'("	Individual components refractive index tables are located in DATA in the following subdirectories:")')
			WRITE(*,'("	* carbonaceous")')
			WRITE(*,'("	* silicates")')
			WRITE(*,'("	* iron")')
			WRITE(*,'("	* iron_sulfides")')
			WRITE(*,'("	* iron")')
			WRITE(*,'("	* ices")')
			WRITE(*,'("")')
			WRITE(*,'("	!!!! WARNING: FOR THE TIME BEING ICES MUST BE DEFINED AS THE LAST COMPONENT in DATA/input/COMPO.dat !!!!")')
			WRITE(*,'("")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Size distribution:")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Set individually for each component according to DATA/input/COMPO.dat")')
			WRITE(*,'("	Should be identical for mixture (mix) while it could be different for added opacities (add)")')
			WRITE(*,'("	Either power law (plaw) or custom size distribution can be used:")')
			WRITE(*,'("	Custom size distribution are taken from DATA/sizedistrib folder")')
			WRITE(*,'("	Filename without extension should be provided in DATA/input/COMPO.dat")')
			WRITE(*,'("	File content is a(microns), n(a)*a, porosity(see below)")')
			WRITE(*,'("")')
			WRITE(*,'("	Other parameters:")')
			WRITE(*,'("	-apow = will be used only if power law")')
			WRITE(*,'("	-na = number of bin sizes")')
			WRITE(*,'("")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Dust Shape:")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	-fmax = volume fraction of vacuum for hollow spheres (DHS)")')
			WRITE(*,'("	-porosity = constant porosity if plaw or if -1 is set in size distribution file")')
			WRITE(*,'("	Otherwise porosity can be defined as a function of size in size distribution file.")')
			WRITE(*,'("")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Other parameters:")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	-lmin = lambda min")')
			WRITE(*,'("	-lmax = lambda max")')
			WRITE(*,'("	-nlam = number of wavelength bins")')
			WRITE(*,'("	-lambda_ref = wavelength in micron to normalize extinction, by default 2.2")')
			WRITE(*,'("	-file or -filename = output filename, by default opacities")')
			WRITE(*,'("	-v = verbose")')
			WRITE(*,'("")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Default parameters :")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("")')
			WRITE(*,'("	You may find default parameters inside the code itself (in SIGMA.f90).")')
			WRITE(*,'("	You can update them to your needs in particular the following ones:")')
			WRITE(*,'("")')
			WRITE(*,'("	verbose   = .false.")')
			WRITE(*,'("	crt       = .true.")')
			WRITE(*,'("	hyperion  = .true.")')
			WRITE(*,'("")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	Outputs:")')
			WRITE(*,'("	========================================================")')
			WRITE(*,'("	By default output files are written in output folder.")')
			WRITE(*,'("	They will be overwritten if a different filename is not given as an argument.")')
			WRITE(*,'("")')
			WRITE(*,'("	SIGMA provides the following outputs:")')
			WRITE(*,'("")')
			WRITE(*,'("	* filename.dat = ASCII file that contains main outputs:")')
			WRITE(*,'("	lambda[micron]   Kabs[cm^2/g]   Ksca[cm^2/g]  Kext[cm^2/g]    albedo    asymmetry_parameter  linear_pol")')
			WRITE(*,'("")')
			WRITE(*,'("	* Kext_opacities.dat = extinction normalized by default by K band")')
			WRITE(*,'("")')
			WRITE(*,'("	OUTPUT RELATED TO SPECIFIC RADIATIVE TRANSFER CODES:")')
			WRITE(*,'("")')
			WRITE(*,'("	1) CRT (Juvela et al.):")')
			WRITE(*,'("	* CRT_opacities.nH2.dust")')
			WRITE(*,'("	* CRT_DSC_opacities.dat")')
			WRITE(*,'("")')
			WRITE(*,'("	2) HYPERION (Robitaille et al.):")')
			WRITE(*,'("	* HYPERION_set_filename.py")')
			WRITE(*,'("	Contains the full parameters needed to build hdf5 dust model that can be directly used by HYPERION")')
			WRITE(*,'("	It includes computation of Rosseland and Planck mean extinctions")')
			STOP
		CASE default
			WRITE(*,'("Option not recognized! SIGMA -help to list all available options.")')
			STOP
	end select
	i = i+1
	CALL getarg(i,tmp)
ENDDO

IF (verbose) THEN
	WRITE(*,'("========================================================")')
	WRITE(*,'("Running SIGMA to compute dust opacities")')
	WRITE(*,'("Copyright (c) 2018, C. Lefèvre, IRAM")')
	WRITE(*,'("E-mail: lefevre@iram.fr")')
	WRITE(*,'("License: one can freely use, modify, or redistribute all parts of the code")')
	WRITE(*,'("SIGMA is a commented and modified version of the DIANA Opacity code (P. Woitke, St Andrews)")')
	WRITE(*,'("========================================================")')
	WRITE(*,'(" Parameters:")')
	WRITE(*,'("  apow    = ",f14.4)')           apow
	WRITE(*,'("  na      = ",i14)')             na
	WRITE(*,'("========================================================")')
ENDIF

ALLOCATE (lam(nlam))
ALLOCATE (p%Kabs(nlam))
ALLOCATE (p%Kext(nlam))
ALLOCATE (p%Ksca(nlam))
ALLOCATE (p%qsca_size(nlam,na))
ALLOCATE (p%qabs_size(nlam,na))
ALLOCATE (p%r_size(na))
ALLOCATE (p%g(nlam))
ALLOCATE (p%F(nlam))
ALLOCATE (p%pol(nlam))
ALLOCATE (Kabs_tot(nlam))
ALLOCATE (Kext_tot(nlam))
ALLOCATE (Ksca_tot(nlam))
ALLOCATE (g_tot(nlam))

DO i = 1,nlam
	lam(i)=10.0_dp**(log10(lam1)+log10(lam2/lam1)*(i-1)/(nlam-1))
ENDDO

IF(verbose) WRITE(*,'("========================================================")')

ALLOCATE (number(nm))
ALLOCATE (dust_type(nm))
ALLOCATE (ref_index(nm))
ALLOCATE (sizedis(nm))
ALLOCATE (vfrac(nm))
ALLOCATE (m_dust(nm))
ALLOCATE (rule(nm))
ALLOCATE (a_min(nm))
ALLOCATE (a_max(nm))

CALL READ_COMPO(nm,number,dust_type,ref_index,vfrac,rule,sizedis,m_dust,a_min,a_max)

!Check number of bin sizes:
nlines = 0
IF (trim(sizedis(1)).NE."plaw") THEN
	OPEN (1, file = "DATA/sizedistrib/"//trim(sizedis(1))//".dat")
	DO
			READ (1,*, END=119)
			nlines = nlines + 1
	END DO
	119 CLOSE (1)
	IF (na.NE.nlines) THEN
		write(*,'("ERROR: ",i4," size bins found in DATA/sizedistrib/",a,".dat")') nlines, trim(sizedis(1))
		write(*,'("You should run SIGMA with -na ",i4)') nlines
		stop
	ENDIF
ENDIF

DO i=1,nm
	IF (trim(rule(i)).EQ."add") THEN
		n_add = n_add+1
	ELSE IF (trim(rule(i)).EQ."mix") THEN
		n_mix = n_mix+1
	ENDIF
END DO

ALLOCATE (dust_type_mix(n_mix))
ALLOCATE (ref_index_mix(n_mix))
ALLOCATE (vfrac_mix(n_mix))
ALLOCATE (sizedistrib_mix(n_mix))
ALLOCATE (a_min_mix(n_mix))
ALLOCATE (a_max_mix(n_mix))

ALLOCATE (dust_type_add(n_add))
ALLOCATE (ref_index_add(n_add))
ALLOCATE (vfrac_add(n_add))
ALLOCATE (sizedistrib_add(n_add))
ALLOCATE (a_min_add(n_add))
ALLOCATE (a_max_add(n_add))


count_mix = 0
count_add = 0
DO i = 1, (n_add+n_mix)
	IF (trim(rule(i)).EQ."mix") THEN
		count_mix                  = count_mix+1
		dust_type_mix(count_mix)   = trim(dust_type(i))
		ref_index_mix(count_mix)   = trim(ref_index(i))
		vfrac_mix(count_mix)       = vfrac(i)
		sizedistrib_mix(count_mix) = trim(sizedis(i))
		a_min_mix(count_mix)       = a_min(i)
		a_max_mix(count_mix)       = a_max(i)
	ELSE
		count_add                  = count_add+1
		dust_type_add(count_add)   = trim(dust_type(i))
		ref_index_add(count_add)   = trim(ref_index(i))
		vfrac_add(count_add)       = vfrac(i)
		sizedistrib_add(count_add) = trim(sizedis(i))
		a_min_add(count_add)       = a_min(i)
		a_max_add(count_add)       = a_max(i)
	ENDIF
ENDDO

!Safety check - for mixture the same size distribution need to be used
DO i=1,n_mix
	IF (trim(sizedistrib_mix(i)).NE.trim(sizedistrib_mix(1))) THEN
		write(*,*) "ERROR The same size distribution needs to be defined for mixture (mix option, sizedistrib in DATA/input/COMPO.dat)."
		stop
  ENDIF
	IF (a_min_mix(i).NE.a_min_mix(1)) THEN
		write(*,*) "ERROR The same minimum size needs to be defined for mixture (mix option, amin in DATA/input/COMPO.dat)."
		stop
  ENDIF
	IF (a_max_mix(i).NE.a_max_mix(1)) THEN
		write(*,*) "ERROR The same maximum size needs to be defined for mixture (mix option, amax in DATA/input/COMPO.dat)."
		stop
	ENDIF
	IF (m_dust(i).NE.m_dust(1).AND.trim(rule(1)).EQ."mix") THEN
		write(*,*) "ERROR The same Mdust/(100*MH) needs to be defined for mixture in DATA/input/COMPO.dat (mix option)"
		stop
	ENDIF
END DO

Kabs_tot = 0.0_dp
Ksca_tot = 0.0_dp
Kext_tot = 0.0_dp
g_tot = 0.0_dp

IF (n_mix.GT.0) THEN
	nm = n_mix
	CALL ComputePart(p,a_min_mix(1),a_max_mix(1),apow,fmax,vfrac_mix, porosity,na,&
	& nm,dust_type_mix,ref_index_mix,sizedistrib_mix(1))
	Kabs_tot = 0.0_dp
	Ksca_tot = 0.0_dp
	Kext_tot = 0.0_dp
	g_tot = 0.0_dp
	DO ilam = 1,nlam
		Kabs_tot(ilam) = Kabs_tot(ilam)+p%Kabs(ilam)*m_dust(1)
		Ksca_tot(ilam) = Ksca_tot(ilam)+p%Ksca(ilam)*m_dust(1)
		Kext_tot(ilam) = Kext_tot(ilam)+p%Kext(ilam)*m_dust(1)
		g_tot(ilam)    = g_tot(ilam)+p%g(ilam)*p%Ksca(ilam)*m_dust(1)
	ENDDO
ENDIF


IF (n_add.GT.0) THEN
	Kabs_tot = 0.0_dp
	Ksca_tot = 0.0_dp
	Kext_tot = 0.0_dp
	g_tot = 0.0_dp
	nm = 1
	DO i = 1,n_add
		CALL ComputePart(p,a_min_add(i),a_max_add(i),apow,fmax,vfrac_add(i), porosity,na,&
		& nm,dust_type_add(i),ref_index_add(i),sizedistrib_add(i))
		DO ilam = 1,nlam
			Kabs_tot(ilam) = Kabs_tot(ilam)+p%Kabs(ilam)*vfrac_add(i)*m_dust(i)
			Ksca_tot(ilam) = Ksca_tot(ilam)+p%Ksca(ilam)*vfrac_add(i)*m_dust(i)
			Kext_tot(ilam) = Kext_tot(ilam)+p%Kext(ilam)*vfrac_add(i)*m_dust(i)
			g_tot(ilam)    = g_tot(ilam)+p%g(ilam)*p%Ksca(ilam)*vfrac_add(i)*m_dust(i)
		ENDDO
	ENDDO
ENDIF

DO ilam = 1,nlam
	g_tot(ilam) = g_tot(ilam)/Ksca_tot(ilam)
ENDDO


CALL WriteOutput(particlefile,nlam,na,Kabs_tot,Ksca_tot,Kext_tot,g_tot,p, lambda_ref)


IF (verbose) WRITE(*,'("========================================================")')
IF (verbose) WRITE(*,'("Done! Have a nice day.")')
IF (verbose) WRITE(*,'("========================================================")')

end

SUBROUTINE WriteOutput (particlefile,nwav,ns,Kabs_tot,Ksca_tot,Kext_tot,g_tot,p,aref)
  use Tools
	IMPLICIT NONE
  TYPE(PARTICLE)                    :: p
	CHARACTER (len=500)               :: particlefile           ! name of output
	CHARACTER (len=500)               :: tablefile              ! name of output
	CHARACTER (len=500)               :: ddsfile                 ! name of output for debris disks simulator (dds)
	CHARACTER (len=500)               :: crtfile                ! name of output for CRT
	CHARACTER (len=500)               :: dscfile                ! name of output for DSC file for CRT
	CHARACTER (len=500)               :: extinctfile            ! name of output for normalized extinction
	CHARACTER (len=500)               :: phasefunctionfile      ! name of output
	CHARACTER (len=500)               :: lineardegpolfile      ! name of output
	CHARACTER (len=500)               :: hyperionfile           ! name of output for hyperion python file
	CHARACTER (len=500)               :: dust_hyperionfile       ! name of output for hyperion dust file
	INTEGER                           :: nwav
	INTEGER                           :: ns
	INTEGER                           :: i, ilam, k
	REAL (KIND=dp)                    :: Kabs_tot(nwav)
	REAL (KIND=dp)                    :: Ksca_tot(nwav)
	REAL (KIND=dp)                    :: Kext_tot(nwav)
	REAL (KIND=dp)                    :: g_tot(nwav)
	LOGICAL                           :: truefalse
	REAL (KIND=dp)                    :: aref
	REAL (KIND=dp)                    :: Kext_ref


	tablefile    = "output/" // trim(particlefile) // ".dat"
	ddsfile  = "output/dds_"// trim(particlefile) // ".dat"
	crtfile      = "output/CRT_" // trim(particlefile) // ".nH2.dust"
	dscfile      = "output/CRT_DSC_"// trim(particlefile) // ".dat"
	hyperionfile = "output/HYPERION_set_"//trim(particlefile)//".py"
	extinctfile  = "output/Kext_"// trim(particlefile) // ".dat"
	phasefunctionfile  = "output/PF_"// trim(particlefile) // ".dat"
	lineardegpolfile  = "output/LDP_"// trim(particlefile) // ".dat"
	dust_hyperionfile = "output/" // trim(particlefile)



  IF(verbose) THEN
		WRITE(*,'("Cross sections output to ascii table:   ",a)') trim(tablefile)
		IF(crt) WRITE(*,'("Optical properties output to ascii file for CRT: ",a)') trim(crtfile)
		IF(crt) WRITE(*,'("Full scattring matrix output to ascii file for CRT: ",a)') trim(dscfile)
		IF(hyperion) WRITE(*,'("Python program to set dust properties for HYPERION: ",a)') trim(hyperionfile)
  ENDIF

  ! Output file as ascii file
	inquire(file=tablefile,exist=truefalse)
	if(truefalse) then
		OPEN(unit=20,file=tablefile)
		IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(tablefile)
		CLOSE(unit=20,status='delete')
	endif
	OPEN(unit=20,file=tablefile,RECL=100000)
	write(20,'("#====================================================================#")')
	write(20,'("# Cross sections computed using SIGMA")')
	write(20,'("#=====================================================================================#")')
	write(20,'("# lambda[micron]   Kabs[cm^2/g]   Ksca[cm^2/g]  Kext[cm^2/g]    albedo    asymmetry_parameter  linear_pol")')
	write(20,'("#=====================================================================================#")')

	do ilam=1,nwav
		write(20,*) lam(ilam),Kabs_tot(ilam),Ksca_tot(ilam), Kext_tot(ilam), Ksca_tot(ilam)/Kext_tot(ilam), g_tot(ilam), p%pol(ilam)
	enddo
	CLOSE(unit=20)

	inquire(file=ddsfile,exist=truefalse)
	if(truefalse) then
		OPEN(unit=25,file=ddsfile)
		IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(ddsfile)
		CLOSE(unit=25,status='delete')
	endif
	OPEN(unit=25,file=ddsfile,RECL=100000)
	write(25,FMT='(a)') "#-------------------------------------------------------------------/"
	write(25,FMT='(a)') "# Tabulated dust properties: Q_abs, Q_sca, Q_ext, Albedo"
	write(25,FMT='(a)') "# as a function of wavelength for different dust grain radii"
	write(25,FMT='(a)') "#"
	write(25,FMT='(a)') "# ===> Please respect file formatting. Insert *no* comment lines!"
	write(25,FMT='(a)') "#-------------------------------------------------------------------/"
	write(25,FMT='(a,i4,a)') " 				", ns, "  # Number of different dust grain radii"
	write(25,FMT='(a,i4,a)') " 				", nwav, "  # Number of wavelengths tabulated for each grain radius"
	write(25,FMT='(a)') "#List of dust grain radii [micron], sorted monotonically"
	do k=1,ns
		write(25,FMT='(D15.3)') p%r_size(k)
	enddo
	write(25,FMT='(a)') "# List of wavelengths [micron], sorted monotonically"
	do ilam=1,nwav
		write(25,FMT='(D15.3)') lam(ilam)
	enddo
	do k=1,ns
		write(25,FMT='(a,D15.3,a)') "# Dust grain radius =    ",p%r_size(k)," micron"
		write(25,FMT='(a)') "# Q_abs           Q_sca             Q_ext           albedo"
		write(25,FMT='(a)') "#"
		do ilam=1,nwav
				write(25,*) p%qabs_size(ilam,k),p%qsca_size(ilam,k), p%qabs_size(ilam,k)+p%qsca_size(ilam,k), &
				& p%qsca_size(ilam,k)/(p%qabs_size(ilam,k)+p%qsca_size(ilam,k))
		enddo
	enddo
	CLOSE(unit=25)

  ! Output files needed to set CRT dust model (Juvela et al.)
	inquire(file=crtfile,exist=truefalse)
	if(truefalse.and.crt) then
		OPEN(unit=30,file=crtfile)
		IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(crtfile)
		CLOSE(unit=30,status='delete')
	endif
	IF(crt) THEN
		OPEN(unit=30,file=crtfile,RECL=100000)
	!*****Header's*****
		WRITE(UNIT=30,FMT='(a)') "eqdust" ! keyword in CRT dust file could be eqdust or simple
		WRITE(UNIT=30,FMT='(E10.3)') 2.0e-7_dp ! N(H2) weighted by respective abundance of each type of grain compare to other types
		WRITE(UNIT=30,FMT='(E10.3)') 1.0e-4_dp !CRT file need an arbitrary local size (a)
		WRITE(UNIT=30,FMT='(i4)') nwav ! number of frequencies
		WRITE(UNIT=30,FMT='(a)') "#"
	!*****Header's end*****
		do ilam=nwav,1,-1
			write(UNIT=30,FMT='(E15.3)',advance='no') clight/lam(ilam)
			write(UNIT=30,FMT='(E15.3)',advance='no') g_tot(ilam)
			write(UNIT=30,FMT='(E15.3)',advance='no') Kabs_tot(ilam)/(pi*1e-15_dp*N_Avo*100)
			write(UNIT=30,FMT='(E15.3)') Ksca_tot(ilam)/(pi*1e-15_dp*N_Avo*100)
		enddo
		CLOSE(unit=30)
	ENDIF

	IF(crt) THEN
		inquire(file=dscfile,exist=truefalse)
		if(truefalse.and.crt) then
			OPEN(unit=40,file=dscfile)
			IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(dscfile)
			CLOSE(unit=40,status='delete')
		endif
		OPEN(unit=40,file=dscfile,RECL=100000)
	!*****Header's*****
		write(40,'(I4)') nwav+1 !first colum is dedicated to angle number
	!*****Header's end*****
		write(40,'(I4)',advance='no') -999
	  DO ilam = 1,nwav-1
			write(40,'(E15.7)',advance="no") clight/lam(ilam)
		END DO
		write(40,'(E15.7)') clight/lam(nwav)
		write(40,'(I4)',advance='no') -999
		DO ilam = 1,nwav-1
			write(40,'(E15.7)',advance="no") g_tot(ilam) !dummy values not used by CRT
		END DO
		write(40,'(E15.7)') g_tot(nwav)
		! DO ilam = 1,nwav
		! 	write(40,'(E15.7)',advance="no") 0.0_dp
		! END DO
		DO i = 1, 180
				write(40,'(I3)',advance="no") i
			DO ilam = 1,nwav-1
				write(40,'(E15.7)',advance="no") p%F(ilam)%F11(i)
			END DO
			  write(40,'(E15.7)') p%F(nwav)%F11(i)
		END DO
		CLOSE(unit=40)
	ENDIF

	IF (hyperion) then
	  ! Output files needed to set hyperion dust model (Robitaille et al.)
		inquire(file=hyperionfile,exist=truefalse)
		if(truefalse.and.hyperion) then
			OPEN(unit=45,file=hyperionfile)
			IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(hyperionfile)
			CLOSE(unit=45,status='delete')
		endif
		OPEN(unit=45,file=hyperionfile,RECL=10000000)
		WRITE(UNIT=45,FMT='(a)') "import os"
		WRITE(UNIT=45,FMT='(a)') "from hyperion.dust import SphericalDust"
		WRITE(UNIT=45,FMT='(a)') "from numpy import *"
		WRITE(UNIT=45,FMT='(a)') "d = SphericalDust()"

		!***** frequencies *****
		write(45,FMT='(a)',advance="no") "d.optical_properties.nu = ["
		DO ilam = nwav,2,-1
			write(45,'(E15.7)',advance="no") clight/lam(ilam)
			write(45,FMT='(a)',advance="no") ","
		END DO
		write(45,'(E15.7)',advance="no") clight/lam(1)
		write(45,'(a)') "]"

		!***** albedo *****
		write(45,FMT='(a)',advance="no") "d.optical_properties.albedo = ["
		DO ilam = nwav,2,-1
			write(45,'(E15.7)',advance="no") Ksca_tot(ilam)/Kext_tot(ilam)
			write(45,FMT='(a)',advance="no") ","
		END DO
		write(45,'(E15.7)',advance="no") Ksca_tot(ilam)/Kext_tot(1)
		write(45,'(a)') "]"

		!***** extinction *****
		write(45,FMT='(a)',advance="no") "d.optical_properties.chi = ["
		DO ilam = nwav,2,-1
			write(45,'(E15.7)',advance="no") Kext_tot(ilam)
			write(45,FMT='(a)',advance="no") ","
		END DO
		write(45,'(E15.7)',advance="no") Kext_tot(1)
		write(45,'(a)') "]"

		!***** scattering angle *****
		write(45,FMT='(a)',advance="no") "d.optical_properties.mu = ["
		DO i = 180, 2, -1
			write(45,FMT='(a)',advance="no") "cos("
			write(45,'(I3)',advance="no") i
			write(45,FMT='(a)',advance="no") "*pi/180.)"
			write(45,FMT='(a)',advance="no") ","
		END DO
		write(45,FMT='(a)',advance="no") "cos(1*pi/180.)"
		write(45,'(a)') "]"

		!***** scattering matrix *****
		!convention of Code & Whitney (1995):
		! P1 (equivalent to S11), P2 (equivalent to S12), P3 (equivalent to S44), and P4 (equivalent to -S34).
		!Each of these variables should be specified as a 2-d array with dimensions (n_nu, n_mu),
		!where n_nu is the number of frequencies, and n_mu is the number of values of the cosine of the scattering angle
		write(45,FMT='(a)') "d.optical_properties.initialize_scattering_matrix()"

		DO i = 180, 2, -1
			write(45,FMT='(a)',advance="no") "d.optical_properties.P1[:,"
			write(45,'(I3)',advance="no") i-1
			write(45,FMT='(a)',advance="no") "] = ["
			DO ilam = nwav,2,-1
				write(45,'(E15.7)',advance="no") p%F(ilam)%F11(i)
				write(45,FMT='(a)',advance="no") ","
			END DO
				write(45,'(E15.7)') p%F(1)%F11(i)
				write(45,'(a)') "]"
		END DO

		DO i = 180, 2, -1
			write(45,FMT='(a)',advance="no") "d.optical_properties.P2[:,"
			write(45,'(I3)',advance="no") i-1
			write(45,FMT='(a)',advance="no") "] = ["
			DO ilam = nwav,2,-1
				write(45,'(E15.7)',advance="no") p%F(ilam)%F12(i)
				write(45,FMT='(a)',advance="no") ","
			END DO
				write(45,'(E15.7)') p%F(1)%F12(i)
				write(45,'(a)') "]"
		END DO

		DO i = 180, 2, -1
			write(45,FMT='(a)',advance="no") "d.optical_properties.P3[:,"
			write(45,'(I3)',advance="no") i-1
			write(45,FMT='(a)',advance="no") "] = ["
			DO ilam = nwav,2,-1
				write(45,'(E15.7)',advance="no") p%F(ilam)%F44(i)
				write(45,FMT='(a)',advance="no") ","
			END DO
				write(45,'(E15.7)') p%F(1)%F44(i)
				write(45,'(a)') "]"
		END DO

		DO i = 180, 2, -1
			write(45,FMT='(a)',advance="no") "d.optical_properties.P4[:,"
			write(45,'(I3)',advance="no") i-1
			write(45,FMT='(a)',advance="no") "] = ["
			DO ilam = nwav,2,-1
				write(45,'(E15.7)',advance="no") -p%F(ilam)%F34(i)
				write(45,FMT='(a)',advance="no") ","
			END DO
				write(45,'(E15.7)') p%F(1)%F34(i)
				write(45,'(a)') "]"
		END DO

		!***** set emissivities to compute Rosseland and Planck *****
		! extrapolation needed:
		! the frequency range should extend almost three orders of magnitude above the peak frequency for the coldest temperature
		! and one order of magnitude below the peak frequency for the hottest temperature.
	  write(45,'(a)') "d.optical_properties.extrapolate_nu(1.e1, 1.e20)"
		write(45,'(a)') "d.set_lte_emissivities(n_temp=1000,temp_min=10,temp_max=100000.)"

		!***** set hdf5 dust model *****
		write(45,'(a)') "d.write('"//trim(dust_hyperionfile)//".hdf5')"
	  !***** plot dust model *****
		write(45,'(a)') "d.plot('"//trim(dust_hyperionfile)//".png')"
		write(45,'(a)') "os.system('csh make_rosseland.csh "//trim(dust_hyperionfile)//"')"
		CLOSE(unit=45)
	ENDIF

		inquire(file=extinctfile,exist=truefalse)
		if(truefalse) then
			OPEN(unit=50,file=extinctfile)
			IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(extinctfile)
			CLOSE(unit=50,status='delete')
		endif

		DO i=1,nwav
			IF (lam(i)>aref-0.05_dp.and.lam(i)<aref+0.05_dp) then
				Kext_ref = Kext_tot(i)
				exit
			ENDIF
	  ENDDO
		Kext_tot = Kext_tot/Kext_ref

	  OPEN(unit=50,file=extinctfile,RECL=100000)
		 do ilam=1,nwav
	 		write(50,*) lam(ilam),Kext_tot(ilam)
	 	enddo
	 	CLOSE(unit=50)


		OPEN(unit=60,file=phasefunctionfile,RECL=100000)
		IF(verbose) write(*,'("Output file ",a," already exists, overwriting")') trim(phasefunctionfile)
		DO ilam = 1,nwav-1
			write(60,'(E15.7)',advance="no") lam(ilam)
		ENDDO
		write(60,'(E15.7)') lam(nwav)
		DO i = 1, 180
				write(60,'(I3)',advance="no") i
			DO ilam = 1,nwav-1
				write(60,'(E15.7)',advance="no") p%F(ilam)%F11(i)
			END DO
				write(60,'(E15.7)') p%F(nwav)%F11(i)
		END DO
		CLOSE(unit=60)
END

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

	SUBROUTINE READ_COMPO(ndust,number,type,refr_index,vtype,rule,sizedis,mdust,amin,amax)
		use Tools
		IMPLICIT NONE
		INTEGER                           :: ndust
		INTEGER                           :: nlines
		INTEGER                           :: i
		INTEGER							              :: number(ndust)                      ! number associated to a component
		CHARACTER*500							        :: type(ndust)										! type of component
		CHARACTER*500							        :: refr_index(ndust)
		INTEGER   							          :: sd_type(ndust)
		CHARACTER*3   							      :: rule(ndust)
		REAL (KIND=dp)						        :: vtype(ndust)                   ! volume fraction of each component
		REAL (KIND=dp)						        :: mdust(ndust)
		CHARACTER*500							        :: sizedis(ndust)
		REAL (KIND=dp)						        :: amin(ndust)
		REAL (KIND=dp)						        :: amax(ndust)

		!Check that nm is correct:
		nlines = 0
		OPEN (1, file = 'DATA/input/COMPO.dat')
		DO
		    READ (1,*, END=109)
		    nlines = nlines + 1
		END DO
		109 CLOSE (1)
		nlines = nlines - 2 !subtraction of command lines
		IF (ndust.NE.nlines) THEN
			write(*,'("ERROR: ",i2," components found in DATA/input/COMPO.dat")') nlines
			write(*,'("You should run SIGMA with -nm ",i2)') nlines
			stop
		ENDIF

		OPEN(99,file="DATA/input/COMPO.dat")
		read (99, fmt=* )
		read (99, fmt=* )
		do i=3, ndust+2
			read (99, fmt=* ) number(i-2), type(i-2), refr_index(i-2), vtype(i-2), rule(i-2), sizedis(i-2), mdust(i-2), amin(i-2), amax(i-2)
		end do
		CLOSE(99)
	END

	SUBROUTINE tellertje(i,n)
	use Tools
	IMPLICIT NONE
	INTEGER                           :: i,n,f

	IF(i.EQ.1.AND.verbose) WRITE(*,'("....................")')
	f=int(20.0_dp*dble(i)/dble(n))

	IF(20.0_dp*real(i-1)/real(n).LT.real(f) &
     &   .and.20.0_dp*real(i+1)/real(n).GT.real(f).and.verbose) then
		WRITE(*,'(".",$)')
		CALL flush(6)
	ENDIF

	IF(i.EQ.n.and.verbose) WRITE(*,*)

	return
	end

	SUBROUTINE ComputePart(p,amin,amax,apow,fmax,vfrac,porosity,nsubgrains,nm,d_type,ref_ind,sizedistrib)
		use Tools
		IMPLICIT NONE

		INTEGER                           :: na                     ! nr of sizes for size distribution
		REAL (KIND=dp)                    :: amin                   ! minimum size
		REAL (KIND=dp)                    :: amax                   ! maximum size
		REAL (KIND=dp)                    :: apow                   ! power law index
		CHARACTER*500                     :: sizedistrib            ! IF 0 power law - ELSE READ from file(s)
		REAL (KIND=dp)                    :: fmax                   ! maximum fraction of vaccum for DHS
		REAL (KIND=dp)                    :: porosity               ! porosity = 1 - filling_factorSI

		INTEGER                           :: nm                     ! nr of grain material in a composite grain

		REAL (KIND=dp)                    :: vfrac(nm)
		CHARACTER*500							        :: d_type(nm)										! type of component
		CHARACTER*500							        :: ref_ind(nm)
		REAL (KIND=dp)                    :: V_ices                  ! volume fraction of ices

		TYPE (PARTICLE)                   :: p

	  INTEGER                           :: ii, i, j, k, l
		INTEGER                           :: MAXMAT                   ! maximum number of grain material
		INTEGER, PARAMETER                :: n_ang              = 180 ! number of angles
		INTEGER                           :: iopac
		INTEGER                           :: nsubgrains
		INTEGER                           :: nlines
		INTEGER                           :: nf
		INTEGER                           :: ns
		INTEGER                           :: nmono
		INTEGER                           :: ilam
		INTEGER                           :: i_ice
		INTEGER                           :: Err
		INTEGER                           :: spheres
		INTEGER                           :: toolarge
		INTEGER                           :: READ_option

		REAL (KIND=dp)                    :: cext
		REAL (KIND=dp)                    :: csca
		REAL (KIND=dp)                    :: cabs
		REAL (KIND=dp)                    :: ksca_bis
		REAL (KIND=dp)                    :: kabs_bis
		REAL (KIND=dp)                    :: kabs_G
		REAL (KIND=dp)                    :: G
		REAL (KIND=dp)                    :: cabs_RGD
		REAL (KIND=dp)                    :: cabs_mono
		REAL (KIND=dp)                    :: cemie_mono
		REAL (KIND=dp)                    :: csmie_mono
		REAL (KIND=dp)                    :: QEXT
		REAL (KIND=dp)                    :: QSCA
		REAL (KIND=dp)                    :: QABS
		REAL (KIND=dp)                    :: GQSC
		REAL (KIND=dp)                    :: totA

	  REAL (KIND=dp), ALLOCATABLE       :: f11(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: f12(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: f22(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: f33(:,:)
	  REAL (KIND=dp), ALLOCATABLE       :: f34(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: f44(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: Mief11(:)
		REAL (KIND=dp), ALLOCATABLE       :: Mief12(:)
		REAL (KIND=dp), ALLOCATABLE       :: Mief22(:)
		REAL (KIND=dp), ALLOCATABLE       :: Mief33(:)
	  REAL (KIND=dp), ALLOCATABLE       :: Mief34(:)
		REAL (KIND=dp), ALLOCATABLE       :: Mief44(:)

		REAL (KIND=dp), ALLOCATABLE       :: mu(:)
		REAL (KIND=dp), ALLOCATABLE       :: M1(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: M2(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: S21(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: D21(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: r(:)
		REAL (KIND=dp), ALLOCATABLE       :: nr(:,:)
		REAL (KIND=dp), ALLOCATABLE       :: pr(:)
		REAL (KIND=dp), ALLOCATABLE       :: f(:)
		REAL (KIND=dp), ALLOCATABLE       :: wf(:)
		REAL (KIND=dp), ALLOCATABLE       :: rho(:)

		REAL (KIND=dp)                    :: rmie
		REAL (KIND=dp)                    :: lmie
		REAL (KIND=dp)                    :: e1mie
		REAL (KIND=dp)                    :: e2mie
		REAL (KIND=dp)                    :: csmie
		REAL (KIND=dp)                    :: cemie
		REAL (KIND=dp)                    :: KR
		REAL (KIND=dp)                    :: theta
		REAL (KIND=dp)                    :: dummy

		REAL (KIND=dp)                    :: maxf
		REAL (KIND=dp)                    :: lambda
		REAL (KIND=dp)                    :: lmin
		REAL (KIND=dp)                    :: lmax
		REAL (KIND=dp)                    :: minlog
		REAL (KIND=dp)                    :: maxlog
		REAL (KIND=dp)                    :: rad
		REAL (KIND=dp)                    :: r1
		REAL (KIND=dp)                    :: r2
		REAL (KIND=dp)                    :: tot, tot2
		REAL (KIND=dp)                    :: pow
		REAL (KIND=dp)                    :: Mass
                REAL (KIND=dp)                    :: Mass2
	        REAL (KIND=dp)                    :: Vol
		REAL (KIND=dp)                    :: rho_av
		REAL (KIND=dp)                    :: rho_av_no_ice
		REAL (KIND=dp)                    :: rho_av_no_por
		REAL (KIND=dp)                    :: rho_av_no_por_no_ice
		REAL (KIND=dp)                    :: rho_ice
		REAL (KIND=dp)                    :: rcore
		REAL (KIND=dp)                    :: wvno
		REAL (KIND=dp)                    :: scale

	  !Effective refractive index
		REAL (KIND=dp)                    :: e1blend
	  REAL (KIND=dp)                    :: e2blend
		REAL (KIND=dp)                    :: e1blend_ice
	  REAL (KIND=dp)                    :: e2blend_ice
	  REAL (KIND=dp)                    :: e1av
		REAL (KIND=dp)                    :: e2av
	  REAL (KIND=dp), ALLOCATABLE       :: e1(:,:,:)
	  REAL (KIND=dp), ALLOCATABLE       :: e2(:,:,:)
		REAL (KIND=dp), ALLOCATABLE       :: e1bis(:,:,:)
		REAL (KIND=dp), ALLOCATABLE       :: e2bis(:,:,:)
		REAL (KIND=dp), ALLOCATABLE       :: e1d(:)
		REAL (KIND=dp), ALLOCATABLE       :: e2d(:)
		REAL (KIND=dp), ALLOCATABLE       :: e1ice(:)
		REAL (KIND=dp), ALLOCATABLE       :: e2ice(:)
		REAL (KIND=dp), ALLOCATABLE       :: frac(:)
		REAL (KIND=dp), ALLOCATABLE       :: frac_record(:)

		COMPLEX (KIND=dp), ALLOCATABLE    :: epsj(:)
		COMPLEX (KIND=dp)                 :: m
		COMPLEX (KIND=dp)                 :: eps_eff
		COMPLEX (KIND=dp)                 :: min
		COMPLEX (KIND=dp)                 :: mav
		COMPLEX (KIND=dp)                 :: alpha

	  CHARACTER (len=3)                 :: meth
		CHARACTER (len=1)                 :: num
		CHARACTER (len=500)               :: input
		CHARACTER (len=500)               :: filename(100)
		CHARACTER (len=500)               :: grid
		CHARACTER (len=500)               :: tmp, tmp2
		CHARACTER (len=500)               :: lnkfile
		CHARACTER (len=500)               :: ice

	  LOGICAL                           :: truefalse
	  LOGICAL                           :: checkparticlefile

		REAL a,da
		INTEGER index_WD, index_type

		WRITE(meth,100)
	100	format('DHS')



	na=180
        MAXMAT = nm+1

	allocate(Mief11(na))
	allocate(Mief12(na))
	allocate(Mief22(na))
	allocate(Mief33(na))
	allocate(Mief34(na))
	allocate(Mief44(na))
	allocate(mu(na))
	allocate(M1(na,2))
	allocate(M2(na,2))
	allocate(S21(na,2))
	allocate(D21(na,2))

	allocate(frac(MAXMAT))
	allocate(epsj(MAXMAT))
	allocate(frac_record(MAXMAT))
	allocate(rho(MAXMAT))
	allocate(f11(nlam,na))
	allocate(f12(nlam,na))
	allocate(f22(nlam,na))
	allocate(f33(nlam,na))
	allocate(f34(nlam,na))
	allocate(f44(nlam,na))

	minlog=log10(amin)
	maxlog=log10(amax)
	pow=-apow
	maxf=fmax

	ns=nsubgrains
	if (amin.EQ.amax) then
		ns = 1
	endif

	allocate(e1(MAXMAT,nlam,ns))
	allocate(e2(MAXMAT,nlam,ns))
	allocate(e1bis(MAXMAT,nlam,ns)) !same but without porosity
	allocate(e2bis(MAXMAT,nlam,ns)) !same but without porosity


	nf=20
	if(maxf.eq.0e0) nf=1
	allocate(r(ns))
	allocate(nr(MAXMAT,ns))
	allocate(pr(ns))
	allocate(f(nf))
	allocate(wf(nf))
	allocate(e1d(nlam))
	allocate(e2d(nlam))
	allocate(e1ice(nlam))
	allocate(e2ice(nlam))




	! ------------------------------------------------------------------------
	! Definition of dust components volume fractions
	! ------------------------------------------------------------------------

	DO i = 1,nm
		frac(i) = vfrac(i)
		frac_record(i) = vfrac(i)
	ENDDO

  V_ices = 0.0_dp
	i_ice = 0

  DO i=1,nm
	  !Loop to deal with ices need to be here
		IF (trim(d_type(i)).EQ."ices") THEN
			ice = ref_ind(i)
			i_ice = i
			if (i.eq.nm) then
				nm = nm-1
			ELSE
				write(*,*) "ERROR : Ice mantles need to be defined as the last component in DATA/input/COMPO.DAT"
				stop
			ENDIF
		ENDIF
	ENDDO


! Normalisation mandatory to take into acount
! the case where the user already normalized all fractions by V_ices (06-may-2020)
! there is a second normalization with variable porosity later
! consistency has to be checked between the two !
	tot=0.0_dp
	do i=1,nm
		tot=tot+frac(i)
	enddo
	frac=frac/tot
	if (i_ice.gt.0) V_ices = vfrac(i_ice)/tot
	!modified 07-may-2020 to deal with the case V_refractory ne 100%

!First read size distribution:
tot=0d0
do l=1,nm
	if(ns.eq.1) then
		r(1)=10d0**((minlog+maxlog)/2d0)!/(1d0-porosity)**(1d0/3d0)
		nr(l,1)=frac(l)
		tot=tot+nr(l,1)*r(1)**3
	else
!Size distribution is defined here
!analytic expression for size distribution law:
		if(trim(sizedistrib).EQ."plaw") then
			do k=1,ns
				r(k)=10d0**(minlog &
			&				+(maxlog-minlog)*real(k-1)/real(ns-1))!/(1d0-porosity)**(1d0/3d0)
					nr(l,k)=r(k)**(pow+1d0)
					if (r(k).lt.amin.or.r(k).gt.amax) nr(l,k) = 0.0_dp
					tot=tot+nr(l,k)*r(k)**3 !normalized
			enddo
		!tabulated version added by CL
		ELSE !single size distribution
			OPEN(109,file="DATA/sizedistrib/"//trim(sizedistrib)//".dat")
			IF(verbose) write(*,*) "Size distribution from: DATA/sizedistrib/"//trim(sizedistrib)//".dat"
			DO k=1, ns
				READ(109,fmt=*) r(k), nr(l,k), pr(k) !input files units must be micron and number
				if (r(k).lt.amin.or.r(k).gt.amax) nr(l,k) = 0.0_dp
				tot=tot+nr(l,k)*r(k)**3
			end DO
			CLOSE(109)
		ENDIF
	endif
	do k=1,ns
		nr(l,k)=frac(l)*nr(l,k)/tot
	enddo

enddo

! Record size distribution as an output file to check

	! OPEN(unit=50,file="output/sizedis.dat",RECL=100000)

	! index_WD = 7 ! 7 = 3.1A, 22=4.4A, 16 = 5.5A, 25 = 5.5B
	! index_type = 1 ! 1 = silicates, 2 = carbonaceous
	! DO k=1,ns
		! a = r(k)*1.0e-4_dp
		! CALL GRAIN_DIST_WD01(index_WD,index_type,a,da) !Rv = 5.5B
		! WRITE(50,*) a*1.0e4_dp, da
		! WRITE(50,*) r(k), nr(1,k)
	! ENDDO
	! CLOSE(unit=50)


	DO i = 1,nm
		filename(i) = trim(d_type(i))//"/"//trim(ref_ind(i))
	ENDDO


IF (verbose)	write(*,'("========================================================")')
IF (verbose)	write(*,'("Refractive index tables used:")')
	do i=1,nm
		call RegridDataLNK_file(filename(i),lam(1:nlam),e1d(1:nlam),e2d(1:nlam),nlam,.true.,rho(i))
		do k=1,ns
			e1(i,1:nlam,k)=e1d(1:nlam)
			e2(i,1:nlam,k)=e2d(1:nlam)
			e1bis(i,1:nlam,k)=e1d(1:nlam)
			e2bis(i,1:nlam,k)=e2d(1:nlam)
		enddo
	end do
	if (i_ice.eq.nm+1) then
		filename(nm+1)= "ices/"//trim(ice)
		call RegridDataLNK_file(filename(nm+1),lam(1:nlam),e1d(1:nlam),e2d(1:nlam),nlam,.true.,rho_ice)
		e1ice(1:nlam)=e1d(1:nlam)
		e2ice(1:nlam)=e2d(1:nlam)
	endif


	deallocate(e1d)
	deallocate(e2d)


	min=dcmplx(1d0,0d0)


  IF(verbose) THEN
		write(*,'("========================================================")')
		write(*,'("Computing particle:")')
		write(*,'("Size: ",f10.3," - ",f10.3," micron")') amin,amax
		write(*,'("Shape: porosity = ",f8.2,"% , fmax = ",f8.3)') porosity*100,fmax
		write(*,'("========================================================")')
	ENDIF

	nm=nm+1
	e1(nm,1:nlam,1:ns)=1.0_dp
	e2(nm,1:nlam,1:ns)=0.0_dp
	rho(nm)=0.0_dp
	IF(trim(sizedistrib).EQ."plaw".OR.pr(1).LE.-1.0_dp) THEN
		frac(nm) = porosity
		frac(1:nm-1)=frac(1:nm-1)*(1.0_dp-porosity)

		! Normalisation mandatory to make the Bruggeman rule converge
		! Please check if the normalisation is still needed here
		tot=0.0_dp
		do i=1,nm
			tot=tot+frac(i)
		enddo
		frac=frac/tot

	ELSE
		frac_record(1:nm) = frac(1:nm)
		frac(1:nm) = 0.0_dp
	ENDIF

	if (V_ices.le.0.0_dp) then
		do k=1,ns
			if(trim(sizedistrib).NE."plaw".AND.pr(1).GT.-1.0_dp) THEN
				frac(1:nm-1)=frac_record(1:nm-1)*(1.0_dp-pr(k))
				frac(nm)=pr(k)
				! Normalization needs to be here if porosity different for each bin size
				! Otherwise Bruggeman rule does not converge
				tot=0.0_dp
				do j=1,nm
					tot=tot+frac(j)
				enddo
				frac=frac/tot
			endif
		do i=1,nlam
			do j=1,nm
				epsj(j)=((dcmplx(e1(j,i,k),e2(j,i,k))))**2
			end do
			IF (nm.eq.2.and.porosity.le.0) then
				continue
			ELSE

				if (nm.gt.2.and.porosity.gt.0.and.geom) then
	        !same without porosity to compute geometrical cross section
					!Minato (2006), Tazaki (2018)
					call brugg(frac/(1.0_dp-porosity),nm-1,epsj,eps_eff)
					e1bis(1,i,k)=dreal(cdsqrt(eps_eff))
					e2bis(1,i,k)=dimag(cdsqrt(eps_eff))
				else
					e1bis(1,i,k)=e1(1,i,k)
					e2bis(1,i,k)=e2(1,i,k)
				endif

				call brugg(frac,nm,epsj, eps_eff)
				e1(1,i,k)=dreal(cdsqrt(eps_eff))
				e2(1,i,k)=dimag(cdsqrt(eps_eff))

				IF(verbose.and.i.eq.1.and.k.eq.1) THEN
					write(*,'("========================================================")')
					write(*,'("No ices included in the dust mixture")')
					write(*,'("========================================================")')
				ENDIF


			ENDIF
			enddo
		enddo
	else
		do k=1,ns
			if(trim(sizedistrib).NE."plaw".AND.pr(1).GT.-1.0_dp) THEN
				frac=frac_record*(1.0_dp-pr(k))
				!print *, nm ! for record: in the past nm was lost at this step
				frac(nm)=pr(k)
				! print*, frac,pr(k)
				! Normalization needs to be here if porosity different for each bin size
				! Otherwise Bruggeman rule does not converge
				tot=0.0_dp
				do j=1,nm
					tot=tot+frac(j)
				enddo
				frac=frac/tot
			endif
		do i=1,nlam
			do j=1,nm
				epsj(j)=((dcmplx(e1(j,i,k),e2(j,i,k))))**2
			end do
			IF (nm.eq.2.and.porosity.le.0) then
				continue
			ELSE
				call brugg(frac,nm,epsj, eps_eff)
				e1blend=dreal(cdsqrt(eps_eff))
				e2blend=dimag(cdsqrt(eps_eff))
			ENDIF

			if (i_ice.eq.nm) then
				call maxgarn_2compo(e1blend,e2blend,e1ice(i),e2ice(i),V_ices,e1blend_ice,e2blend_ice)
				e1(1,i,k)=e1blend_ice
				e2(1,i,k)=e2blend_ice
				IF(verbose.and.i.eq.1.and.k.eq.1) THEN
					write(*,'("========================================================")')
					write(*,'("Maxwell Garnett rule is used to include ice in dust mixture")')
					write(*,'("========================================================")')
				ENDIF
			else
				e1(1,i,k)=e1blend
				e2(1,i,k)=e2blend
				IF(verbose.and.i.eq.1.and.k.eq.1) THEN
					write(*,'("========================================================")')
					write(*,'("Bruggeman rule is used to include ice in dust mixture")')
					write(*,'("========================================================")')
				ENDIF
			endif

			enddo
		enddo
	endif

	rho_av=0.0_dp
	rho_av_no_por=0.0_dp

	! Tested on 13/02/2020
	! Remained to be checked: differences between:
	! MG: non-porous ICES
	! BR: porous ICES
	do i=1,nm
		! This part is to compute rho without ice
		! from 1 to nm is only refractories + porosity
		rho_av        = rho_av + frac(i)*rho(i)
		rho_av_no_por = rho_av_no_por + frac(i) / (1.0_dp - porosity) * rho(i)
	enddo
	rho_av_no_ice = rho_av
	rho_av_no_por_no_ice = rho_av_no_por
	! variation for ice computed with MG rule (default)
	! i_ice defined before nm = nm+1 if porosity > 0
	if(i_ice.eq.nm) THEN
		rho_av = rho_av / (1.0_dp + V_ices) + V_ices / (1.0_dp + V_ices) * rho_ice
		rho_av_no_por = rho_av_no_por / (1.0_dp + V_ices) + V_ices / (1.0_dp + V_ices) * rho_ice
	endif
	rho(1) = rho_av
	rho(2) = rho_av_no_por
	nm     = 1 !the composite aggregate is now defined
	IF(verbose) THEN
		write(*,'("Average bulk density = ",f8.3, " g/cm3")') rho_av
		if (porosity.gt.0) write(*,'("Average bulk density without porosity = ",f8.3, " g/cm3")') rho_av_no_por
		if (V_ices.gt.0) write(*,'("Average bulk density without ice = ",f8.3, " g/cm3")') rho_av_no_ice
		if (porosity.gt.0.and.V_ices.gt.0) write(*,'("Average bulk density without ice &
		& and without porosity = ",f8.3, " g/cm3")') rho_av_no_por_no_ice
		write(*,'("========================================================")')
	ENDIF

	do i=1,nlam
		do j=1,n_ang
			f11(i,j)=0d0
			f12(i,j)=0d0
			f22(i,j)=0d0
			f33(i,j)=0d0
			f34(i,j)=0d0
			f44(i,j)=0d0
		enddo
	enddo

	if(nf.gt.1.and.maxf.gt.0.01e0) then
		call gauleg2(0.01e0,maxf,f(1:nf),wf(1:nf),nf)
	else if(maxf.eq.0e0) then
		f(1:nf)=0d0
		wf(1:nf)=1d0/real(nf)
	else
		f(1)=maxf
		wf(1)=1d0
	endif





	do ilam=1,nlam
		call tellertje(ilam,nlam)
		csca=0d0
		cabs=0d0
		ksca_bis=0d0
		kabs_bis=0d0
		cext=0d0
		Mass=0d0
                Mass2=0d0
		Vol=0d0

		do i=1,n_ang/2
			theta=(real(i)-0.5)/real(n_ang/2)*pi/2d0
			mu(i)=cos(theta)
		enddo


		do l=1,nm
		if(frac(l).eq.0d0) goto 10
			do k=1,ns
				r1=r(k)*(1.0_dp+V_ices)**(1.0_dp/3.0_dp)
				Err=0
				spheres=0
				toolarge=0
				do i=1,nf
					rad=r1/(1d0-f(i))**(1d0/3d0)
					m=dcmplx(e1(l,ilam,k),-e2(l,ilam,k)) !porosity added by CL
					wvno=2d0*pi/lam(ilam)

				if(f(i).eq.0d0) then
					spheres=1
					goto 20
				endif
				if(r1*wvno.gt.10000d0) then
					toolarge=1
					goto 20
				endif
				if(meth(1:3).eq.'DHS') then
		                        rcore=rad*f(i)**(1d0/3d0)
					call DMiLay(RCORE, rad, WVNO, m, min, MU, &
		     &                   NA/2, QEXT, QSCA, QABS, GQSC, &
		     &                   M1, M2, S21, D21, NA ,Err)
		                else
		                        rcore=rad*0.999

					call DMiLay(RCORE, rad, WVNO, min, m, MU, &
		     &                   NA/2, QEXT, QSCA, QABS, GQSC, &
		     &                   M1, M2, S21, D21, NA ,Err)
				endif
20		if(Err.eq.1.or.spheres.eq.1.or.toolarge.eq.1) then
				rad=r1
				rcore=rad
				rmie=rad
				lmie=lam(ilam)
				e1mie=e1(l,ilam,k)
				e2mie=e2(l,ilam,k)
				if(Err.eq.1.or.i.eq.1) then
					if(rmie/lmie.lt.5000d0) then
						call MeerhoffMie(rmie,lmie,e1mie,e2mie,csmie,cemie &
	     &								,Mief11,Mief12,Mief33,Mief34,n_ang)
					else
						call MeerhoffMie(rmie,rmie/5000d0,e1mie,e2mie,csmie,cemie &
	     &								,Mief11,Mief12,Mief33,Mief34,n_ang)
					endif
				endif

				Mief22=Mief11
				Mief44=Mief33

			else
				cemie=qext*pi*rad**2
				csmie=qsca*pi*rad**2
				do j=1,n_ang/2
					Mief11(j)=(M2(j,1)+M1(j,1))/csmie/wvno**2*2d0*pi
					Mief12(j)=(M2(j,1)-M1(j,1))/csmie/wvno**2*2d0*pi
					Mief22(j)=(M2(j,1)+M1(j,1))/csmie/wvno**2*2d0*pi
					Mief33(j)=(S21(j,1))/csmie/wvno**2*2d0*pi
					Mief34(j)=(-D21(j,1))/csmie/wvno**2*2d0*pi
					Mief44(j)=(S21(j,1))/csmie/wvno**2*2d0*pi
					Mief11(n_ang-j+1)=(M2(j,2)+M1(j,2))/csmie/wvno**2*2d0*pi
					Mief12(n_ang-j+1)=(M2(j,2)-M1(j,2))/csmie/wvno**2*2d0*pi
					Mief22(n_ang-j+1)=(M2(j,2)+M1(j,2))/csmie/wvno**2*2d0*pi
					Mief33(n_ang-j+1)=(S21(j,2))/csmie/wvno**2*2d0*pi
					Mief34(n_ang-j+1)=(-D21(j,2))/csmie/wvno**2*2d0*pi
					Mief44(n_ang-j+1)=(S21(j,2))/csmie/wvno**2*2d0*pi
				enddo
			endif



!	make sure the scattering matrix is properly normalized by adjusting the forward peak.
		tot=0d0
		tot2=0d0
		do j=1,n_ang
			tot=tot+Mief11(j)*sin(pi*(real(j)-0.5)/real(n_ang))
			tot2=tot2+sin(pi*(real(j)-0.5)/real(n_ang))
		enddo
		Mief11(1)=Mief11(1)+(tot2-tot)/sin(pi*(0.5)/real(n_ang))
		if(Mief11(1).lt.0d0) Mief11(1)=0d0

		do j=1,n_ang
			f11(ilam,j)=f11(ilam,j)+wf(i)*nr(l,k)*Mief11(j)*csmie
			f12(ilam,j)=f12(ilam,j)+wf(i)*nr(l,k)*Mief12(j)*csmie
			f22(ilam,j)=f22(ilam,j)+wf(i)*nr(l,k)*Mief22(j)*csmie
			f33(ilam,j)=f33(ilam,j)+wf(i)*nr(l,k)*Mief33(j)*csmie
			f34(ilam,j)=f34(ilam,j)+wf(i)*nr(l,k)*Mief34(j)*csmie
			f44(ilam,j)=f44(ilam,j)+wf(i)*nr(l,k)*Mief44(j)*csmie
		enddo
		cext=cext+wf(i)*nr(l,k)*cemie
		csca=csca+wf(i)*nr(l,k)*csmie
	  cabs=cabs+wf(i)*nr(l,k)*(cemie-csmie)
		Mass=Mass+wf(i)*nr(l,k)*rho(l)*4d0*pi*r1**3/3d0
    if (i.eq.1) Mass2=Mass2+nr(l,k)*rho(l)*4d0*pi*r1**3/3d0
		Vol=Vol+wf(i)*nr(l,k)*4d0*pi*r1**3/3d0
		p%r_size(k)=rad
		p%qabs_size(ilam,k)=(cemie-csmie)/(pi*rad**2)
		p%qsca_size(ilam,k)=csmie/(pi*rad**2)
	enddo
	enddo

	! Replace by geometrical cross section
	if (geom) then
		nmono = nint(r1**3/a0**3)
		if ( nmono.ge.16 ) then
			if (Df.gt.2.0_dp) then
				G = (4.27d0*nmono**(-0.315d0)*exp(-1.74d0/nmono**(0.243d0)))*nmono*pi*(a0)**2
			else
				G = (0.352d0+0.566d0*nmono**(-0.138d0))*nmono*pi*(a0)**2
			endif
		else
			G = (12.5d0*nmono**(-0.315d0)*exp(-2.53d0/nmono**(0.0920d0)))*nmono**(2.0_dp/3.0_dp)*pi*(a0)**2
		end if
		call MeerhoffMie(a0,lam(ilam),e1bis(1,ilam,1),e2bis(1,ilam,1),csmie_mono,cemie_mono &
	&								,Mief11,Mief12,Mief33,Mief34,n_ang)
		cabs_mono = cemie_mono-csmie_mono
		cabs_RGD = cabs_mono*nmono

	        kabs_G = G*(1d0-exp(-cabs_RGD/G))*1d4/(Mass2*rho(2)/rho(1))

		kabs_bis=max(kabs_G, cabs*1d4/(Mass))
		!kabs_bis=kabs_G  !Uncomment only to compute pure Minato contribution
		ksca_bis=1d4*cext/(Mass)-kabs_bis
	endif

10	continue
	enddo

	p%rho=Mass/Vol
	!if (ilam.eq.1) print *, "Mass, Vol, r1, rho =", Mass, Vol, r1, p%rho
	! print needed to check ice normalization
	p%mass= Mass
	if (geom) then
		p%Kabs(ilam)=kabs_bis
		p%Ksca(ilam)=ksca_bis
	else
		!mdust_with_ice/mdust_without_ice  = 1 + V_ices * rho_ice / rho_without_ice
		p%Kabs(ilam)=1d4*cabs/(Mass)*(1.0_dp+V_ices*rho_ice/rho_av_no_ice) !correction added on 05-may-2020
		p%Ksca(ilam)=1d4*csca/(Mass)*(1.0_dp+V_ices*rho_ice/rho_av_no_ice)  !correction added on 05-may-2020
		!if (ilam.eq.1) print *, "corr factor", (1.0_dp+V_ices*rho_ice/rho_av_no_ice)
		!print needed to check ice normalization
	endif
	p%Kext(ilam)=1d4*cext/(Mass)*(1.0_dp+V_ices*rho_ice/rho_av_no_ice) !correction added on 05-may-2020
	!p%r_size(k) = ((3d0*Vol)/(4d0*pi))**(1d0/3d0) !commented out on 17/02/2020 - to be checked with variable porosity
	p%F(ilam)%F11(1:180)=f11(ilam,1:180)/csca
	p%F(ilam)%F12(1:180)=f12(ilam,1:180)/csca
	p%F(ilam)%F22(1:180)=f22(ilam,1:180)/csca
	p%F(ilam)%F33(1:180)=f33(ilam,1:180)/csca
	p%F(ilam)%F34(1:180)=f34(ilam,1:180)/csca
	p%F(ilam)%F44(1:180)=f44(ilam,1:180)/csca
	tot = 0.0_dp
	do i=1,180
		p%g(ilam)=p%g(ilam)+p%F(ilam)%F11(i)*cos(pi*(real(i)-0.5)/180d0) &
	 &					*sin(pi*(real(i)-0.5)/180d0)
		tot=tot+p%F(ilam)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
		p%pol(ilam)=p%pol(ilam)-p%F(ilam)%F12(i)/p%F(ilam)%F11(i)*sin(pi*(real(i)-0.5)/180d0)/(115d0)
	enddo
	p%g(ilam)=p%g(ilam)/tot
	enddo


#ifdef USE_FITSIO
	call ParticleFITS(p,r,nr(1:nm,1:ns),nm,ns,rho_av,amin,amax,apow,fmax,porosity,frac,rho,filename)
#endif


	deallocate(e1)
	deallocate(e2)

	deallocate(Mief11)
	deallocate(Mief12)
	deallocate(Mief22)
	deallocate(Mief33)
	deallocate(Mief34)
	deallocate(Mief44)
	deallocate(mu)
	deallocate(M1)
	deallocate(M2)
	deallocate(S21)
	deallocate(D21)

	deallocate(frac)
	deallocate(rho)
	deallocate(f11)
	deallocate(f12)
	deallocate(f22)
	deallocate(f33)
	deallocate(f34)
	deallocate(f44)

	deallocate(r)
	deallocate(nr)
	deallocate(f)
	deallocate(wf)


	return
	end



#ifdef USE_FITSIO
	subroutine ParticleFITS(p,r,nr,nm,na,rho_av,amin,amax,apow,fmax,porosity,frac,rho,lnkfiles)
	use Tools
	IMPLICIT NONE

	character*500 lnkfiles(nm)
	real*8 amin,amax,apow,fmax,porosity
	real frac(nm),rho(nm)
	logical blend
	character*6 word

	type(particle) p
	integer nm,na,i,j,ii,iopac,nm2
	real r(na),nr(nm,na)
	real a0,a1,a2,a3,rho_av,rmin,rmax
	real*8,allocatable :: array(:,:,:)

	  integer status,unit,blocksize,bitpix,naxis,naxes(3)
	  integer group,fpixel,nelements
	  logical simple,extend,truefalse

	inquire(file=particlefile,exist=truefalse)
	if(truefalse) then
		write(*,'("FITS file already exists, overwriting")')
		OPEN(unit=90,file=particlefile)
		CLOSE(unit=90,status='delete')
	endif

	blend=.true.

	a0=0d0
	a1=0d0
	a2=0d0
	a3=0d0
	rmin=r(1)
	rmax=r(1)
	nm2=nm
	if(blend) nm2=1
	do i=1,nm2
	do j=1,na
		a0=a0+nr(i,j)
		a1=a1+nr(i,j)*r(j)
		a2=a2+nr(i,j)*r(j)**2
		a3=a3+nr(i,j)*r(j)**3
		if(r(j).lt.rmin) rmin=r(j)
		if(r(j).gt.rmax) rmax=r(j)
	enddo
	enddo
	a1=a1/a0
	a2=a2/a0
	a3=a3/a0

	p%rv=sqrt(a2)*1d-4
	p%rvmin=amin*1d-4
	p%rvmax=amax*1d-4

	  status=0
!	 Get an unused Logical Unit Number to use to create the FITS file
	  call ftgiou(unit,status)
!	 create the new empty FITS file
	  blocksize=1
	  call ftinit(unit,particlefile,blocksize,status)

	  simple=.true.
	  extend=.true.
	group=1
	fpixel=1

	bitpix=-64
	naxis=2
	naxes(1)=nlam
	naxes(2)=4
	nelements=naxes(1)*naxes(2)
	allocate(array(nlam,4,1))

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write read_optional keywords to the header

	call ftpkye(unit,'r_min',real(amin),8,'[micron]',status)
	call ftpkye(unit,'r_max',real(amax),8,'[micron]',status)
	call ftpkye(unit,'r_pow',real(apow),8,'',status)
	call ftpkye(unit,'f_max',real(fmax),8,'',status)

	call ftpkye(unit,'a1',real(a1),8,'[micron]',status)
	call ftpkye(unit,'a2',real(a2),8,'[micron^2]',status)
	call ftpkye(unit,'a3',real(a3),8,'[micron^3]',status)
	call ftpkye(unit,'density',real(rho_av),8,'[g/cm^3]',status)

	if(blend) call ftpkye(unit,'porosity',real(porosity),8,'[g/cm^3]',status)

	do i=1,nm
		write(word,'("file",i0.2)') i
		call ftpkys(unit,word,trim(lnkfiles(i)),'',status)
	enddo
	do i=1,nm
		write(word,'("frac",i0.2)') i
		call ftpkye(unit,word,real(frac(i)),8,'[volume fraction]',status)
	enddo
	do i=1,nm
		write(word,'("rho",i0.2)') i
		call ftpkye(unit,word,real(rho(i)),8,'[g/cm^3]',status)
	enddo

	call ftpkyj(unit,'n_radii',na,' ',status)
	call ftpkyj(unit,'n_mat',nm,' ',status)



	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: opacities
	!------------------------------------------------------------------------------

	do i=1,nlam
		array(i,1,1)=lam(i)
		array(i,2,1)=p%Kext(i)
		array(i,3,1)=p%Kabs(i)
		array(i,4,1)=p%Ksca(i)
	enddo

	call ftpprd(unit,group,fpixel,nelements,array(1:nlam,1:4,1),status)

	deallocate(array)

	!------------------------------------------------------------------------------
	! HDU 1: Temperature
	!------------------------------------------------------------------------------
	bitpix=-64
	naxis=3
	naxes(1)=nlam
	naxes(2)=6
	naxes(3)=180
	nelements=naxes(1)*naxes(2)*naxes(3)

	allocate(array(nlam,6,180))

	! create new hdu
	call ftcrhd(unit, status)

	!  Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	do i=1,nlam
		do j=1,180
			array(i,1,j)=p%F(i)%F11(j)
			array(i,2,j)=p%F(i)%F12(j)
			array(i,3,j)=p%F(i)%F22(j)
			array(i,4,j)=p%F(i)%F33(j)
			array(i,5,j)=p%F(i)%F34(j)
			array(i,6,j)=p%F(i)%F44(j)
		enddo
	enddo

	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,array,status)

	deallocate(array)

	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   print*,'error in export to fits file',status
	end if


	return
	end
#endif


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

	subroutine brugg(f, nm, e, epsavg)

	!***********************************************************************
	!  This subroutine calculates the average dielectric function.
	!  It uses a generalized version of the Bruggeman dielectric mixing function.
	!  This routine was first coded by:
	!
	!  July 2000 : Benjamin T. Johnson, Atmospheric and Oceanic Sciences Dept.
	!              University of Wisconsin - Madison
	!              jbenjam@aos.wisc.edu
	!
	!  Feb. 2002 : Modifications by Michael A. Walters, SSEC, UW-Madison
	!              walters@rain.aos.wisc.edu
	!              Rewritten to use complex variables in input/output.
	!              Converted to Fortran 90
	!
	!  Jul 2018  : Modified by C. Lefèvre, IRAM, France
	!              lefevre@iram.fr
	!              Generalization to nm grain material and rewritten to accept complex arrays
	!              Integrated to SIGMA code to compute dust properties
	!
  !  Jul 2019   : Generalization to n components by Michiel Min and Charlène Lefèvre
	!  The subroutine will:
	!  1. Accept the parameters f, eps, nm
	!
	!  2. Calculate the average or mixed dielectric constant (epsavg), and
	!     return it back to the calling program.
	!
	!  Variable and parameter descriptions:
	!     nm      = number of grain material
	!     f    = volume fraction of each component
	!     e     = dielectric constant of each component
	!     epsavg  = averaged dielectric constant
	!**********************************************************************

		use Tools
		implicit none
		INTEGER               :: nm, i, k, l, m
		INTEGER               :: j(nm+1)
		REAL(KIND =dp)        :: f(nm)
		COMPLEX (KIND=dp)     :: e(nm)
		COMPLEX (KIND=dp)     :: epsavg
		COMPLEX (KIND=dp)     :: c(nm+1)
		COMPLEX (KIND=dp)     :: prod
		COMPLEX (KIND=dp)     :: x(nm)
		COMPLEX (KIND=dp)     :: roots(nm)
		COMPLEX (KIND=dp)     :: total !to check the result of Bruggeman rule
		LOGICAL polish
		polish=.false.

		c=0d0
		do i=1,nm
			x=-e/2d0
			x(i)=e(i)

			c(nm+1)=c(nm+1)+f(i)
			do k=1,nm
				do l=1,k
					j(l)=l
				enddo
				j(k+1)=nm+1
	1			continue
					prod=1.0_dp
					do l=1,k
						prod=prod*x(j(l))
					enddo
					c(nm-k+1)=c(nm-k+1)+f(i)*prod*(-1.0_dp)**k
					do l=1,k
						if((j(l)+1).lt.j(l+1)) then
							j(l)=j(l)+1
							do m=1,l-1
								j(m)=m
							enddo
							goto 1
						endif
					enddo
				continue
			enddo
		enddo

		call zroots(c,nm,roots,polish)
		do i=1,nm
			if(real(roots(i)).gt.0d0.and.dimag(roots(i)).gt.0d0) THEN
				epsavg=roots(i)
			else if (roots(i).eq.roots(1).and.dimag(roots(i)).lt.0d0) then
				write(*,*) "ERROR Bruggeman rule could cannot converge: no positive solution for effective refractive index: "
				write(*,*) roots(i)
				write(*,*) "Please try with a restricted range by using -lmin and -lmax"
				write(*,FMT='(a,f6.1,a,f6.1)') "currently lmin = ", lam(1), ", lmax = ", lam(nlam)
				stop
			endif
		enddo

		total = 0.0_dp
		do i = 1,nm
			total = total + (e(i)-epsavg)/(e(i)+2.0_dp*epsavg)*f(i)
		enddo

		if (abs(dreal(total)).gt.1e-15_dp.or.abs(dimag(total)).gt.1e-15_dp) THEN
			write(*,*) "ERROR Bruggeman rule did not converge"
			write(*,*) total
			stop
		ENDIF
		return
	end subroutine brugg


SUBROUTINE zroots(a,m,roots,polish)
	use Tools
	INTEGER m,MAXM
	REAL(KIND=dp) EPS
	COMPLEX(KIND=dp) a(m+1),roots(m)
	LOGICAL polish
	PARAMETER (EPS=1.e-15,MAXM=101) !A small number and maximum anticipated value of m+1.
	!USES laguer
	!Given the degree m and the complex coefficients a(1:m+1) of the polynomial  m+1 a(i)xi−1, i=1
	!this routine successively calls laguer and finds all m complex roots. The logical variable polish should be input as .true. if polishing (also by Laguerre’s method) is desired, .false. if the roots will be subsequently polished by other means.
	INTEGER i,j,jj,its
	COMPLEX(KIND=dp) ad(MAXM),x,b,c
	do j=1,m+1
		ad(j)=a(j)
	enddo
	!Copy of coefficients for successive deflation.
	do j=m,1,-1
	! Loop over each root to be found.
	! Start at zero to favor convergence to smallest remaining root.
		x=cmplx(0.,0.)
		call laguer(ad,j,x,its) !Find the root.
		if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
		roots(j)=x
		b=ad(j+1)
		do jj=j,1,-1
			c=ad(jj)
			ad(jj)=b
			b=x*b+c
		enddo
	enddo

	if (polish) then
		do j=1,m
	  !Polish the roots using the undeflated coefficients.
			call laguer(a,m,roots(j),its)
		enddo
	endif
  !adapted for SIGMA by CL - only root with positive real part is kept
	do j=1,m
		if (dreal(roots(j)).gt.0.0_dp) roots(1)=roots(j)
	enddo
	return
	END

SUBROUTINE laguer(a,m,x,its)
	use Tools
	INTEGER m,its,MAXIT,MR,MT
	REAL (KIND=dp) EPSS
	COMPLEX (KIND=dp) a(m+1),x
	PARAMETER (EPSS=1.e-15_dp,MR=8,MT=10,MAXIT=MT*MR)
	! Given the degree m and the complex coefficients a(1:m+1) of the polynomial
	! and given a complex value x, this routine improves x by Laguerre’s method until it con- verges, within the achievable roundoff limit, to a root of the given polynomial. The number of iterations taken is returned as its.
	! Parameters: EPSS is the estimated fractional roundoff error. We try to break (rare) limit cycles with MR different fractional values, once every MT steps, for MAXIT total allowed iterations.
	INTEGER iter,j
	REAL (KIND=dp) abx,abp,abm,err,frac(MR)
	COMPLEX (KIND=dp) dx,x1,b,d,f,g,h,sq,gp,gm,g2
	SAVE frac
	DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
	do iter=1,MAXIT
		its=iter
		b=a(m+1)
		err=abs(b)
		d=cmplx(0.,0.)
		f=cmplx(0.,0.)
		abx=abs(x)
		do j=m,1,-1
			f=x*f+d
			d=x*d+b
			b=x*b+a(j)
			err=abs(b)+abx*err
		enddo
		err=EPSS*err
		if(abs(b).le.err) then
			return
		else
	! Fractions used to break a limit cycle. Loop over iterations up to allowed maximum.
	! Efficient computation of the polynomial and its first two derivatives.
	! Estimate of roundoff error in evaluating polynomial. We are on the root.
	! The generic case: use Laguerre’s formula.
		g=d/b
		g2=g*g
		h=g2-2.*f/b
		sq=sqrt((m-1)*(m*h-g2))
		gp=g+sq
		gm=g-sq
		abp=abs(gp)
		abm=abs(gm)
		if(abp.lt.abm) gp=gm
		if (max(abp,abm).gt.0.) then
			dx=m/gp
		else
			dx=exp(cmplx(log(1.+abx),float(iter)))
		endif
	endif
	x1=x-dx
	if(x.eq.x1)return
	if (mod(iter,MT).ne.0) then
		x=x1
	else
		x=x-dx*frac(iter/MT)
	endif
	enddo
	!pause "too many iterations in laguer"
	END


	subroutine Blender(e1in,e2in,abun,nm,e1out,e2out)
	use Tools
	IMPLICIT NONE
	integer nm,j,iter
	REAL (KIND=dp) e1in(nm),e2in(nm),e1out,e2out,abun(nm)
	complex (KIND=dp) mm,m(nm),me,sum

	mm=dcmplx(1d0,0d0)
	do j=1,nm
		m(j)=dcmplx(e1in(j),e2in(j))
	enddo
	do iter=1,100
		sum=0d0
		do j=1,nm
! polarizability for homogeneous spheres
			sum=sum+((m(j)**2-mm**2)/(m(j)**2+2d0*mm**2))*abun(j)
! test for CDE polarizability - non conclusive yet
!      sum=sum+(m(j)**2*zlog(m(j)**2)/(m(j)**2-mm**2)-mm**2)*abun(j)*2d0/3d0
		enddo
		me=(2d0*sum+1d0)/(1d0-sum)
		me=mm*cdsqrt(me)
		mm=me
	enddo


	e1out=dreal(me)
	e2out=dimag(me)


! LLL mixing rule (not preferred)
!	me=0d0
!	do j=1,nm
!		me=me+m(j)**(2d0/3d0)*abun(j)
!	enddo
!	me=me**(3d0/2d0)
!	e1out=dreal(me)
!	e2out=dimag(me)

	return
	end

	subroutine maxgarn_2compo(e1in,e2in,e1icein,e2icein,abun,e1out,e2out)
	use Tools
	IMPLICIT NONE
	integer nm,j,iter
	real (KIND=dp) e1in,e2in,e1out,e2out,e1icein,e2icein,abun !up, down
	complex (KIND=dp) m1,m2,me,sum

	!Maxwell Garnet Rules is used
	m2=dcmplx(e1icein,e2icein) !mantle = coated material = matrix
	m1=dcmplx(e1in,e2in) !inner core consider as inclusions
	me = m2**2*( (2d0*m2**2+m1**2-2d0*(1d0-abun)*(m2**2-m1**2)) /(2d0*m2**2+m1**2+(1d0-abun)*(m2**2-m1**2) ))
  !from Mukai 1986: optical constants of the mixture of ices

	me=cdsqrt(me)
	e1out=dreal(me)
	e2out=dimag(me)

	return
	end



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      SUBROUTINE gauleg2(x1,x2,x,w,n)
      use Tools
      INTEGER n
      REAL (KIND=dp) x1,x2,x(n),w(n)
      REAL (KIND=dp) EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      REAL (KIND=dp) p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(pi*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

! ---------------------------------------------------------------------
      function cdlog10(x)
	IMPLICIT NONE
	complex*16 x, cdlog10

	cdlog10=log(x)/log(10d0)

	return
      end function cdlog10

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
