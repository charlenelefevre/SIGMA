# makefile for SIGMA (with comments!)
# MacOSX 10.10.1 gfortran gcc version 5.0.0 20141005
#

# compiler= FC, flags = FFlags
# linker= LINKER, flags= LDFLAGS, libraries=LIBS
FC	  = gfortran
LINKER	  = gfortran

ifeq ($(ifort),true)
  FC	  = ifort
  LINKER  = ifort
endif

# array boundary check
ifeq ($(debug),true)
  ifeq ($(ifort),true)
    DEBUGGING = -check all -traceback -check bounds -O0 -g -check -fpe1
  else	
    DEBUGGING = -fbounds-check -fbacktrace
  endif
endif

# Platform specific compilation options
ifeq ($(ifort),true)
  FLAG_ALL      = -O3 -g -extend-source -zero -prec-div $(DEBUGGING)
  FLAG_LINUX    = -xHOST -fpp
  FLAG_MAC      = -xHOST -opt-prefetch -static-intel -fpp
else
  FLAG_ALL      = -O3 -g -fdefault-double-8 -fdefault-real-8 $(DEBUGGING)
  FLAG_LINUX    = -cpp
  FLAG_MAC      = -m64 -cpp
endif

ifeq ($(fitsio),true)
  FLAG_FITS		= -DUSE_FITSIO
  LIBS_FITS		= -lcfitsio -L/home/lefevre/miniconda2/lib/
endif

  FFLAGS  = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
  LDFLAGS = $(FLAG_ALL) $(FLAG_LINUX) $(FLAG_FITS)
  LIBS    = $(LIBS_FITS)


# files to make
OBJS	= SIGMA.o \
	  SIGMA_extra.o \
	  RefInd.o \
	  dmilay_f95.o \
	  RefIndData.o \
          grain_dist.o

# program name and install location
PROGRAM       = SIGMA
DEST	      = ${HOME}/bin

# make actions 
all:		$(PROGRAM)
clean:;		rm -f $(OBJS) $(PROGRAM) *.mod *.i
install:	$(PROGRAM)
		mv $(PROGRAM) $(DEST)

# how to compile program 
.SUFFIXES : .o .f .f90 .F

.f.o:
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

.F.o:
	$(FC) $(LDFLAGS) -c $<

$(PROGRAM):     $(OBJS)
		$(LINKER) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

# recompile everything if SIGMA.f90 has changed 
$(OBJS):	SIGMA.f90


