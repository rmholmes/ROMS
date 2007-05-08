# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2007 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
# Include file for Solaris F95 compiler on SUN
# -------------------------------------------------------------------------
#
# ARPACK_LIBDIR  ARPACK libary directory
# FC             Name of the fortran compiler to use
# FFLAGS         Flags to the fortran compiler
# CPP            Name of the C-preprocessor
# CPPFLAGS       Flags to the C-preprocessor
# CLEAN          Name of cleaning executable after C-preprocessing
# NETCDF_INCDIR  NetCDF include directory
# NETCDF_LIBDIR  NetCDF libary directory
# LD             Program to load the objects into an executable
# LDFLAGS        Flags to the loader
# RANLIB         Name of ranlib command
# MDEPFLAGS      Flags for sfmakedepend  (-s if you keep .f files)
#
# First the defaults
#
               FC := f95
           FFLAGS := -u -U
              CPP := /usr/lib/cpp
         CPPFLAGS := -P
               LD := $(FC)
          LDFLAGS :=
               AR := ar
          ARFLAGS := r
            MKDIR := mkdir -p
               RM := rm -f
           RANLIB := ranlib
	     PERL := perl
             TEST := test

        MDEPFLAGS := --cpp --fext=f90 --file=- --objdir=$(SCRATCH_DIR)

#
# Library locations, can be overridden by environment variables.
#

ifdef USE_LARGE
           FFLAGS += -xarch=v9
       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib64
else
       MCT_LIBDIR ?= /usr/local/mct/lib
    NETCDF_INCDIR ?= /usr/local/include
    NETCDF_LIBDIR ?= /usr/local/lib
endif
         CPPFLAGS += -I$(NETCDF_INCDIR)
             LIBS := -L$(NETCDF_LIBDIR) -lnetcdf

ifdef USE_ARPACK
 ifdef USE_MPI
   PARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(PARPACK_LIBDIR) -lparpack
 endif
    ARPACK_LIBDIR ?= /usr/local/lib
             LIBS += -L$(ARPACK_LIBDIR) -larpack
endif

ifdef USE_MPI
         CPPFLAGS += -DMPI
               FC := tmf90
               LD := $(FC)
             LIBS += -lmpi
endif

ifdef USE_OpenMP
         CPPFLAGS += -D_OPENMP
           FFLAGS += -openmp
endif

ifdef USE_DEBUG
           FFLAGS += -g -C
else
           FFLAGS += -fast
endif

             LIBS += -lnsl

ifdef SWAN_COUPLE
           FFLAGS += -fixed -I/usr/local/mct/include
             LIBS += -L$(MCT_LIBDIR) -lmct -lmpeu
endif

#
# Set free form format in source files to allow long string for
# local directory and compilation flags inside the code.
#

$(SCRATCH_DIR)/mod_ncparam.o: FFLAGS += -free
$(SCRATCH_DIR)/mod_strings.o: FFLAGS += -free

#
# Supress free format in SWAN source files since there are comments
# beyond column 72.
#

ifdef SWAN_COUPLE

$(SCRATCH_DIR)/ocpcre.o: FFLAGS += -fixed
$(SCRATCH_DIR)/ocpids.o: FFLAGS += -fixed
$(SCRATCH_DIR)/ocpmix.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swancom1.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swancom2.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swancom3.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swancom4.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swancom5.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanmain.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanout1.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanout2.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanparll.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanpre1.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanpre2.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swanser.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swmod1.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swmod2.o: FFLAGS += -fixed
$(SCRATCH_DIR)/swmod3.o: FFLAGS += -fixed

endif
