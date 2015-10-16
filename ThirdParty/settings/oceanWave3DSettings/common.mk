# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D
LIBNAME  = libOceanWave3D.so

# Installation directory
#INSTALLDIR = $(PWD)/../bin
INSTALLDIR = $(WAVES_APPBIN)
LIBINSTALLDIR = $(WAVES_LIBBIN)

# Build directory where object files are stored 
BUILDDIR = $(PWD)/../build


FC       = gfortran
LIBDIRS  = -L$(PWD)/../lib 
LINLIB   = -ltmglib_gfortran -llapack_gfortran  -lskit_gfortran -lblas
DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
OPTFLAGS = -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
SHLIBFLAGS  = -shared -O2 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all


