LOCAL_FC             = gfortran
LOCAL_DWDTOOLS_INC   = /newhome/hes134/apps/icontools/include
LOCAL_DWDTOOLS       = -L/newhome/hes134/apps/icontools/lib
LOCAL_NETCDF         = -L/usr/lib64
LOCAL_INCLUDE        = -I${LOCAL_DWDTOOLS_INC} -I/usr/lib64/gfortran/modules/
LOCAL_LIBS			 = -licontools -lcbasic_f -lcbasic -lcdi_f2003 -lcdi -lftnbasic_f
LOCAL_OMP            = -fopenmp
LOCAL_DEBUG          = -ffpe-trap=zero,overflow,invalid -finit-real=nan -finit-integer=-2147483648 -finit-character=127 -g
LOCAL_FOPT           = -march=native -O0 -g -ffast-math $(LOCAL_OMP) -fbacktrace
LOCAL_WARN           = -std=f2003 -fall-intrinsics -Wall

FFLAGS = ${LOCAL_FOPT} $(LOCAL_DEBUG) $(LOCAL_WARN) ${ADDITIONAL_FLAGS} \
	-cpp -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -DNOMPI -I /usr/lib64/gfortran/modules -I/local/apps/eccodes/2.19.0/include -fbacktrace -g -L/local/apps/eccodes/2.19.0/lib -leccodes_f90 -leccodes -lm -L /usr/lib64  -ljasper -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lsz -lgfortran

SUFFIXES = .o .f90

EXE = ex04_icon_point_query interpol_delaunay

all:  ${EXE}

ifeq ($(target),ddm06)
  $(info Machine: ddm06)
  FC = gcc
  FFLAGS = -cpp -I /usr/lib64/gfortran/modules -I/local/apps/eccodes/2.19.0/include -fbacktrace -g 
  LDFLAGS = $(FFLAGS) -L/local/apps/eccodes/2.19.0/lib -leccodes_f90 -leccodes -lm -L /usr/lib64  -ljasper -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lsz -lgfortran
endif

interpol_delaunay: interpol_delaunay.o
	@echo "04 Building \"interpol_delaunay\""
	$(LOCAL_FC)  -o interpol_delaunay interpol_delaunay.o $(FFLAGS) \
	-L../lib -L/usr/lib64 $(LOCAL_DWDTOOLS) $(LOCAL_NETCDF) \
	$(LOCAL_LIBS) -lnetcdf -lnetcdff -leccodes  -lhdf5 -ljasper

ex04_icon_point_query: ex04_icon_point_query.o
	@echo "04 Building \"icon_point_query\""
	$(LOCAL_FC)  -o ex04_icon_point_query ex04_icon_point_query.o $(FFLAGS) \
	-L../lib -L/usr/lib64 $(LOCAL_DWDTOOLS) $(LOCAL_NETCDF) \
	$(LOCAL_LIBS) -lnetcdf -lnetcdff -leccodes  -lhdf5 -ljasper

%.o:  %.f90
	$(LOCAL_FC) -c -o $*.o  $(FFLAGS) $(LOCAL_INCLUDE)  $<
 
clean:
	rm -f *.o ${EXE}
