AUTOMAKE_OPTIONS = foreign subdir-objects

SOURCES = src/main.F90 src/schemes/schemes_module.F90 src/schemes/abbdf3_module.F90 src/schemes/rkcn2_module.F90 src/comm_module.F90 src/transform/transform_module.F90 src/params_module.F90 src/fields_module.F90 src/IO/output_module.F90 src/IO/input_module.F90 src/grid/grid_module.F90 src/grid/grid_ops_module.F90 src/schemes/equations_module.F90
if MPI
XSOURCES = src/transform/transpose_pencil_module.F90
else
XSOURCES = 
endif

bin_PROGRAMS = CaTSM
CaTSM_SOURCES = $(SOURCES) $(XSOURCES)

AM_LDFLAGS = $(LIBFFTW3) $(LIBNC)

src/main.o: src/schemes/schemes_module.o src/IO/output_module.o src/IO/input_module.o
src/params_module.o: src/comm_module.o src/grid/grid_module.o
src/fields_module.o: src/grid/grid_ops_module.o src/grid/grid_module.o

src/schemes/schemes_module.o: src/schemes/abbdf3_module.o src/schemes/rkcn2_module.o src/grid/grid_module.o src/transform/transform_module.o src/params_module.o src/fields_module.o src/comm_module.o
src/schemes/abbdf3_module.o: src/grid/grid_module.o src/comm_module.o src/params_module.o src/fields_module.o src/schemes/equations_module.o
src/schemes/rkcn2_module.o: src/grid/grid_module.o src/comm_module.o src/params_module.o src/fields_module.o src/schemes/equations_module.o
src/schemes/equations_module.o: src/grid/grid_module.o src/params_module.o src/fields_module.o

src/grid/grid_ops_module.o: src/transform/transform_module.o src/grid/grid_module.o
src/grid/grid_module.o: src/transform/transform_module.o src/comm_module.o

src/IO/output_module.o: src/params_module.o src/comm_module.o src/fields_module.o src/IO/output_netcdf.F90 src/grid/grid_ops_module.o
src/IO/input_module.o: src/fields_module.o src/comm_module.o src/IO/output_module.o src/IO/input_netcdf.F90

src/transform/transform_module.o: src/transform/transform_serial.F90

clean-local: 
	rm -f *.mod src/*.mod
