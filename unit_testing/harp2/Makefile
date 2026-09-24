# Builds and runs the HARP2 simulator unit test (COSP v2.1.4cesm layout).
#   make            build test_harp2_simulator
#   make test       build and run it
F90      ?= gfortran
F90FLAGS ?= -O2 -ffree-line-length-none -fcheck=bounds -finit-real=nan

ROOT = ../..
SRC  = $(ROOT)/model-interface/cosp_kinds.F90 \
       $(ROOT)/model-interface/cosp_errorHandling.F90 \
       $(ROOT)/src/cosp_constants.F90 \
       $(ROOT)/src/cosp_config.F90 \
       $(ROOT)/src/cosp_stats.F90 \
       $(ROOT)/src/simulator/HARP2_simulator/harp2_simulator.F90 \
       $(ROOT)/src/simulator/cosp_harp2_interface.F90 \
       test_harp2_simulator.F90
LUT  = $(ROOT)/src/simulator/HARP2_simulator/harp2_lut_670nm.txt

test_harp2_simulator: $(SRC)
	$(F90) $(F90FLAGS) $(SRC) -o $@

test: test_harp2_simulator
	./test_harp2_simulator $(LUT)

clean:
	rm -f *.o *.mod test_harp2_simulator
