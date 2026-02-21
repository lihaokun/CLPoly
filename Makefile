.PHONY: all clean
all: libs

## Load Previous Configuration ####################################################################

-include config.mk

## Configurable options ###########################################################################

# Top-level build directory (internal artifacts)
CLPoly_BUILD_DIR      ?= _build
# Public library output directory (stable external link path: -Llib -lclpoly)
CLPoly_LIB_DIR        ?= lib

# Sets of compile flags for different build types
CLPoly_REL    ?= -O3 -D NDEBUG
CLPoly_DEB    ?= -g -D DEBUG
CLPoly_FPIC   ?= -fpic

# GNU Standard Install Prefix
CLPoly_prefix ?= /usr/local

config:
	@( echo 'CLPoly_BUILD_DIR?=$(CLPoly_BUILD_DIR)'   ; \
	   echo 'CLPoly_LIB_DIR?=$(CLPoly_LIB_DIR)'       ; \
	   echo 'CLPoly_REL?=$(CLPoly_REL)'               ; \
	   echo 'CLPoly_DEB?=$(CLPoly_DEB)'               ; \
	   echo 'CLPoly_FPIC?=$(CLPoly_FPIC)'             ; \
	   echo 'CLPoly_prefix?=$(CLPoly_prefix)'          ) > config.mk

CXX      = g++
IPATHS   = -I./
DEPFLAGS = -MMD -MP

Numberlib  = -lgmpxx -lgmp
FLINT_LIBS = -lflint
NTL_LIBS   = -lntl -lm -lpthread

## Internal path shorthands #######################################################################

B       = $(CLPoly_BUILD_DIR)
OBJ_DEB = $(B)/debug/obj
OBJ_REL = $(B)/release/obj
BIN_DEB = $(B)/debug/bin
BIN_REL = $(B)/release/bin

## Library source files (excluding deprecated interval.cc) ########################################

clpoly_cc  = $(filter-out clpoly/interval.cc,$(wildcard clpoly/*.cc))
clpoly_d_o = $(clpoly_cc:clpoly/%.cc=$(OBJ_DEB)/clpoly/%.o)
clpoly_r_o = $(clpoly_cc:clpoly/%.cc=$(OBJ_REL)/clpoly/%.o)

## Library object compilation #####################################################################

$(OBJ_REL)/clpoly/%.o: clpoly/%.cc
	mkdir -p $(dir $@)
	$(CXX) $(CLPoly_REL) $(CLPoly_FPIC) $(DEPFLAGS) $(IPATHS) -c $< -o $@

$(OBJ_DEB)/clpoly/%.o: clpoly/%.cc
	mkdir -p $(dir $@)
	$(CXX) $(CLPoly_DEB) $(CLPoly_FPIC) $(DEPFLAGS) $(IPATHS) -c $< -o $@

# Prevent Make from deleting .o as intermediate files
.PRECIOUS: $(OBJ_REL)/clpoly/%.o $(OBJ_DEB)/clpoly/%.o

## Static libraries → lib/ #######################################################################

$(CLPoly_LIB_DIR)/libclpoly.a: $(clpoly_r_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $^

$(CLPoly_LIB_DIR)/debug/libclpoly.a: $(clpoly_d_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $^

## Shared libraries → lib/ #######################################################################

$(CLPoly_LIB_DIR)/libclpoly.so: $(clpoly_r_o)
	mkdir -p $(dir $@)
	$(CXX) -shared -o $@ $^ $(Numberlib)

$(CLPoly_LIB_DIR)/debug/libclpoly.so: $(clpoly_d_o)
	mkdir -p $(dir $@)
	$(CXX) -shared -o $@ $^ $(Numberlib)

## Aggregate library target ######################################################################

.PHONY: libs
libs: $(CLPoly_LIB_DIR)/libclpoly.a \
      $(CLPoly_LIB_DIR)/libclpoly.so \
      $(CLPoly_LIB_DIR)/debug/libclpoly.a \
      $(CLPoly_LIB_DIR)/debug/libclpoly.so

## Test binaries (debug build) → _build/debug/bin/ ###############################################
## `make test/test_name` builds test/test_name.cc → _build/debug/bin/test_name

test/%: test/%.cc $(CLPoly_LIB_DIR)/debug/libclpoly.a
	mkdir -p $(BIN_DEB)
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $(BIN_DEB)/$* $(CLPoly_LIB_DIR)/debug/libclpoly.a $(Numberlib)

## Cross-library correctness tests (explicit rules take priority) ################################

$(BIN_DEB)/test_crosscheck_flint: test/test_crosscheck_flint.cc $(CLPoly_LIB_DIR)/debug/libclpoly.a
	mkdir -p $(BIN_DEB)
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $@ $(CLPoly_LIB_DIR)/debug/libclpoly.a $(Numberlib) $(FLINT_LIBS)

$(BIN_DEB)/test_crosscheck_ntl: test/test_crosscheck_ntl.cc $(CLPoly_LIB_DIR)/debug/libclpoly.a
	mkdir -p $(BIN_DEB)
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $@ $(CLPoly_LIB_DIR)/debug/libclpoly.a $(Numberlib) $(NTL_LIBS)

.PHONY: crosscheck
crosscheck: $(BIN_DEB)/test_crosscheck_flint $(BIN_DEB)/test_crosscheck_ntl
	$(BIN_DEB)/test_crosscheck_flint
	$(BIN_DEB)/test_crosscheck_ntl

## Performance benchmarks (release build) → _build/release/bin/ ##################################

$(BIN_REL)/bench_clpoly: test/bench_clpoly.cc $(CLPoly_LIB_DIR)/libclpoly.a
	mkdir -p $(BIN_REL)
	$(CXX) $(CLPoly_REL) $(DEPFLAGS) $(IPATHS) $< -o $@ $(CLPoly_LIB_DIR)/libclpoly.a $(Numberlib)

$(BIN_REL)/bench_comparative: test/bench_comparative.cc $(CLPoly_LIB_DIR)/libclpoly.a
	mkdir -p $(BIN_REL)
	$(CXX) $(CLPoly_REL) $(DEPFLAGS) $(IPATHS) $< -o $@ $(CLPoly_LIB_DIR)/libclpoly.a $(Numberlib) $(FLINT_LIBS) $(NTL_LIBS)

$(BIN_REL)/test_stress_factorize: test/test_stress_factorize.cc $(CLPoly_LIB_DIR)/libclpoly.a
	mkdir -p $(BIN_REL)
	$(CXX) $(CLPoly_REL) $(DEPFLAGS) $(IPATHS) $< -o $@ $(CLPoly_LIB_DIR)/libclpoly.a $(Numberlib)

.PHONY: stress
stress: $(BIN_REL)/test_stress_factorize
	$(BIN_REL)/test_stress_factorize

.PHONY: bench bench-clpoly bench-comparative
bench-clpoly: $(BIN_REL)/bench_clpoly
	$(BIN_REL)/bench_clpoly

bench-comparative: $(BIN_REL)/bench_comparative
	$(BIN_REL)/bench_comparative

bench: bench-clpoly bench-comparative

## Auto-generated dependencies ####################################################################

-include $(clpoly_d_o:.o=.d)
-include $(clpoly_r_o:.o=.d)
-include $(wildcard $(BIN_DEB)/*.d)
-include $(wildcard $(BIN_REL)/*.d)

## Clean ##########################################################################################

clean:
	rm -rf $(CLPoly_BUILD_DIR)
	rm -f $(CLPoly_LIB_DIR)/libclpoly.a $(CLPoly_LIB_DIR)/libclpoly.so
	rm -rf $(CLPoly_LIB_DIR)/debug
	@# Legacy cleanup (old build layout)
	rm -rf build/
	rm -f test/*.d
