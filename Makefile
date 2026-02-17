.PHONY: all clean
all: make_lib_a
## Load Previous Configuration ####################################################################

-include config.mk

## Configurable options ###########################################################################

# Directory to store object files, libraries, executables, and dependencies:
CLPoly_BUILD_DIR      ?= build
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

## Library source files (excluding deprecated interval.cc) ########################################

clpoly_cc  = $(filter-out clpoly/interval.cc,$(wildcard clpoly/*.cc))
clpoly_d_o = $(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/debug/clpoly/%.o)
clpoly_r_o = $(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/release/clpoly/%.o)

## Library object compilation #####################################################################

$(CLPoly_BUILD_DIR)/release/%.o: %.cc
	mkdir -p $(CLPoly_BUILD_DIR)/release/clpoly
	$(CXX) $(CLPoly_REL) $(CLPoly_FPIC) $(DEPFLAGS) $(IPATHS) -c $< -o $@

$(CLPoly_BUILD_DIR)/debug/%.o: %.cc
	mkdir -p $(CLPoly_BUILD_DIR)/debug/clpoly
	$(CXX) $(CLPoly_DEB) $(CLPoly_FPIC) $(DEPFLAGS) $(IPATHS) -c $< -o $@

# Prevent Make from deleting .o as intermediate files
.PRECIOUS: $(CLPoly_BUILD_DIR)/release/%.o $(CLPoly_BUILD_DIR)/debug/%.o

## Static library archiving #######################################################################

%/lib/release/libclpoly.a: $(clpoly_r_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_r_o)

%/lib/debug/libclpoly.a: $(clpoly_d_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_d_o)

make_lib_a: $(CLPoly_LIB_DIR)/clpoly/libclpoly.a $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a

$(CLPoly_LIB_DIR)/clpoly/libclpoly.a: $(CLPoly_BUILD_DIR)/lib/release/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/clpoly
	cp $< $@

$(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a: $(CLPoly_BUILD_DIR)/lib/debug/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/debug/clpoly
	cp $< $@

## Generic test target (debug build) ##############################################################

%: %.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $@ $(Numberlib) -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

## Cross-library correctness tests ################################################################

test/test_crosscheck_flint: test/test_crosscheck_flint.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $@ $(Numberlib) $(FLINT_LIBS) -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

test/test_crosscheck_ntl: test/test_crosscheck_ntl.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a
	$(CXX) $(CLPoly_DEB) $(DEPFLAGS) $(IPATHS) $< -o $@ $(Numberlib) $(NTL_LIBS) -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

.PHONY: crosscheck
crosscheck: test/test_crosscheck_flint test/test_crosscheck_ntl
	./test/test_crosscheck_flint
	./test/test_crosscheck_ntl

## Performance benchmarks (release build) #########################################################

test/bench_clpoly: test/bench_clpoly.cc $(CLPoly_LIB_DIR)/clpoly/libclpoly.a
	$(CXX) $(CLPoly_REL) $(DEPFLAGS) $(IPATHS) $< -o $@ $(Numberlib) -L$(CLPoly_LIB_DIR)/clpoly -lclpoly

test/bench_comparative: test/bench_comparative.cc $(CLPoly_LIB_DIR)/clpoly/libclpoly.a
	$(CXX) $(CLPoly_REL) $(DEPFLAGS) $(IPATHS) $< -o $@ $(Numberlib) $(FLINT_LIBS) $(NTL_LIBS) -L$(CLPoly_LIB_DIR)/clpoly -lclpoly

.PHONY: bench bench-clpoly bench-comparative
bench-clpoly: test/bench_clpoly
	./test/bench_clpoly

bench-comparative: test/bench_comparative
	./test/bench_comparative

bench: bench-clpoly bench-comparative

## Auto-generated dependencies ####################################################################

-include $(clpoly_d_o:.o=.d)
-include $(clpoly_r_o:.o=.d)
-include $(wildcard test/*.d)

## Clean ##########################################################################################

clean:
	rm -rf $(CLPoly_BUILD_DIR)/debug $(CLPoly_BUILD_DIR)/release $(CLPoly_BUILD_DIR)/lib
	rm -rf $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(CLPoly_LIB_DIR)/clpoly/libclpoly.a
	rm -f test/*.d
