.PHONY: clean
all:make_lib_a
## Load Previous Configuration ####################################################################

-include config.mk

## Configurable options ###########################################################################

# Directory to store object files, libraries, executables, and dependencies:
CLPoly_BUILD_DIR      ?= build
CLPoly_LIB_DIR   ?=lib
# Include debug-symbols in release builds
CLPoly_RELSYM ?= -g

# Sets of compile flags for different build types
CLPoly_REL    ?= -O3 -D NDEBUG
CLPoly_DEB    ?=  -g -D DEBUG 
CLPoly_PRF    ?= -O3 -D NDEBUG
CLPoly_FPIC   ?= -fpic

# GNU Standard Install Prefix
CLPoly_prefix         ?= /usr/local

config:
	@( echo 'CLPoly_BUILD_DIR?=$(CLPoly_BUILD_DIR)'           ; \
	   echo 'CLPoly_LIB_DIR?=$(CLPoly_LIB_DIR)'           ; \
	   echo 'CLPoly_RELSYM?=$(CLPoly_RELSYM)' ; \
	   echo 'CLPoly_REL?=$(CLPoly_REL)'       ; \
	   echo 'CLPoly_DEB?=$(CLPoly_DEB)'       ; \
	   echo 'CLPoly_PRF?=$(CLPoly_PRF)'       ; \
	   echo 'CLPoly_FPIC?=$(CLPoly_FPIC)'     ; \
	   echo 'CLPoly_prefix?=$(CLPoly_prefix)'                 ) > config.mk

CXX=g++ 
IPATHS=-I./ 
CFLAGS=-O3 -flto -DNDEBUG  

Numberlib= -lgmpxx  -lgmp 
interval_lib=-lmpria 
clpoly_interval_hh=$(wildcard clpoly/*.hh)
interval_hh=clpoly/interval.hh
clpoly_hh=$(filter-out $(interval_hh),$(clpoly_interval_hh))
clpoly_interval_cc=$(wildcard clpoly/*.cc)
interval_cc=clpoly/interval.cc
clpoly_cc=$(filter-out $(interval_cc),$(clpoly_interval_cc))
clpoly_d_o=$(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/debug/clpoly/%.o)
clpoly_r_o=$(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/release/clpoly/%.o)
$(CLPoly_BUILD_DIR)/release/%.o:%.cc $(clpoly_hh)
	mkdir -p $(CLPoly_BUILD_DIR)/release/clpoly
	$(CXX) $(CLPoly_REL) $(CLPoly_FPIC) $(IPATHS) -c $< -o $@ 

#$(Numberlib)
$(CLPoly_BUILD_DIR)/debug/%.o:%.cc $(clpoly_hh)
	mkdir -p $(CLPoly_BUILD_DIR)/debug/clpoly
	$(CXX) $(CLPoly_DEB) $(CLPoly_FPIC) $(IPATHS) -c $< -o $@
#$(Numberlib)
%/lib/debug/libclpoly.a:$(clpoly_d_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_d_o)
	
%/lib/release/libclpoly.a:$(clpoly_r_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_r_o)

make_lib_a:$(CLPoly_LIB_DIR)/clpoly/libclpoly.a  $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a

$(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a:$(CLPoly_BUILD_DIR)/lib/debug/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/debug/clpoly
	cp $(CLPoly_BUILD_DIR)/lib/debug/libclpoly.a $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a

$(CLPoly_LIB_DIR)/clpoly/libclpoly.a:$(CLPoly_BUILD_DIR)/lib/release/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/clpoly
	cp $(CLPoly_BUILD_DIR)/lib/release/libclpoly.a $(CLPoly_LIB_DIR)/clpoly/libclpoly.a

# ---- Optional: cross-library correctness tests ----
FLINT_LIBS = -lflint
NTL_LIBS   = -lntl -lm -lpthread

test/test_crosscheck_flint: test/test_crosscheck_flint.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(clpoly_hh)
	$(CXX) $(CLPoly_DEB) $(IPATHS) $< -o $@ $(Numberlib) $(FLINT_LIBS) -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

test/test_crosscheck_ntl: test/test_crosscheck_ntl.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(clpoly_hh)
	$(CXX) $(CLPoly_DEB) $(IPATHS) $< -o $@ $(Numberlib) $(NTL_LIBS) -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

.PHONY: crosscheck
crosscheck: test/test_crosscheck_flint test/test_crosscheck_ntl
	./test/test_crosscheck_flint
	./test/test_crosscheck_ntl

%:%.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(clpoly_hh)
	$(CXX) $(CLPoly_DEB) $(IPATHS)  $< -o $@ $(Numberlib)  -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

clean:
	rm -r $(clpoly_d_o) $(clpoly_r_o) $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(CLPoly_LIB_DIR)/clpoly/libclpoly.a