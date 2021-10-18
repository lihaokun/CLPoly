.PHONY: clean
all:
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

clpoly_hh=$(wildcard clpoly/*.hh)
CXX=g++ 
IPATHS=-I./ 
CFLAGS=-O3 -flto -DNDEBUG  

Numberlib= -lgmpxx -lgmp  
clpoly_hh=$(wildcard clpoly/*.hh)
clpoly_cc=$(wildcard clpoly/*.cc)
clpoly_d_o=$(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/debug/clpoly/%.o)
clpoly_r_o=$(clpoly_cc:clpoly/%.cc=$(CLPoly_BUILD_DIR)/release/clpoly/%.o)
$(CLPoly_BUILD_DIR)/release/%.o:%.cc $(clpoly_hh)
	mkdir -p $(CLPoly_BUILD_DIR)/release/clpoly
	$(CXX) $(CLPoly_REL) $(IPATHS) -c $< -o $@ $(Numberlib)
$(CLPoly_BUILD_DIR)/debug/%.o:%.cc $(clpoly_hh)
	mkdir -p $(CLPoly_BUILD_DIR)/debug/clpoly
	$(CXX) $(CLPoly_DEB) $(IPATHS) -c $< -o $@ $(Numberlib)
%/lib/debug/libclpoly.a:$(clpoly_d_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_d_o)
	
%/lib/release/libclpoly.a:$(clpoly_r_o)
	mkdir -p $(dir $@)
	ar -rsv $@ $(clpoly_r_o)

$(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a:$(CLPoly_BUILD_DIR)/lib/debug/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/debug/clpoly
	cp $(CLPoly_BUILD_DIR)/lib/debug/libclpoly.a $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a

$(CLPoly_LIB_DIR)/release/clpoly/libclpoly.a:$(CLPoly_BUILD_DIR)/lib/release/libclpoly.a
	mkdir -p $(CLPoly_LIB_DIR)/release/clpoly
	cp $(CLPoly_BUILD_DIR)/lib/release/libclpoly.a $(CLPoly_LIB_DIR)/release/clpoly/libclpoly.a
%:%.cc $(CLPoly_LIB_DIR)/debug/clpoly/libclpoly.a $(clpoly_hh)
	$(CXX) $(CLPoly_DEB) $(IPATHS)  $< -o $@ $(Numberlib)  -L$(CLPoly_LIB_DIR)/debug/clpoly -lclpoly

clean:
	rm -r $(clpoly_o) $(CLPoly_BUILD_DIR)/lib/libclpoly.a $(CLPoly_LIB_DIR)/clpoly/libclpoly.a