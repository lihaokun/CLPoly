# B2B C++ driver build
# 使用方式：make -C proof/b2b cpp_driver
#
# 注意：编译 with -DNDEBUG 关 C++ assert（防止 abort 整个 driver）

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -DNDEBUG -Wall -Wextra
GMP_CFLAGS ?= $(shell pkg-config --cflags gmp 2>/dev/null)
GMP_LIBS   ?= $(shell pkg-config --libs gmpxx 2>/dev/null)
PACKAGE_PREFIX ?= $(shell brew --prefix 2>/dev/null)
PACKAGE_CFLAGS = $(if $(strip $(PACKAGE_PREFIX)),-I$(PACKAGE_PREFIX)/include,)
INCLUDES  = -I third_party -I ../../ $(GMP_CFLAGS) $(PACKAGE_CFLAGS)

DRIVER_SRCS = driver/cpp_driver.cc types/b2b_types.cc registry/cpp_registry.cc
DRIVER_OBJS = $(DRIVER_SRCS:.cc=.o)

# The registry uses header-defined factorization routines, whose only
# out-of-line CLPoly dependencies are the univariate polynomial and variable
# implementations.  Link those sources directly so the B2B gate does not
# depend on first producing a platform-specific full-library shared object.
CLPOLY_SRCS = ../../clpoly/upolynomial.cc ../../clpoly/variable.cc

TEST_SRCS = types/test_types.cc types/b2b_types.cc
TEST_OBJS = $(TEST_SRCS:.cc=.o)

LIBS       = $(if $(strip $(GMP_LIBS)),$(GMP_LIBS),-lgmpxx -lgmp)

.PHONY: all clean test_types

all: cpp_driver

cpp_driver: $(DRIVER_OBJS) $(CLPOLY_SRCS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(DRIVER_OBJS) $(CLPOLY_SRCS) $(LIBS) -o $@

test_types: $(TEST_OBJS) $(CLPOLY_SRCS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(INCLUDES) $(TEST_OBJS) $(CLPOLY_SRCS) $(LIBS) -o $@
	./$@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f cpp_driver test_types driver/*.o types/*.o registry/*.o
