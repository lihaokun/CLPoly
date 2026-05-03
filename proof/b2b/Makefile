# B2B C++ driver build
# 使用方式：make -C proof/b2b cpp_driver
#
# 注意：编译 with -DNDEBUG 关 C++ assert（防止 abort 整个 driver）

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -DNDEBUG -Wall -Wextra
INCLUDES  = -I third_party -I ../../

DRIVER_SRCS = driver/cpp_driver.cc types/b2b_types.cc registry/cpp_registry.cc
DRIVER_OBJS = $(DRIVER_SRCS:.cc=.o)

TEST_SRCS = types/test_types.cc types/b2b_types.cc
TEST_OBJS = $(TEST_SRCS:.cc=.o)

CLPOLY_LIB = ../../lib/libclpoly.so
LIBS       = -L../../lib -lclpoly -lgmp -Wl,-rpath,$(abspath ../../lib)

.PHONY: all clean test_types

all: cpp_driver

cpp_driver: $(DRIVER_OBJS) $(CLPOLY_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DRIVER_OBJS) $(LIBS) -o $@

test_types: $(TEST_OBJS) $(CLPOLY_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEST_OBJS) $(LIBS) -o $@
	./$@

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f cpp_driver test_types driver/*.o types/*.o registry/*.o
