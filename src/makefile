CXX=g++
IPATHS=-I./
CFLAGS=-O3 -flto 
objects= variable.o variable_test.o variable_test atomic_polynomial_test.o
atomic_polynomial_test:atomic_polynomial_test.o variable.o 
	$(CXX) $(CFLAGS)   variable.o atomic_polynomial_test.o  -o atomic_polynomial_test
atomic_polynomial_test.o:test/atomic_polynomial_test.cc atomic_polynomial.hh basic.hh
	$(CXX) $(CFLAGS) $(IPATHS)  test/atomic_polynomial_test.cc  -c
variable_test:variable_test.o variable.o
	$(CXX) $(CFLAGS)   variable.o variable_test.o  -o variable_test
variable.o: variable.cc variable.hh
	$(CXX) $(CFLAGS)  $(IPATHS) variable.cc  -c
variable_test.o:test/variable_test.cc variable.hh
	$(CXX) $(CFLAGS)  $(IPATHS) test/variable_test.cc  -c
.PHONY : clean
clean :
	-rm $(objects)