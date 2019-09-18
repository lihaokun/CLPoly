CXX=g++ 
IPATHS=-I./ 
CFLAGS=-O3 -flto  
objects= variable.o variable_test.o variable_test atomic_polynomial_test.o atomic_polynomial_test monomial.o monomial_test.o monomial_test polynomial_test.o polynomial_test
polynomial_test:variable.o monomial.o polynomial_test.o
	$(CXX) $(CFLAGS) variable.o monomial.o polynomial_test.o -o polynomial_test -lgmpxx -lgmp
polynomial_test.o:test/polynomial_test.cc polynomial.hh monomial.hh basic.hh variable.hh  atomic_polynomial.hh number.hh
	$(CXX) $(CFLAGS) $(IPATHS) test/polynomial_test.cc -c
monomial_test:variable.o monomial.o monomial_test.o
	$(CXX) $(CFLAGS) variable.o monomial.o monomial_test.o -o monomial_test
monomial_test.o:test/monomial_test.cc monomial.hh basic.hh variable.hh
	$(CXX) $(CFLAGS) $(IPATHS) test/monomial_test.cc -c
monomial.o:monomial.cc monomial.hh basic.hh variable.hh
	$(CXX) $(CFLAGS) $(IPATHS) monomial.cc -c
atomic_polynomial_test:atomic_polynomial_test.o variable.o  monomial.o
	$(CXX) $(CFLAGS)   monomial.o variable.o atomic_polynomial_test.o  -o atomic_polynomial_test -lgmpxx -lgmp
atomic_polynomial_test.o:test/atomic_polynomial_test.cc atomic_polynomial.hh basic.hh monomial.hh  variable.hh number.hh
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