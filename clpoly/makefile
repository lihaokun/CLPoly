CXX=g++ 
IPATHS=-I./ 
CFLAGS=-O3 -flto -DNDEBUG  
Numberlib= -lgmpxx -lgmp  
objects= variable_test basic_polynomial_test   monomial_test polynomial_test
polynomial_test:test/polynomial_test.cc polynomial.hh monomial.hh basic_monomial.hh  basic.hh variable.hh  basic_polynomial.hh number.hh
	$(CXX) $(CFLAGS) $(IPATHS) test/polynomial_test.cc -o polynomial_test $(Numberlib)
# polynomial_test.o:test/polynomial_test.cc polynomial.hh monomial.hh basic_monomial.hh  basic.hh variable.hh  basic_polynomial.hh number.hh
# 	$(CXX) $(CFLAGS) $(IPATHS) test/polynomial_test.cc -c
monomial_test:test/monomial_test.cc monomial.hh basic_monomial.hh  basic.hh variable.hh
	$(CXX) $(CFLAGS)  $(IPATHS) test/monomial_test.cc -o monomial_test
# monomial_test.o:test/monomial_test.cc monomial.hh basic_monomial.hh  basic.hh variable.hh
# 	$(CXX) $(CFLAGS) $(IPATHS) test/monomial_test.cc -c
basic_polynomial_test:test/basic_polynomial_test.cc basic_polynomial.hh basic.hh monomial.hh basic_monomial.hh   variable.hh number.hh 
	$(CXX) $(CFLAGS)  $(IPATHS)  test/basic_polynomial_test.cc  -o basic_polynomial_test $(Numberlib)
# atomic_polynomial_test.o:test/atomic_polynomial_test.cc basic_polynomial.hh basic.hh monomial.hh basic_monomial.hh   variable.hh number.hh
# 	$(CXX) $(CFLAGS) $(IPATHS)  test/atomic_polynomial_test.cc  -c
variable_test:test/variable_test.cc variable.hh
	$(CXX) $(CFLAGS) $(IPATHS) test/variable_test.cc  -o variable_test 
# variable.o: variable.cc variable.hh
# 	$(CXX) $(CFLAGS)  $(IPATHS) variable.cc  -c 
# variable_test.o:test/variable_test.cc variable.hh
# 	$(CXX) $(CFLAGS)  $(IPATHS) test/variable_test.cc  -c
.PHONY : clean
clean :
	-rm $(objects)