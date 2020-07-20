CXX=g++ 
IPATHS=-I../ 
CFLAGS=-O3 -flto -DNDEBUG  
Numberlib= -lgmpxx -lgmp  
objects= variable_test basic_polynomial_test   monomial_test polynomial_test
clpoly_hh=../clpoly/polynomial.hh ../clpoly/monomial.hh ../clpoly/basic_monomial.hh  ../clpoly/basic.hh ../clpoly/variable.hh  ../clpoly/basic_polynomial.hh ../clpoly/number.hh ../clpoly/monomial_order.hh
polynomial_test:polynomial_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) polynomial_test.cc -o polynomial_test $(Numberlib)
# polynomial_test.o:polynomial_test.cc polynomial.hh monomial.hh basic_monomial.hh  basic.hh variable.hh  basic_polynomial.hh number.hh
# 	$(CXX) $(CFLAGS) $(IPATHS) polynomial_test.cc -c
monomial_test:monomial_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS)  $(IPATHS) monomial_test.cc -o monomial_test
# monomial_test.o:monomial_test.cc monomial.hh basic_monomial.hh  basic.hh variable.hh
# 	$(CXX) $(CFLAGS) $(IPATHS) monomial_test.cc -c
basic_polynomial_test:basic_polynomial_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS)  $(IPATHS)  basic_polynomial_test.cc  -o basic_polynomial_test $(Numberlib)
# atomic_polynomial_test.o:atomic_polynomial_test.cc basic_polynomial.hh basic.hh monomial.hh basic_monomial.hh   variable.hh number.hh
# 	$(CXX) $(CFLAGS) $(IPATHS)  atomic_polynomial_test.cc  -c
variable_test:variable_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) variable_test.cc  -o variable_test 
# variable.o: variable.cc variable.hh
# 	$(CXX) $(CFLAGS)  $(IPATHS) variable.cc  -c 
# variable_test.o:variable_test.cc variable.hh
# 	$(CXX) $(CFLAGS)  $(IPATHS) variable_test.cc  -c
.PHONY : clean
clean :
	-rm $(objects)