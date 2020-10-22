CXX=g++ 
IPATHS=-I../ 
CFLAGS=-O3 -flto -DNDEBUG  
Numberlib= -lgmpxx -lgmp  
objects= variable_test basic_polynomial_test   monomial_test polynomial_test
clpoly_hh=../clpoly/associatedgraph.hh ../clpoly/basic.hh ../clpoly/basic_monomial.hh ../clpoly/basic_polynomial.hh ../clpoly/monomial.hh ../clpoly/monomial_order.hh ../clpoly/number.hh ../clpoly/polynomial_gcd.hh ../clpoly/polynomial.hh ../clpoly/polynomial_.hh ../clpoly/random.hh ../clpoly/upolynomial.hh ../clpoly/variable.hh
polynomial_GCD_test:polynomial_GCD_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) polynomial_GCD_test.cc -o polynomial_GCD_test $(Numberlib)
polynomial_res_test:polynomial_res_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) polynomial_res_test.cc -o polynomial_res_test $(Numberlib)
polynomial_prem_test:polynomial_prem_test.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) polynomial_prem_test.cc -o polynomial_prem_test $(Numberlib)
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