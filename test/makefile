.PHONY : clean
all:
CXX=g++ 
IPATHS=-I../ 
CFLAGS=-O3 -flto -DNDEBUG  
Numberlib= -lgmpxx -lgmp  
objects= $(wildcard *test)
clpoly_hh=$(wildcard ../clpoly/*.hh)
%:%.cc $(clpoly_hh)
	$(CXX) $(CFLAGS) $(IPATHS) $< -o $@ $(Numberlib)


clean :
	-rm $(objects)