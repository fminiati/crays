OBJ = oxcr.o CRTransport.o

CC     = /opt/local/bin/gcc-mp-10
CXX    = /opt/local/bin/g++-mp-10 
OPTS   = -std=c++17 -O2 #-g -Wall
INCL   = -I/opt/local/include
LIB    = -L/opt/local/lib
LIBS   = $(LIB) -lm -lstdc++
EXEDIR = ./

.cpp.o:
	$(CXX) -c $(INCL) $(OPTS) $<

.c.o:
	$(CC) $(OPTS) -c $<

ox: $(OBJ)
	$(CXX) ${OPTS} -o $(EXEDIR)oxcr.ex $(OBJ) $(LIBS)

clean: 
	\rm *.o *.ex

#------------------------------------------------------------------------------
# list of dependencies
#------------------------------------------------------------------------------

oxcr.o: ../exe/oxcr.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../exe/oxcr.cpp

CRTransport.o: CRTransport.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c CRTransport.cpp
