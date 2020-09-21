#
# Copyright (C) 2020 Francesco Miniati <francesco.miniati@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
OBJ = oxcr.o CRTransport.o

CC     = /opt/local/bin/gcc-mp-10
CXX    = /opt/local/bin/g++-mp-10 
OPTS   = -std=c++17 -O2 #-g -Wall
INCL   = -I/opt/local/include -I./../src
LIB    = -L/opt/local/lib
LIBS   = $(LIB) -lm -lstdc++
EXEDIR = ../exe/

.cpp.o:
	$(CXX) -c $(INCL) $(OPTS) $<

.c.o:
	$(CC) $(OPTS) -c $<

crays: $(OBJ)
	$(CXX) ${OPTS} -o $(EXEDIR)oxcr.ex $(OBJ) $(LIBS)

clean: 
	\rm *.o ../exe/*.ex

#------------------------------------------------------------------------------
# list of dependencies
#------------------------------------------------------------------------------

oxcr.o: ../exe/oxcr.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../exe/oxcr.cpp

CRTransport.o: ../src/CRTransport.cpp
	$(CXX) $(EMI_INC) $(INCL) $(OPTS) -c ../src/CRTransport.cpp
