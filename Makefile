CXX=g++
CXXFLAGS=-static-libstdc++ -I lib \
	-DGITCOMMIT=\"CXX-$(shell git rev-parse HEAD)-diff$(shell git diff --name-only | wc -l)\" \
	-DCXXVERSION="\"$(shell $(CXX) --version | head -n1)\"" \
	-DDESC="\"$(desc)\""

OUTFLAGS=$(CXXFLAGS) -O3 -Wall
GPRFLAGS=$(CXXFLAGS) -O2 -pg

LIBFLAGS=-I lib/edlib-1.2.1/edlib/include lib/edlib-1.2.1/edlib/src/edlib.cpp -fopenmp

all: bin/len-revname-sort.out bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out \
                              bin/fastaseqs-to-distmatrix.gpr

bin/fastaseqs-to-distmatrix.gpr : src/fastaseqs-to-distmatrix.cpp
	$(CXX) src/fastaseqs-to-distmatrix.cpp -o bin/fastaseqs-to-distmatrix.gpr $(GPRFLAGS) $(LIBFLAGS)

bin/fastaseqs-to-distmatrix.out : src/fastaseqs-to-distmatrix.cpp
	$(CXX) src/fastaseqs-to-distmatrix.cpp -o bin/fastaseqs-to-distmatrix.out $(OUTFLAGS) $(LIBFLAGS)

bin/linsetcover.out             : src/linsetcover.cpp
	$(CXX) src/linsetcover.cpp             -o bin/linsetcover.out             $(OUTFLAGS)

bin/setcover-ords-to-hdrs.out   : src/setcover-ords-to-hdrs.cpp
	$(CXX) src/setcover-ords-to-hdrs.cpp   -o bin/setcover-ords-to-hdrs.out   $(OUTFLAGS)

bin/setcover-ords-to-fasta.out  : src/setcover-ords-to-fasta.cpp
	$(CXX) src/setcover-ords-to-fasta.cpp  -o bin/setcover-ords-to-fasta.out  $(OUTFLAGS)

bin/len-revname-sort.out        : src/len-revname-sort.cpp
	$(CXX) src/len-revname-sort.cpp        -o bin/len-revname-sort.out        $(OUTFLAGS)

.PHONY clean:
	rm bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out bin/len-revname-sort.out \
	   bin/fastaseqs-to-distmatrix.gpr
	
