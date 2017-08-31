CXX=g++
CXXFLAGS=-O3 -I lib -DGITCOMMIT=\"$(shell git rev-parse HEAD)\" -DCXXVERSION="\"$(shell $(CXX) --version | head -n1)\""

all: bin/len-revname-sort.out bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out

bin/fastaseqs-to-distmatrix.out : src/fastaseqs-to-distmatrix.cpp
	$(CXX) src/fastaseqs-to-distmatrix.cpp -o bin/fastaseqs-to-distmatrix.out $(CXXFLAGS) lib/edlib.cpp -fopenmp

bin/linsetcover.out             : src/linsetcover.cpp
	$(CXX) src/linsetcover.cpp             -o bin/linsetcover.out             $(CXXFLAGS)

bin/setcover-ords-to-hdrs.out   : src/setcover-ords-to-hdrs.cpp
	$(CXX) src/setcover-ords-to-hdrs.cpp   -o bin/setcover-ords-to-hdrs.out   $(CXXFLAGS)

bin/setcover-ords-to-fasta.out  : src/setcover-ords-to-fasta.cpp
	$(CXX) src/setcover-ords-to-fasta.cpp  -o bin/setcover-ords-to-fasta.out  $(CXXFLAGS)

bin/len-revname-sort.out        : src/len-revname-sort.cpp
	$(CXX) src/len-revname-sort.cpp        -o bin/len-revname-sort.out        $(CXXFLAGS)

.PHONY clean:
	rm bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out bin/len-revname-sort.out

