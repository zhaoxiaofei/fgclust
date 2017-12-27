CXX=g++
CXXFLAGS=-O3 -static-libstdc++ -I lib -DGITCOMMIT=\"CXX-$(shell git rev-parse HEAD)-diff$(shell git diff --name-only | wc -l)\" -DCXXVERSION="\"$(shell $(CXX) --version | head -n1)\"" -DDESC="\"$(desc)\""
GPRFLAGS=-O2 -static-libstdc++ -I lib -DGITCOMMIT=\"GPR-$(shell git rev-parse HEAD)-diff$(shell git diff --name-only | wc -l)\" -DCXXVERSION="\"$(shell $(CXX) --version | head -n1)\"" -DDESC="\"$(desc)\"" -pg

all: bin/len-revname-sort.out bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out \
                              bin/fastaseqs-to-distmatrix.gpr

bin/fastaseqs-to-distmatrix.gpr : src/fastaseqs-to-distmatrix.cpp
	$(CXX) src/fastaseqs-to-distmatrix.cpp -o bin/fastaseqs-to-distmatrix.gpr $(GPRFLAGS) lib/edlib.cpp -fopenmp

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
	rm bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out bin/len-revname-sort.out \
	   bin/fastaseqs-to-distmatrix.gpr
	
