
all: bin/len-revname-sort.out bin/fastaseqs-to-distmatrix.out bin/linsetcover.out bin/setcover-ords-to-hdrs.out bin/setcover-ords-to-fasta.out

bin/fastaseqs-to-distmatrix.out : src/fastaseqs-to-distmatrix.cpp
	g++ src/fastaseqs-to-distmatrix.cpp -o bin/fastaseqs-to-distmatrix.out -O3 -I lib lib/edlib.cpp -fopenmp

bin/linsetcover.out : src/linsetcover.cpp
	g++ src/linsetcover.cpp -o bin/linsetcover.out -O3 

bin/setcover-ords-to-hdrs.out : src/setcover-ords-to-hdrs.cpp
	g++ src/setcover-ords-to-hdrs.cpp -o bin/setcover-ords-to-hdrs.out -O3 -I lib

bin/setcover-ords-to-fasta.out : src/setcover-ords-to-fasta.cpp
	g++ src/setcover-ords-to-fasta.cpp -o bin/setcover-ords-to-fasta.out -O3 -I lib

bin/len-revname-sort.out : src/len-revname-sort.cpp
	g++ src/len-revname-sort.cpp -o bin/len-revname-sort.out -O3 -I lib

