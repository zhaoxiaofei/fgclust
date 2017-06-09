#include "kseq.h"

#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include <unistd.h>

KSEQ_INIT(int, read)

int main(int argc, char **argv) {
    
    std::vector<std::pair<std::string, std::string>> fastarecords;
    
    FILE *fastafile = fopen(argv[1], "r");
    kseq_t *kseq = kseq_init(fileno(fastafile));
    while ( kseq_read(kseq) >= 0 ) { 
        fastarecords.push_back(std::make_pair(kseq->name.s, kseq->seq.s));
    }
    kseq_destroy(kseq);
    fclose(fastafile);
    
    std::string line;
    int dec = 1;
    for (size_t i = 0; i < fastarecords.size(); i++) {
        std::getline(std::cin, line);
        std::stringstream ss(line);
        int inner, outer, sim;
        ss >> inner; //assert(inner > 0);
        ss >> outer; //assert(outer > 0);
        ss >> sim;
        if (0 == outer) {
            dec = 0;
        }
        inner -= dec;
        outer -= dec;
        assert(inner >= 0);
        assert(outer >= 0);
        assert(sim > 0);
        std::cout << inner << "\t" << outer << "\t" << sim << std::endl;
    }
}

