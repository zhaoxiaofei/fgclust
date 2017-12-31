#include "kseq.h"

#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

#include <unistd.h>

KSEQ_INIT(int, read)

int main(int argc, char **argv) {
    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;

    std::vector<std::tuple<std::string, std::string>> fastarecords;
    
    FILE *fastafile = fopen(argv[1], "r");
    kseq_t *kseq = kseq_init(fileno(fastafile));
    while ( kseq_read(kseq) >= 0 ) { 
        fastarecords.push_back(std::make_tuple(kseq->name.s, kseq->comment.s));
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
        std::cout << std::get<0>(fastarecords.at(inner)) << "\t" << std::get<0>(fastarecords.at(outer)) << "\t" << sim << "\t" << std::get<1>(fastarecords.at(outer)) << std::endl;
    }
}

