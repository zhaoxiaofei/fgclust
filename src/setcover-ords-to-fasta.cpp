#include "kseq.h"

#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include <unistd.h>

KSEQ_INIT(int, read)

int main(int argc, char **argv) {
    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;

    std::vector<std::pair<std::string, std::string>> fastarecords;
    
    FILE *fastafile = fopen(argv[1], "r");
    kseq_t *kseq = kseq_init(fileno(fastafile));
    while ( kseq_read(kseq) >= 0 ) { 
        fastarecords.push_back(std::make_pair(kseq->name.s, kseq->seq.s));
    }
    kseq_destroy(kseq);
    fclose(fastafile);
    
    std::vector<std::vector<int>> inner_to_outers(fastarecords.size(), std::vector<int>());
    std::string line;
    int dec = 1;
    for (size_t i = 0; i < fastarecords.size(); i++) {
        std::getline(std::cin, line);
        std::stringstream ss(line);
        int inner, outer, sim;
        ss >> inner;
        ss >> outer;
        if (0 == outer) {
            dec = 0;
        }
        inner -= dec;
        outer -= dec;
        inner_to_outers[inner].push_back(outer);
    }

    for (int i = 0; i < fastarecords.size(); i++) {
        if (0 < inner_to_outers[i].size()) {
            std::cout << ">" << fastarecords[i].first;
            for (auto outer : inner_to_outers[i] ) {
                std::cout << "," << fastarecords[outer].first;
            }
            std::cout << std::endl << fastarecords[i].second << std::endl;
        }
    }
}

