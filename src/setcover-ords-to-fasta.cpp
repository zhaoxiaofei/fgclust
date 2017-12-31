#include "kseq.h"

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
    int is_out_fa_hdr_in_tsv = false;
    for (int i = 1; i + 1 < argc; i+=2) {
        if (!strcmp("--is-output-fasta-header-in-tsv", argv[i])) {
            is_out_fa_hdr_in_tsv = atoi(argv[i+1]);
        } else {
            abort();
        }
    }

    std::vector<std::tuple<std::string, std::string, std::string>> fastarecords;
    
    FILE *fastafile = fopen(argv[1], "r");
    kseq_t *kseq = kseq_init(fileno(fastafile));
    while ( kseq_read(kseq) >= 0 ) { 
        fastarecords.push_back(std::make_tuple(kseq->name.s, kseq->seq.s, kseq->comment.s));
    }
    kseq_destroy(kseq);
    fclose(fastafile);
    
    std::vector<std::vector<int>> inner_to_outers(fastarecords.size(), std::vector<int>());
    std::string line;
    int dec = 1;
    for (size_t i = 0; i < fastarecords.size(); i++) {
        std::getline(std::cin, line);
        std::stringstream ss(line);
        int inner, outer;
        ss >> inner;
        ss >> outer;
        if (0 == outer) {
            dec = 0;
        }
        inner -= dec;
        outer -= dec;
        inner_to_outers[inner].push_back(outer);
    }

    for (unsigned inner = 0; inner < fastarecords.size(); inner++) {
        if (0 < inner_to_outers[inner].size()) {
            std::cout << ">" << std::get<0>(fastarecords[inner]);
            if (is_out_fa_hdr_in_tsv) {
                for (auto outer : inner_to_outers[inner] ) {
                    std::cout << "\t" << std::get<0>(fastarecords[outer]);
                }
            } else {
                std::cout << "\t" << std::get<2>(fastarecords[inner]);
            }
            std::cout << std::endl << std::get<1>(fastarecords[inner]) << std::endl;
        }
    }
}

