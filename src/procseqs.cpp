#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include <unistd.h>

#include "kseq.h"

KSEQ_INIT(int, read)

class Sequence {
public:
    std::string name;
    std::string seq;
    std::string comment;
    Sequence(const kseq_t *kseq);
};

Sequence::Sequence(const kseq_t *kseq) {
    this->name = std::string(kseq->name.s);
    this->seq = std::string(kseq->seq.s);
    //std::reverse(this->seq.begin(), this->seq.end());
    if (kseq->comment.l) {
        this->comment = std::string(kseq->comment.s);
    } else {
        this->comment = "";
    }
}

struct 
{
    std::hash<std::string> ptr_hash;
    bool operator()(Sequence *a, Sequence *b)
    {
        int alen = a->seq.length();
        int blen = b->seq.length();
        if (alen != blen) { return alen > blen; }
        auto ahash = ptr_hash(a->name);
        auto bhash = ptr_hash(b->name);
        if (ahash != bhash) { return ahash < bhash; }
        else { return a->name < b->name; }
    }
}
customLess;

int main(int argc, char **argv) {
    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;

    int PROCSEQS_ORDER = 1;
    for (int i = 1; i < argc; i += 2) {
        if (i+1 < argc && !strcmp("--procseqs-order", argv[i])) {
            PROCSEQS_ORDER = atoi(argv[i+1]);
        } else {
            std::cerr << "Program : " << argv[0] << std::endl;
            std::cerr << "  version " << GITCOMMIT << " compiled by " << CXXVERSION << std::endl;
            std::cerr << "Command-line arguments with [default-values]:" << std::endl;
            std::cerr << "  --procseqs-order\t:1 and 2 mean by pseudorandom order and by decreasing sequence length, respectively. [" << PROCSEQS_ORDER << "]" << std::endl;
            exit(-1);
        }
    }
    kseq_t *kseq = kseq_init(fileno(stdin));
    std::vector<Sequence*> seqs;
    while (kseq_read(kseq) >= 0)
    {
        Sequence *sequence = new Sequence(kseq);
        seqs.push_back(sequence);
        if (!(seqs.size() & (seqs.size() - 1))) { std::cerr << "sort : processed " << seqs.size() << " sequences." << std::endl; }
    }
    kseq_destroy(kseq);
    if (seqs.size() == 0) {
        std::cerr << "ERROR: the input fasta file has no sequence!" << std::endl;
        return -1;
    }
    if (1 == PROCSEQS_ORDER) {
        std::mt19937 g(7);
        std::shuffle(seqs.begin(), seqs.end(), g);
    } else if (2 == PROCSEQS_ORDER) {
        std::sort(seqs.begin(), seqs.end(), customLess);
    }
    for (auto seq : seqs)
    {
        //std::reverse(seq.seq.begin(), seq.seq.end());
        std::cout << ">" << seq->name << " " << seq->comment << std::endl;
        std::cout << seq->seq << std::endl;
    }
}

