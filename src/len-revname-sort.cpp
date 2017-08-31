#include <algorithm>
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
    bool operator()(Sequence a, Sequence b)
    {
        int alen = a.seq.length();
        int blen = b.seq.length();
        if (alen != blen) { return alen > blen; }
        else { return a.name < b.name; }
    }
}
customLess;

int main(int argc, char **argv) {
    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;

    bool israndom = false;
    for (int i = 1; i < argc; i++) {
        if (!strcmp("--israndom", argv[i])) {
            israndom = true;
        }
    }
    kseq_t *kseq = kseq_init(fileno(stdin));
    std::vector<Sequence> seqs;
    while (kseq_read(kseq) >= 0)
    {
        Sequence sequence(kseq);
        seqs.push_back(sequence);
        if (!(seqs.size() & (seqs.size() - 1))) { std::cerr << "sort : processed " << seqs.size() << " sequences." << std::endl; }
    }
    kseq_destroy(kseq);
    if (israndom) {
        std::mt19937 g(7);
        std::shuffle(seqs.begin(), seqs.end(), g);
    } else {
        std::sort(seqs.begin(), seqs.end(), customLess);
    }
    for (auto seq : seqs)
    {
        //std::reverse(seq.seq.begin(), seq.seq.end());
        std::cout << ">" << seq.name << " " << seq.comment << std::endl;
        std::cout << seq.seq << std::endl;
    }
}

