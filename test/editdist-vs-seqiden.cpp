#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "kseq.h"
#include "edlib.h"
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

#include <cassert>
#include <iostream>     // std::cout
#include <algorithm>
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include <map>
#include <unordered_map>

KSEQ_INIT(int, read)

const int PERC_SIM = 0;

static const int calc_perc_seq_sim_editdist(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    int maxEditDist = strlen(s1);// * (100 - PERC_SIM) / 100 * 2;
    EdlibAlignResult result;
    result = edlibAlign(s1, strlen(s1), s2, strlen(s2),
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_PATH));
    int editdist = result.editDistance;
    int alnlen = result.alignmentLength;
    //std::cerr << editdist << " from " <<  alnlen << std::endl;
    //char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    //std::cerr << cigar << std::endl;
    edlibFreeAlignResult(result);    

    if (-1 == editdist) { return 0; }
    int ret = 100 * ( (alnlen * 1 + strlen(s1) * 0) - editdist) / ( strlen(s1) /*strlen(s1) + 18 */);
    assert (ret >= 0);
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

static const int calc_perc_seq_sim_seqident(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    const int s1len = (int)strlen(s1);
    const int s2len = (int)strlen(s2);
    
    parasail_result_t *result = parasail_sg_stats_striped_16(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    #if 0
    int matches = 0;
    for (int i = 0; i < s2len; i++) {
        if (result->matches_row[i] > matches) {
            matches = result->matches_row[i];
        }
    }
    #endif
    // int alnscore = result->score;
    int matches = result->matches;
    // int similar = result->similar;
    parasail_result_free(result);

    int ret = 100     * matches  / s1len;
        // ret = 100 / 2 * alnscore / s1len;
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

int main(int argc, char **argv) {
    
    kseq_t *kseq = kseq_init(fileno(stdin));
    int niter = 0;
    char prevseq[50000];
    std::vector<std::pair<int, int>> diffs;
    while (kseq_read(kseq) >= 0) {
        if (niter % 1000 == 0) {
            strcpy(prevseq, kseq->seq.s);
        } else if (niter % 1000 == 1) {
            int measure1 = calc_perc_seq_sim_editdist(kseq->seq.s, prevseq);
            int measure2 = calc_perc_seq_sim_seqident(kseq->seq.s, prevseq);
            if (measure1 > 45 && measure2 > 45) diffs.push_back(std::make_pair(measure1, measure2));
        }
        niter++;
    }
    kseq_destroy(kseq);
    for (auto diff : diffs) {
        std::cout << diff.first << "\t" <<diff.second << std::endl;
    }
}

