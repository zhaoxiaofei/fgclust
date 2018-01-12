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
#include <random>
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include <map>
#include <tuple>
#include <unordered_map>

KSEQ_INIT(int, read)

const int PERC_SIM = 0;

static const int calc_perc_seq_sim_editdist(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    int maxEditDist = strlen(s1);// * (100 - PERC_SIM) / 100 * 2;
    EdlibAlignResult result;
    result = edlibAlign(s1, strlen(s1), s2, strlen(s2),
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    int editdist = result.editDistance;
    int alnlen = result.alignmentLength;
    //std::cerr << editdist << " from " <<  alnlen << std::endl;
    //char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    //std::cerr << cigar << std::endl;
    edlibFreeAlignResult(result);    

    if (-1 == editdist) { return 0; }
    int ret = 1000 * ( (alnlen * 1 + strlen(s1) * 0) - editdist) / ( strlen(s1) /*strlen(s1) + 18 */);
    assert (ret >= 0);
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

static const int calc_perc_seq_sim_seqident(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    const int s1len = (int)strlen(s1);
    const int s2len = (int)strlen(s2);
    
    parasail_result_t *result = parasail_nw_stats_striped_16(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    int alnscore = result->score;
    int matches = result->matches;
    // int similar = result->similar;
    int length = result->length;
    parasail_result_free(result);
    int ret;
    // ret = 100     * matches  / length;
    // ret = 100 / 2 * alnscore / s1len;
    ret = 1000     * matches / s1len;
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

int main(int argc, char **argv) { 
    std::vector<std::tuple<std::string, std::string>> name_seq_vec;
    kseq_t *kseq = kseq_init(fileno(stdin));
    std::vector<std::tuple<int, int, std::string, std::string>> diffs;
    while (kseq_read(kseq) >= 0) {
        name_seq_vec.push_back(std::make_tuple(std::string(kseq->name.s), std::string(kseq->seq.s)));
    }
    kseq_destroy(kseq);
    int step = name_seq_vec.size() / (10*1000);
    for (int i = 0; i < name_seq_vec.size(); i+=step) {
        auto name_seq = name_seq_vec[i];
        for (int j = 1; i - j >= 0; j *= 2) {
            auto name_seq2 = name_seq_vec[i-j];
            int measure1 = calc_perc_seq_sim_editdist(std::get<1>(name_seq).c_str(), std::get<1>(name_seq2).c_str());
            int measure2 = calc_perc_seq_sim_seqident(std::get<1>(name_seq).c_str(), std::get<1>(name_seq2).c_str());
            diffs.push_back(std::make_tuple(measure1, measure2, std::get<0>(name_seq), std::get<0>(name_seq2)));
        }
        for (int j = 1; i + j < name_seq_vec.size(); j *= 2) {
            auto name_seq2 = name_seq_vec[i+j]; 
            int measure1 = calc_perc_seq_sim_editdist(std::get<1>(name_seq).c_str(), std::get<1>(name_seq2).c_str());
            int measure2 = calc_perc_seq_sim_seqident(std::get<1>(name_seq).c_str(), std::get<1>(name_seq2).c_str());
            diffs.push_back(std::make_tuple(measure1, measure2, std::get<0>(name_seq), std::get<0>(name_seq2)));
        }
    }
    std::shuffle(diffs.begin(), diffs.end(), std::default_random_engine(-1));
    int cardinalities[10000];
    memset(cardinalities, 0, 10000 * sizeof(*cardinalities));
    for (auto diff : diffs) {
        cardinalities[std::get<0>(diff)]++;
        std::cout << std::get<0>(diff) << "\t" << std::get<1>(diff) << "\t" << std::get<2>(diff) << "\t" << std::get<3>(diff) << "\t" 
                  << cardinalities[std::get<0>(diff)] << std::endl;
    }
}

