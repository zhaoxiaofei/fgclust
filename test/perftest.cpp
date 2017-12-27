#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "edlib.h"
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

#include <seqan/seeds.h>
#include <seqan/align.h>
#include <seqan/sequence.h>

#include <cassert>
#include <iostream>     // std::cout
#include <algorithm>
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include <map>
#include <unordered_map>

#define MIN(a, b) ( (a) < (b) ? (a) : (b) )


typedef uint32_t sign_t;
const uint64_t SIGN_BASE  = 10007L; // 1009L
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L;
uint64_t SIGN_POWER = 0; // derived constan

int SIGN_SIZE = 8;

using namespace seqan;

const int PERC_SIM = 80;

void sign_INIT() {
    SIGN_POWER = 1;
    for (int i = 0; i < SIGN_SIZE; i++) {
        SIGN_POWER = (SIGN_POWER * SIGN_BASE) % SIGN_MOD;
    }
}

const uint64_t sign_init(const char *beg) {
    uint64_t ret = 0;

    for (int i = 0; i < SIGN_SIZE; i++) {
        ret *= SIGN_BASE;
        ret += (uint64_t) beg[i];
        ret %= SIGN_MOD;
    }
    return ret;
}

const uint64_t sign_update(uint64_t hash, char prv, char nxt) {
    hash = hash * SIGN_BASE + (uint64_t) nxt;
    hash += SIGN_MOD;
    hash -= (SIGN_POWER * (uint64_t)prv) % SIGN_MOD;
    return hash % SIGN_MOD;
}

static inline const int calc_perc_seq_sim_seed(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    
    std::map<sign_t, sign_t> sign_to_pos1_map;
    uint64_t sign = sign_init(s1);
    sign_to_pos1_map.insert(std::make_pair((sign_t)sign, (sign_t)0));
    for (int i = SIGN_SIZE; i < strlen(s1); i++) {
        sign = sign_update(sign, s1[i-SIGN_SIZE], s1[i]);
        if ((i - SIGN_SIZE + 1) % 10 == 0) sign_to_pos1_map.insert(std::make_pair((sign_t)sign, (sign_t)(i-SIGN_SIZE+1)));
    }
    
    //std::cerr << "Preseed" << std::endl;

    typedef Seed<Simple>    TSeed;
    typedef SeedSet<TSeed> TSeedSet;
    typedef Iterator<TSeedSet>::Type TIter;
    TSeedSet seedSet;
    
    sign = sign_init(s2);
    if (sign_to_pos1_map.find((sign_t)sign) != sign_to_pos1_map.end()) {
        addSeed(seedSet, TSeed(sign_to_pos1_map.find((sign_t)sign)->second, 0, SIGN_SIZE), Single());
    }
    for (int i = SIGN_SIZE; i < strlen(s2); i+=1) {
        sign = sign_update(sign, s2[i-SIGN_SIZE], s2[i]);
        if (sign_to_pos1_map.find((sign_t)sign) != sign_to_pos1_map.end()) {
            if ( ! addSeed(seedSet, TSeed( sign_to_pos1_map.find((sign_t)sign)->second, i-SIGN_SIZE+1, SIGN_SIZE), 1, Merge()) ) {
                addSeed(seedSet, TSeed( sign_to_pos1_map.find((sign_t)sign)->second, i-SIGN_SIZE+1, SIGN_SIZE), Single());
            }
        }
    }
    
    //std::cerr << "PreChain" << std::endl;
    #if 0
    for (TIter it = begin(seedSet, Standard()); it != end(seedSet, Standard()); ++it) {
        //appendValue(seedChain, *it);
        #if 1
        std::cerr << "(" << beginPositionH(*it) << ", " << endPositionH(*it)
                  << ", " << beginPositionV(*it) << ", " << endPositionV(*it)
                  << ", " << lowerDiagonal(*it) << ", " << upperDiagonal(*it)
                  << ")\n";
        #endif
    }
    #endif
    
    String<TSeed> seedChain;
    chainSeedsGlobally(seedChain, seedSet, SparseChaining());
    // std::cerr << "PostChain" << std::endl;

    //std::cerr << "Postseed" << std::endl;

    Peptide seq1 = s1;
    Peptide seq2 = s2;
    Align<Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
    Blosum62 scoringScheme(-1, -12);
    int score = bandedChainAlignment(align, seedChain, scoringScheme, 5);
    AlignmentStats stats;
    computeAlignmentStats(stats, align, scoringScheme);
    int 
    ret = 100 * stats.numMatches / strlen(s1);
    // ret = 100 * score / 2        / strlen(s1); 
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}



static inline const int calc_perc_seq_sim_seqan(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    Peptide seq1 = s1;
    Peptide seq2 = s2;
    Align<Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
    Blosum62 scoringScheme(-1, -12);
    int score = globalAlignment(align, scoringScheme, AlignConfig<false, true, true, false>(), -10, 10);
    AlignmentStats stats;
    computeAlignmentStats(stats, align, scoringScheme);
    int 
    ret = 100 * stats.numMatches / strlen(s1);
    // ret = 100 * score / 2        / strlen(s1); 
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

static const int calc_perc_seq_sim_editdist(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    int maxEditDist = strlen(s1) * (100 - PERC_SIM) / 100 * 2;
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
    int ret = 100 * (alnlen - editdist) / ( strlen(s1) /*strlen(s1) + 18 */);
    assert (ret >= 0);
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

static const int calc_perc_seq_sim_seqident(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    const int s1len = (int)strlen(s1);
    const int s2len = (int)strlen(s2);
    
    parasail_result_t *result;
    #if 0
    result = parasail_sw_striped_16(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_sw_striped_32(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    } 
    #endif
    #if 1
    if (MIN(s1len, s2len) < 4000) {
        result = parasail_sw_striped_16(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    }
    else {
        result = parasail_sw_striped_32(s1, s1len, s2, s2len, 11, 1, &parasail_blosum62);
    }
    #endif
    assert(!result->saturated);
    #if 0
    int matches = 0;
    for (int i = 0; i < s2len; i++) {
        if (result->matches_row[i] > matches) {
            matches = result->matches_row[i];
        }
    }
    #endif
    int alnscore = result->score;
    int matches = result->matches;
    // int similar = result->similar;
    parasail_result_free(result);

    int 
    ret = 100     * matches  / s1len;
    ret = 100     * alnscore / s1len;
    if (ret <= PERC_SIM) {return 0;}
    return ret;
}

int main(int argc, char **argv) {
    sign_INIT(); 
    #if 0
    std::map<int, int> maps [1000*1000];
    for (int i = 0; i < 1000*1000; i++) {
        assert( maps[i].size() == 0);
    }
    #endif

    int niters = 1000 * 1000 * 10;
    int i; 
    time_t beg, end;

    #if 0
    std::vector<int> vec;
    
    for (i = 0; i < niters; i++) {
        vec.push_back(i);
    }
    std::random_shuffle(vec.begin(), vec.end());
    
    std::map<int, int> omap;
    std::unordered_map<int, int> umap;

    time(&beg);
    for (i = 0; i < niters; i++) {
        omap.insert(std::make_pair(*(vec.begin()+i), i));
    }
    time(&end);
    std::cout << niters << " omap insert " << difftime(end, beg) << std::endl;
    
    time(&beg);
    for (i = 0; i < niters; i++) {
        auto it = omap.find(i);
        assert(it->second >= 0);
    }
    time(&end); 
    std::cout << niters << " omap find " << difftime(end, beg) << std::endl;

    time(&beg);
    for (i = 0; i < niters; i++) {
        auto it = omap.erase(i);
        assert(it > 0);
    }
    time(&end); 
    std::cout << niters << " omap erase " << difftime(end, beg) << std::endl;

    
    time(&beg);
    for (i = 0; i < niters; i++) {
        umap.insert(std::make_pair(*(vec.begin()+i), i));
    }
    time(&end);
    std::cout << niters << " umap insert " << difftime(end, beg) << std::endl;
    
    time(&beg);
    for (i = 0; i < niters; i++) {
        auto it = umap.find(i);
        assert(it->second >=0);
    }
    time(&end); 
    std::cout << niters << " umap find " << difftime(end, beg) << std::endl;

    time(&beg);
    for (i = 0; i < niters; i++) {
        auto it = umap.erase(i);
        assert(it > 0);
    }
    time(&end); 
    std::cout << niters << " umap erase " << difftime(end, beg) << std::endl;
    #endif 
    
    srand(-1);

    #if 1
    char aatable[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    char *seqs[100];
    seqs[0] = (char*)malloc(79999+1);
    for (int j = 0; j < 79999; j++) {
        seqs[0][j] = aatable[rand() % 20];
    }
    seqs[0][79999] = '\0';
    for (i = 1; i < 100; i++) {
        seqs[i] = (char*)malloc(79999+1);
        strcpy(seqs[i], seqs[0]);
        for (int j = 0; j < 79999; j++) {
            if (rand() % 10 == 0) {
                seqs[i][j] = aatable[rand() % 20];
            }
        }
    }

    niters = 1000 / 100; // 5000 / 2000;
    
    {
        srand(-1);
        time(&beg);
        int res = 0;
        for (i = 0; i < niters; i++) {
            int psim = calc_perc_seq_sim_seed(seqs[rand()%100], seqs[0]); 
            //std::cerr << seqs[j] << std::endl;
            //std::cerr << seqs[k] << std::endl;
            res += psim;
        }
        time(&end);
        std::cout << niters << " seed " << difftime(end, beg) << " avgdist=" << (double)res / (double)niters << std::endl;
    }

    {
        srand(-1);
        time(&beg);
        int res = 0;
        for (i = 0; i < niters; i++) {
            int psim = calc_perc_seq_sim_seqan(seqs[rand()%100], seqs[0]); 
            //std::cerr << seqs[j] << std::endl;
            //std::cerr << seqs[k] << std::endl;
            res += psim;
        }
        time(&end);
        std::cout << niters << " seqan " << difftime(end, beg) << " avgdist=" << (double)res / (double)niters << std::endl;
    }


    {
        srand(-1);
        time(&beg);
        int res = 0;
        for (i = 0; i < niters; i++) {
            int dist = calc_perc_seq_sim_editdist(seqs[rand()%100], seqs[0]);
            //std::cerr << seqs[j] << std::endl;
            //std::cerr << seqs[k] << std::endl;
            res += dist;
        }
        time(&end);
        std::cout << niters << " editdist " << difftime(end, beg) << " avgdist=" << (double)res / (double)niters << std::endl;
    } 
    {
        srand(-1);
        time(&beg);
        int res = 0;
        for (i = 0; i < niters; i++) {
            int psim = calc_perc_seq_sim_seqident(seqs[rand()%100], seqs[0]); 
            //std::cerr << seqs[j] << std::endl;
            //std::cerr << seqs[k] << std::endl;
            res += psim;
        }
        time(&end);
        std::cout << niters << " seqident " << difftime(end, beg) << " avgdist=" << (double)res / (double)niters << std::endl;
    }
    #endif
}

