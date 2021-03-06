//#include "edlib.h"
#include "kseq.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <unistd.h>

#include "edlib.h"

#ifndef INDEXWORD
#define INDEXWORD 1
#endif

#ifndef SHORTWORD
#define SHORTWORD 0
#endif
#ifndef ENTROHASH
#define ENTROHASH 0
#endif
#ifndef VARSIGNS
#define VARSIGNS 0
#endif
#ifndef ORDER_AWARE_FILT
#define ORDER_AWARE_FILT 0
#endif

#ifndef NUM_SIGNATURES
#if SHORTWORD
#define NUM_SIGNATURES (100)
#elif ENTROHASH
#define NUM_SIGNATURES (256)
#elif VARSIGN
#define NUM_SIGNATURES (64)
#else
#define NUM_SIGNATURES (32)
#endif
#endif

void *xmalloc(size_t size) {
    void *ret = malloc(size);
    if (NULL == ret) {
        fprintf(stderr, "malloc failed!\n");
        abort();
    }
    return ret;
}

void *xrealloc(void *ptr, size_t size) {
    void *ret = realloc(ptr, size);
    if (NULL == ret) {
        fprintf(stderr, "realloc failed!\n");
        abort();
    }
    return ret;
}

const auto MIN(const auto a, const auto b) { return a < b ? a : b; }
const auto MAX(const auto a, const auto b) { return a > b ? a : b; }

const auto SQUARE(const auto v) { return v * v; }
void SWAP(auto &a, auto &b) {auto tmp = a; a = b; b = tmp; }

KSEQ_INIT(int, read)

// constants

const unsigned int BATCH_COVRAT = 64;
const unsigned int BATCHSIZE_INI = 256; // 15999;
const unsigned int BATCHSIZE_INC = 16; //8; // 2; // 15999;

const uint64_t PRIME_BASE = 48271L;
const uint64_t SIGN_BASE  = 48271L; 
const uint64_t PRIME_MOD  = (0x1L << 31L) - 1L; 
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L; 

// variables that are initialized once from constants and other variables

uint64_t SEED_POWERS[30] = {0};
uint64_t SIGN_POWER = 0;
char ALPHA_TYPE_TO_CHAR_TO_REDUCED[3][256]; // derived from ALPHA_TYPE_TO_SIZE
int SIGN_SD_CNTMINS[NUM_SIGNATURES+1];

int SEED_HIGH_LEN_PERC = 0;

// variables that are initialized from command line args

int ALPHA_TYPE_TO_SIZE[] = {10, 10, 10};

int ATTEMPT_BASE = 0;
double ATTEMPT_INI = 0.05; //50;
double ATTEMPT_INC = 0.05; //50;
double ATTEMPT_MAX = 0.05; //50;

double COV_SRC_ADA = -1; // 0.005;
unsigned int COV_SRC_MAX = 0; // 5+1; // 8; // 5;
unsigned int COV_SNK_MAX = INT_MAX;

unsigned int DBFILT_MINSEED = INT_MAX; // 1000*1000; lower -> more filtering, more time saving later
int DBFILT_SUBSAMP = 10*1000; // 50*50; // 800; // lower -> less filtering accuracy, less time
// int DBFILT_ATTEMPT = 50; // 10; // lower -> less filtering accuracy, less time // not useful because same as ATTEMPT_* 
int DBFILT_TIMEFAC = 10*1000; // (50*50-1)/100; // 5; // lower -> less filtering accuracy

unsigned int IDXENTRY_ITMAX = 0; // similar to http://bio-bwa.sourceforge.net/bwa.shtml parameter -c

int LEN_PERC_SRC = -1;
int LEN_PERC_SNK = -1;

uint64_t SEED_N_PER_SEQ = 0;
double SEED_EVALUE = 10; // 1; // 20;
double SEED_BORDER = 0.02; // not used
int SEED_LENGTH = 10; // can be overriden after determination of db size 
int SEED_MINCNT = 40; // 30; // cannot be overriden after determination of db size

int SEQTYPE = 0; // guessed from input by default

int SIGN_CHCOV_MAX = 8;
int SIGN_LENGTH = -1;
int SIGN_CNTMIN = 1; 
int SIGN_ZVAL = 3;

int SIM_PERC = -1;
int SIM_BASE = -1; // 0; // 25;
int SIM_DIFF = -1;
int SIM_ZVAL = 0; // 10 * 100; // 11; // real z-score threshold is about 75% of this value due to innacuracy of normal approximation of binomial

bool ZVAL_AS_SIM = false;

// derived vars
uint64_t DBENTRY_CNT = 0; // can be overriden after determination of db size

void showparams() {
    std::cerr << "Parameter values are as follows: " << std::endl; 
    
    std::cerr << "\tBATCHSIZE_INI = " << BATCHSIZE_INI << std::endl;
    std::cerr << "\tNUM_SIGNATURES = " << NUM_SIGNATURES << std::endl;

    std::cerr << "\tALPHA_TYPE_TO_SIZE[0] = ALPHASIZE_SEED = " << ALPHA_TYPE_TO_SIZE[0] << std::endl;
    std::cerr << "\tALPHA_TYPE_TO_SIZE[1] = ALPHASIZE_SIGN = " << ALPHA_TYPE_TO_SIZE[1] << std::endl;
    std::cerr << "\tALPHA_TYPE_TO_SIZE[2] = ALPHASIZE_SIM  = " << ALPHA_TYPE_TO_SIZE[2] << std::endl;

    std::cerr << "\tATTEMPT_INI = " << ATTEMPT_INI << std::endl;
    std::cerr << "\tATTEMPT_INC = " << ATTEMPT_INC << std::endl;
    std::cerr << "\tATTEMPT_MAX = " << ATTEMPT_MAX << std::endl;
    
    std::cerr << "\tCOV_SRC_ADA = " << COV_SRC_ADA << std::endl;
    std::cerr << "\tCOV_SRC_MAX = "   << COV_SRC_MAX   << std::endl;
    std::cerr << "\tCOV_SNK_MAX = "   << COV_SNK_MAX   << std::endl;
   
    std::cerr << "\tDBFILT_MINSEED = " << DBFILT_MINSEED << std::endl;
    std::cerr << "\tDBFILT_SUBSAMP = " << DBFILT_SUBSAMP << std::endl; 
    std::cerr << "\tDBFILT_TIMEFAC = " << DBFILT_TIMEFAC << std::endl;
    //std::cerr << "\tDBFILT_ATTEMPT = " << DBFILT_ATTEMPT << std::endl;
    
    std::cerr << "\tIDXENTRY_ITMAX = " << IDXENTRY_ITMAX << std::endl;
    
    std::cerr << "\tLEN_PERC_SRC = " << (int)LEN_PERC_SRC << std::endl;
    std::cerr << "\tLEN_PERC_SNK = " << (int)LEN_PERC_SNK << std::endl;

    std::cerr << "\tSEED_N_PER_SEQ = " << SEED_N_PER_SEQ << std::endl;
    std::cerr << "\tSEED_EVALUE = "    << SEED_EVALUE    << std::endl;
    std::cerr << "\tSEED_BORDER = "    << SEED_BORDER    << std::endl;
    std::cerr << "\tSEED_LENGTH = "    << SEED_LENGTH    << std::endl;
    std::cerr << "\tSEED_MINCNT = "    << SEED_MINCNT    << std::endl;
    
    std::cerr << "\tSEQTYPE = " << SEQTYPE << std::endl;

    std::cerr << "\tSIGN_CHCOV_MAX = " << SIGN_CHCOV_MAX << std::endl;
    std::cerr << "\tSIGN_LENGTH = " << SIGN_LENGTH       << std::endl;
    std::cerr << "\tSIGN_CNTMIN = " << SIGN_CNTMIN       << std::endl;
    std::cerr << "\tSIGN_ZVAL = "   << SIGN_ZVAL         << std::endl;

    std::cerr << "\tSIM_PERC = " << (int)SIM_PERC << std::endl;
    std::cerr << "\tSIM_BASE = " << (int)SIM_BASE << std::endl;
    std::cerr << "\tSIM_DIFF = " << (int)SIM_DIFF << std::endl;
    std::cerr << "\tSIM_ZVAL = " << (int)SIM_ZVAL << std::endl;

    std::cerr << "\tZVAL_AS_SIM = " << ZVAL_AS_SIM << std::endl;    
}

void show_usage(const int argc, const char *const *const argv) {
    std::cerr << "Program : " << argv[0] << std::endl;
    std::cerr << "  version " << GITCOMMIT << " compiled by " << CXXVERSION << std::endl;
    std::cerr << "Command-line arguments with [default-values]:" << std::endl;
    
    std::cerr << "  --alphasize-seed\t: alphabet size for computation of seed ["                << ALPHA_TYPE_TO_SIZE[0] << "]" << std::endl;
    std::cerr << "  --alphasize-sign\t: alphabet size for computation of sign ["                << ALPHA_TYPE_TO_SIZE[1] << "]" << std::endl;
    std::cerr << "  --alphasize-sim \t: alphabet size for computation of sequence similarity ["  << ALPHA_TYPE_TO_SIZE[2] << "]" << std::endl;
    
    std::cerr << "  --attempt-base  \t: initial number of attempts. ["                             << ATTEMPT_BASE << "]" << std::endl;
    std::cerr << "  --attempt-ini   \t: initial percent of attempts. ["                            << ATTEMPT_INI << "]" << std::endl;
    std::cerr << "  --attempt-inc   \t: percent of attempts incremented per true positive hits. [" << ATTEMPT_INC << "]" << std::endl;
    std::cerr << "  --attempt-max   \t: percent of attempts capped at this maximum value. ["       << ATTEMPT_MAX << "]" << std::endl;
    
    std::cerr << "  --cov-src-ada   \t: fraction of change to cov-src-max per observation. <0 means no change. [" << COV_SRC_ADA << "]" << std::endl;
    std::cerr << "  --cov-snk-max   \t: max number of times that the covered sequence can be covered.["           << COV_SNK_MAX << "]" << std::endl;
    std::cerr << "  --cov-src-max   \t: max number of times that the covering sequence can be coverered. ["       << COV_SRC_MAX << "]" << std::endl;
   
    std::cerr << "  --dbfilt-minseed\t: minimum number of times a seed occurs to trigger seed pruning.["           << DBFILT_MINSEED << "]" << std::endl;
    std::cerr << "  --dbfilt-subsamp\t: number of sequence pairs (SP) subsampled for seed pruning .["              << DBFILT_SUBSAMP << "]" << std::endl;
    // std::cerr << "  --dbfilt-attempt\t: number of attempts on most likely (minhash) similar SP for seed pruning.[" << DBFILT_ATTEMPT << "]" << std::endl;
    std::cerr << "  --dbfilt-truepos\t: a seed occurring this times more needs 1\% increase in true positive SP (map 0 to 0).[" << DBFILT_TIMEFAC << "]" << std::endl;

    std::cerr << "  --idxentry-itmax\t: maximum number of times a db-index entry is iterated by a sequence. [" << IDXENTRY_ITMAX << "]" << std::endl;

    std::cerr << "  --seqtype       \t: 1, 2, and 3 mean input is protein, RNA, and DNA, respectively (auto detected). [" << SEQTYPE << "]" << std::endl;

    std::cerr << "  --len-perc-src  \t: A cannot cover B if len(A) * len-perc-src > len(B) * 100. [" << LEN_PERC_SRC << "]" << std::endl; 
    std::cerr << "  --len-perc-snk  \t: A cannot cover B if len(B) * len-perc-snk > len(A) * 100. [" << LEN_PERC_SNK << "]" << std::endl; 
    
    std::cerr << "  --seed-n-per-seq\t: number of seed hashtable entries per sequence. [" << SEED_N_PER_SEQ << "]" << std::endl;
    std::cerr << "  --seed-evalue   \t: evalue for seed hit. 0 or negative value means do not use this parameter [" << SEED_EVALUE << "]" << std::endl;
    std::cerr << "  --seed-border   \t: fraction of sequences that are on the borderline of clustering cutoff (not used). ["   << SEED_BORDER << "]" << std::endl;
    std::cerr << "  --seed-length   \t: length of an indexed seed. Overwritten by nonzero seed-evalue. ["           << SEED_LENGTH << "]" << std::endl;
    std::cerr << "  --seed-mincnt   \t: minimum number of seeds per sequence. Overwritten by nonzero seed-evalue [" << SEED_MINCNT << "]" << std::endl;
    
    std::cerr << "  --sign-chcov-max\t: max ratio of sequence length to number of minhash values ["         << SIGN_CHCOV_MAX      << "]" << std::endl;
    std::cerr << "  --sign-length   \t: length of k-mers for computing minhash values. ["                      << SIGN_LENGTH         << "]" << std::endl;
    std::cerr << "  --sign-cntmin   \t: minimum number of minhash values to trigger sequence search [" << SIGN_CNTMIN << "]" << std::endl; 
    std::cerr << "  --sign-zval     \t: minimum number of minhash values to trigger sequence search [" << SIGN_ZVAL << "]" << std::endl; 

    std::cerr << "  --sim-zval      \t: percent of standard deviations above similarity by chance for a sequence to cover another ["   << SIM_ZVAL << "]" << std::endl;
    std::cerr << "                  \t: the effective value is 0.2 to 0.4 times less than this set value due to normal approximation of skewed binomial." << std::endl;
    std::cerr << "  --sim-base      \t: A covers B only if (len(B) - edit-distance(A, B)) >= sim_perc * len(B) + sim_base [sim_base="  << SIM_BASE << "]" << std::endl;
    std::cerr << "  --sim-perc      \t:   where extra residues at both ends of A are not penalized. [sim_perc="                        << SIM_PERC << "]" << std::endl;
    std::cerr << "  --sim-diff      \t: percent of deletions in the covering sequence that are ignored, 0 means skip alignment ["      << SIM_DIFF << "]" << std::endl;

    std::cerr << "  --zval-as-sim   \t: used the sim-zval as similarity threshold [" << ZVAL_AS_SIM   << "]" << std::endl;

    std::cerr << "Note: illegal values such as 0 and 1 (depending on context) mean dependence to other parameters or to the input." << std::endl;
    exit(-1);
}

void alphareduce(const char *const strarg, const int reducetype) {
    const char *str = strarg;
    for (; *str; str++) {
        ALPHA_TYPE_TO_CHAR_TO_REDUCED[reducetype][(int)*str] = *strarg;
    }
}

void hash_sign_INIT() {
    int i;
    for (int seedlen = SEED_LENGTH; seedlen <= SEED_LENGTH + 1; seedlen += 1) {
        SEED_POWERS[seedlen] = 1;
        for (i = 0; i < seedlen; i++) {
            SEED_POWERS[seedlen] = (SEED_POWERS[seedlen] * PRIME_BASE) % PRIME_MOD;
        }
    }
    SIGN_POWER = 1;
    for (i = 0; i < SIGN_LENGTH; i++) { 
        SIGN_POWER = (SIGN_POWER * SIGN_BASE) % SIGN_MOD;
    }
}

const uint64_t hash_init(const char *beg, const int seedlen) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < seedlen; i++) {
        ret *= PRIME_BASE;
        ret += (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)beg[i]];
        ret %= PRIME_MOD; 
    }
    return ret;
}

const uint64_t hash_update(uint64_t hash, char prv, char nxt, const int seedlen) {
    hash = hash * PRIME_BASE + (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)nxt];
    hash += PRIME_MOD;
    hash -= (SEED_POWERS[seedlen] * (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)prv]) % PRIME_MOD;
    return hash % PRIME_MOD;
}

const uint64_t sign_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < SIGN_LENGTH; i++) {
        ret *= SIGN_BASE;
        ret += (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][(int)beg[i]];
        ret %= SIGN_MOD; 
    }
    return ret;
}

const uint64_t sign_update(uint64_t hash, char prv, char nxt) {
    hash = hash * SIGN_BASE + (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][(int)nxt];
    hash += SIGN_MOD;
    hash -= (SIGN_POWER * (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][(int)prv]) % SIGN_MOD;
    return hash % SIGN_MOD;
}

typedef struct {
    uint32_t bufsize;    // 4 bytes
    uint32_t size;       // 4 bytes
    uint32_t *seqidxs;   // 8 bytes ++
}
seed_t; // 16 bytes++

typedef struct {
    char *name;
    char *seq;
#if ENTROHASH
    uint64_t compressedsign;
    uint64_t compressedsign2;
    uint64_t compressedsign3;
    uint64_t compressedsign4;
#elif VARSIGNS 
    uint16_t *varsigns;
#elif !SHORTWORD
    uint16_t signatures[NUM_SIGNATURES];
#endif
    uint32_t coveredcnt;
    // uint32_t coveringcnt;
    uint32_t seqlen;
}
seq_t; // 16+16*4 bytes +++

const unsigned int seqlen_to_n_varsigns(unsigned int seqlen) {
    return MIN(MAX((unsigned int)NUM_SIGNATURES, seqlen / (unsigned int)SIGN_CHCOV_MAX), seqlen - (unsigned int)SIGN_LENGTH + 1);
}

const bool fail_len_cutoff(const seq_t *src, const seq_t *snk) {
    return src->seqlen * LEN_PERC_SRC > snk->seqlen * 100 || snk->seqlen * LEN_PERC_SNK > src->seqlen * 100;
}

static const int comp_perc_seq_sim_editdist(const seq_t *src, const seq_t *snk, const double matchprob, const double *chprobs) {
    
    int maxEditDist = snk->seqlen - ceil((1 - DBL_EPSILON) * (sqrt(SQUARE((double)SIM_BASE) + SQUARE((double)(snk->seqlen * (SIM_PERC - SIM_DIFF)) / 100.0))));
    if (maxEditDist < 0) { return -1; }
    
    unsigned int i;
    for (i = 0; i < snk->seqlen; i++) {
        if (chprobs[(int)snk->seq[i]]) { break; }
    }
    if (snk->seqlen == i) { return -2; }

    EdlibAlignResult result;
    int editdist;
    if (SIM_DIFF) {
        result = edlibAlign(snk->seq, snk->seqlen, src->seq, src->seqlen,
                            edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        int match = result.alignmentLength - result.editDistance;
        editdist = (int)snk->seqlen - match;
    } else {
        result = edlibAlign(snk->seq, snk->seqlen, src->seq, src->seqlen,
                            edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
        editdist = result.editDistance;
    }
    edlibFreeAlignResult(result);

    assert(editdist >= -1);
    assert(editdist <= (int)snk->seqlen || fprintf(stderr, ">%s\n%s\n>%s\n%s\neditdist=%d>seq1len\n", snk->name, snk->seq, src->name, src->seq, editdist));
    // assert(editdist <= (int)src->seqlen || fprintf(stderr, ">%s\n%s\n>%s\n%s\neditdist=%d>srclen\n", seq1->name, seq1>seq, src->name, src->seq, editdist));
    
    if (-1 == editdist) { return 0; }
    if (SIM_ZVAL) {
        
        double mean = snk->seqlen * matchprob;
        double var  = snk->seqlen * matchprob * (1 - matchprob);

        double zval = ((double)(snk->seqlen - editdist) - mean) /  sqrt(var);
        //fprintf(stderr, "seq1=%s,src=%s,mean=%f,var=%f,zval=%f,seqlen1=%d,editdist=%d\n", seq1->name, src->name, mean, var, zval, seq1->seqlen, editdist);
        
        if (zval < SIM_ZVAL / 100.0) { return 0; }
        if (ZVAL_AS_SIM) { return MIN(100*2, 100 + floor(zval)); } // 100 (100*2) are the min (max) scores for a hit by z-value, respectively
    }
    int ret = 100 * (snk->seqlen - editdist) / snk->seqlen;
    return ret;
}

void comp_store_short_words(std::vector<uint64_t> &signs, const seq_t *seq_ptr) {
    if ((int)seq_ptr->seqlen >= (int)SIGN_LENGTH) {
        signs.reserve((int)seq_ptr->seqlen - (int)SIGN_LENGTH + 1);
        uint64_t sign = sign_init(seq_ptr->seq);
        signs.push_back(sign);
        for (unsigned int i = SIGN_LENGTH; i < seq_ptr->seqlen; i += 1) {
            sign = sign_update(sign, seq_ptr->seq[i-SIGN_LENGTH], seq_ptr->seq[i]);
            signs.push_back(sign);
        }
        std::sort(signs.rbegin(), signs.rend());
    }
}

int comp_n_shared_shortwords(std::vector<uint64_t> &src_shortwords, const seq_t *src, const seq_t *snk) {
    if (0 == src_shortwords.size()) {
        comp_store_short_words(src_shortwords, src);
    }
    std::vector<uint64_t> snk_shortwords;
    comp_store_short_words(snk_shortwords, snk);
    
    auto it1 = src_shortwords.begin();
    auto it2 = snk_shortwords.begin();
    int ret = 0;
    while (it1 != src_shortwords.end() && it2 != snk_shortwords.end()) {
        if ((*it1) == (*it2)) {
            ret++;
            it1++;
            it2++;
        } else if ((*it1) < (*it2)) {
            it2++;
        } else if ((*it1) > (*it2)) {
            it1++;
        }
    }
    int ret2 = ret * NUM_SIGNATURES / snk_shortwords.size();
    assert (ret2 <= NUM_SIGNATURES);
    return ret2;
}

void comp_store_words(std::vector<uint32_t> &words, const seq_t *seq_ptr) {
    if ((int)seq_ptr->seqlen >= (int)SIGN_LENGTH) {
        std::vector<uint16_t> signs;
        signs.reserve((int)seq_ptr->seqlen - (int)SIGN_LENGTH + 1);
        uint64_t sign = sign_init(seq_ptr->seq);
        signs.push_back((uint16_t)sign);
        for (unsigned int i = SIGN_LENGTH; i < seq_ptr->seqlen; i += 1) {
            sign = sign_update(sign, seq_ptr->seq[i-SIGN_LENGTH], seq_ptr->seq[i]);
            signs.push_back((uint16_t)sign);
        }
        std::sort(signs.begin(), signs.end());
        unsigned int check = 0;
        uint16_t prevsign = 0;
        uint16_t prevsigncnt = 0;
        for (auto sign : signs) {
            if (sign != prevsign) {
                if (prevsigncnt) { words.push_back(65536 * prevsign + prevsigncnt); check += prevsigncnt; }
                prevsign = sign;
                prevsigncnt = 1;
            } else {
                prevsigncnt++;
            }
        }
        if (prevsigncnt) { words.push_back(65536 * prevsign + prevsigncnt); check += prevsigncnt; }
        assert(check <= seq_ptr->seqlen);
        assert(check < 65536);
    }
}

int comp_n_shared_signatures(std::vector<uint64_t> &src_shortwords,  const seq_t *src, const seq_t *snk) {

#if SHORTWORD
    return comp_n_shared_shortwords(src_shortwords, src, snk);
#elif ENTROHASH
    uint64_t r1 = src->compressedsign  ^ snk->compressedsign;
    uint64_t r2 = src->compressedsign2 ^ snk->compressedsign2;
    uint64_t r3 = src->compressedsign3 ^ snk->compressedsign3;
    uint64_t r4 = src->compressedsign4 ^ snk->compressedsign4;
    uint64_t ra = = r1 & r2 & r3 & r4;
    return __builtin_popcountll(ra);

    int ret1 = __builtin_popcountll(src->compressedsign  ^ snk->compressedsign );
    int ret2 = __builtin_popcountll(src->compressedsign2 ^ snk->compressedsign2);
    int ret3 = __builtin_popcountll(src->compressedsign3 ^ snk->compressedsign3);
    int ret4 = __builtin_popcountll(src->compressedsign4 ^ snk->compressedsign4);
    // assert (ret <= NUM_SIGNATURES);
    return ret1 & ret2 & ret3 & ret4;
    return MIN(MIN(ret2, ret3), ret4);
    return ret1 + ret2 + ret3 + ret4;
    // return MAX(ret - 32, 0);
#elif VARSIGNS 
    unsigned long i = 0;
    unsigned long j = 0;
    const unsigned long ni = seqlen_to_n_varsigns(src->seqlen);
    const unsigned long nj = seqlen_to_n_varsigns(snk->seqlen);
    unsigned long minsigns = (unsigned long) MAX((unsigned int)NUM_SIGNATURES, MIN(src->seqlen, snk->seqlen) / 128);
    const unsigned long imax = MIN(ni, minsigns * 100UL*9999UL / ((unsigned long)LEN_PERC_SRC * 9999UL + 1UL));
    const unsigned long jmax = MIN(nj, minsigns);
    unsigned long ret = 0;
    while (i != imax && j != jmax) {
        if (src->varsigns[i] == snk->varsigns[j]) {
            ret++;
            i++; 
            j++;
        } else if (src->varsigns[i] < snk->varsigns[j]) {
            j++;
        } else {
            i++;
        }
    }
    
    // return ret / (j * snk->seqlen / nj / src->seqlen * ni);
    unsigned long ret2 = (unsigned long)NUM_SIGNATURES * (ret * nj * (unsigned long)src->seqlen) / (ni * (unsigned long)snk->seqlen * (unsigned long)j + 1UL);
    
    // std::cerr << src->name << " and " << snk->name << " share " << ret << " " << ret2 << "/" << NUM_SIGNATURES << " minhash-values out of " << imax << " and " << jmax << std::endl;
    
    return (int)MIN(ret2, (unsigned long)NUM_SIGNATURES);
#else
    // The value 75 is an empirical threshold that works well in practice. O
    // Once a seq A is less than 75% long than another seq B, counting the same number of minhash features does not work.
    if (src->seqlen * 75 > snk->seqlen * 100 || snk->seqlen * 75 > src->seqlen * 100) {
        return comp_n_shared_shortwords(src_shortwords, src, snk); 
    }
    int i = 0;
    int j = 0;
    int ret = 0;
    while (i != NUM_SIGNATURES && j != NUM_SIGNATURES) {
        if (src->signatures[i] == snk->signatures[j]) {
            ret++; 
            i++; 
            j++;
        } else if (src->signatures[i] < snk->signatures[j]) {
            j++;
        } else if (src->signatures[i] > snk->signatures[j]) {
            i++;
        } else {
            abort();
        }
    }
    return ret;
#endif
}

typedef struct {
    seq_t *data;
    uint32_t bufsize;
    uint32_t size;
}
seq_arrlist_t;

seq_arrlist_t seq_arrlist;

int comp_n_shared_words(const std::vector<uint32_t> &src_shortwords, const uint32_t snkidx, uint16_t *snkhcounts, uint32_t *snkseqidxs) {
    seq_t *snk = &seq_arrlist.data[snkidx];
    if ((int)snk->seqlen >= SIGN_LENGTH) {
        uint64_t sign = sign_init(snk->seq);
        if (snkidx == snkseqidxs[sign%65536U]) {
            snkhcounts[sign%65536U]++;
        } else {
            snkhcounts[sign%65536U] = 1;
            snkseqidxs[sign%65536U] = snkidx;
        }
        for (unsigned int i = SIGN_LENGTH; i < snk->seqlen; i += 1) {
            sign = sign_update(sign, snk->seq[i-SIGN_LENGTH], snk->seq[i]);
            if (snkidx == snkseqidxs[sign%65536U]) {
                snkhcounts[sign%65536U]++;
            } else {
                snkhcounts[sign%65536U] = 1;
                snkseqidxs[sign%65536U] = snkidx;
            }
        }
    }
    int same = 0;
    for (uint32_t srcsign : src_shortwords) {
        uint16_t idx = (uint16_t)(srcsign / 65536);
        uint16_t cnt = (uint16_t)(srcsign % 65536);
        // std::cerr << "Cnt = " << cnt << std::endl;
        if (snkidx == snkseqidxs[idx]) {
            same += MIN(snkhcounts[idx], cnt);
        }
    }
    return same;
}



seed_t *seeds;

void seed_add(uint64_t hash, uint32_t seqidx) {
    seed_t *seed = &seeds[hash % DBENTRY_CNT];
    if (seed->size && seed->seqidxs[seed->size-1] == seqidx) { return; }
    if (seed->size == seed->bufsize) {
        seed->bufsize = MAX(seed->bufsize * 2, seed->bufsize + 1);
        seed->seqidxs = (uint32_t*) xrealloc(seed->seqidxs, seed->bufsize * sizeof(uint32_t));
    }
    seed->seqidxs[seed->size] = seqidx;
    seed->size++;
}

void seq_arrlist_init() {
    seq_arrlist.data = (seq_t*) xmalloc(4 * sizeof(seq_t));
    seq_arrlist.bufsize = 4;
    seq_arrlist.size = 0;
}

void seq_arrlist_add(const kseq_t *kseq) {
    if (seq_arrlist.size == seq_arrlist.bufsize) {
        seq_arrlist.bufsize *= 2;
        seq_arrlist.data = (seq_t*)xrealloc(seq_arrlist.data, seq_arrlist.bufsize * sizeof(seq_t));
    }
    assert(kseq->name.l > 0);
    assert(kseq->seq.l > 0);
    seq_arrlist.data[seq_arrlist.size].name = (char*)xmalloc(kseq->name.l + 1);
    seq_arrlist.data[seq_arrlist.size].seq = (char*)xmalloc(kseq->seq.l + 1);
    strcpy(seq_arrlist.data[seq_arrlist.size].name, kseq->name.s);
    // for (size_t i = 0; i <= kseq->seq.l; i++) { seq_arrlist.data[seq_arrlist.size].seq[i] = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][kseq->seq.s[i]]; }
    strcpy(seq_arrlist.data[seq_arrlist.size].seq, kseq->seq.s);
    assert( seq_arrlist.data[seq_arrlist.size].seq[kseq->seq.l] == '\0' || !fprintf(stderr, "The string '%s' is not null-terminated\n", seq_arrlist.data[seq_arrlist.size].seq));
    seq_arrlist.data[seq_arrlist.size].seqlen = kseq->seq.l;
    seq_arrlist.data[seq_arrlist.size].coveredcnt = 0;
    // seq_arrlist.data[seq_arrlist.size].coveringcnt = 0;
    seq_arrlist.size++;
}

void PARAMS_init(const int argc, const char *const *const argv) {
    int nb_cnt = 0;
    int aa_cnt = 0;
    int baseTcnt = 0;
    int baseUcnt = 0;
    for (unsigned int i = 0; i < MIN(seq_arrlist.size, BATCHSIZE_INI); i++) {
        for (unsigned int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
            if (NULL != strchr("ACGTUacgtu", seq_arrlist.data[i].seq[j])) {
                nb_cnt++;
                if (NULL != strchr("Tt", seq_arrlist.data[i].seq[j])) {
                    baseTcnt++;
                }
                if (NULL != strchr("Uu", seq_arrlist.data[i].seq[j])) {
                    baseUcnt++;
                }

            } else {
                aa_cnt++;
            }
        }
    }
    if (aa_cnt * 4 > nb_cnt) {
        SEQTYPE = 1;
    } else {
        if (baseTcnt < baseUcnt) {
            SEQTYPE = 2;
        } else {
            SEQTYPE = 3;
        }
    }
    
    std::vector<int> are_args_parsed(argc);
    std::fill(are_args_parsed.begin(), are_args_parsed.end(), 0);

    for (int i = 1; i+1 < argc; i += 2) {
        int is_arg_parsed = 1;
        if      (!strcmp("--procseqs-order", argv[i])) { /* pass */ } 
        else if (!strcmp("--setcover-depth", argv[i])) { /* pass */ }
        else if (!strcmp("--seqtype", argv[i])) { SEQTYPE = atoi(argv[i+1]); } 
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }

    LEN_PERC_SRC = LEN_PERC_SNK = 80;
    SIM_BASE = 25;
    if (1 == SEQTYPE) {
        SIM_PERC = 50;
        //SEED_N_PER_SEQ = 10;
    } else if (2 == SEQTYPE) {
        SIM_PERC = 70;
        //SEED_N_PER_SEQ = 10;// * 2;
    } else if (3 == SEQTYPE) {
        LEN_PERC_SRC = LEN_PERC_SNK = 0;
        SIM_PERC = 90;
        SIM_BASE = 0;
        //SEED_N_PER_SEQ = 10;// * 2;
    } else {
        std::cerr << "The value of " << SEQTYPE << " is invalid for the parameter --seqtype." << std::endl;
        show_usage(argc, argv);
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            ALPHA_TYPE_TO_CHAR_TO_REDUCED[j][i] = (char)i;
        }
    }
    
    for (int i = 1; i+1 < argc; i += 2) {
        int is_arg_parsed = 1;
        
        if      (!strcmp("--alphasize-seed", argv[i])) { ALPHA_TYPE_TO_SIZE[0] = atoi(argv[i+1]); } 
        else if (!strcmp("--alphasize-sign", argv[i])) { ALPHA_TYPE_TO_SIZE[1] = atoi(argv[i+1]); } 
        else if (!strcmp("--alphasize-sim",  argv[i])) { ALPHA_TYPE_TO_SIZE[2] = atoi(argv[i+1]); } 
        
        else if (!strcmp("--attempt-ini",    argv[i])) { ATTEMPT_INI           = atof(argv[i+1]); } 
        else if (!strcmp("--attempt-inc",    argv[i])) { ATTEMPT_INC           = atof(argv[i+1]); } 
        else if (!strcmp("--attempt-max",    argv[i])) { ATTEMPT_MAX           = atof(argv[i+1]); } 
        
        else if (!strcmp("--cov-src-ada",    argv[i])) { COV_SRC_ADA           = atof(argv[i+1]); }
        else if (!strcmp("--cov-snk-max",    argv[i])) { COV_SNK_MAX           = atoi(argv[i+1]); } 
        
        else if (!strcmp("--dbfilt-minseed", argv[i])) { DBFILT_MINSEED        = atoi(argv[i+1]); } 
        else if (!strcmp("--dbfilt-subsamp", argv[i])) { DBFILT_SUBSAMP        = atoi(argv[i+1]); } 
        // else if (!strcmp("--dbfilt-attempt", argv[i])) { DBFILT_ATTEMPT        = atoi(argv[i+1]); } 
        else if (!strcmp("--dbfilt-timefac", argv[i])) { DBFILT_TIMEFAC        = atoi(argv[i+1]); } 
                
        else if (!strcmp("--len-perc-src",   argv[i])) { LEN_PERC_SRC          = atoi(argv[i+1]); } 
        else if (!strcmp("--len-perc-snk",   argv[i])) { LEN_PERC_SNK          = atoi(argv[i+1]); } 
        
        else if (!strcmp("--seed-n-per-seq", argv[i])) { SEED_N_PER_SEQ        = atoi(argv[i+1]);} 
        else if (!strcmp("--seed-border",    argv[i])) { SEED_BORDER           = atof(argv[i+1]); } 
        else if (!strcmp("--seed-length",    argv[i])) { SEED_LENGTH           = atoi(argv[i+1]); } 
        else if (!strcmp("--sign-chcov-max", argv[i])) { SIGN_CHCOV_MAX        = atoi(argv[i+1]); }
                
        else if (!strcmp("--sim-zval",       argv[i])) { SIM_ZVAL              = atoi(argv[i+1]); } 
        else if (!strcmp("--sim-perc",       argv[i])) { SIM_PERC              = atoi(argv[i+1]); } 
        else if (!strcmp("--sim-base",       argv[i])) { SIM_BASE              = atoi(argv[i+1]); } 
        
        else if (!strcmp("--zval-as-sim",    argv[i])) { ZVAL_AS_SIM           = atoi(argv[i+1]); } 
        
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }
    LEN_PERC_SNK = MAX(LEN_PERC_SNK, SIM_PERC);
    
    if (1 != SEQTYPE) {
        SIGN_LENGTH = (SIM_PERC + 900) / (200 - SIM_PERC); // heuristic values from cd-hit manual page
    } else {
        SIGN_LENGTH = (SIM_PERC + 360) / (150 - SIM_PERC); // heuristic values from cd-hit manual page
    }

    if (3 == SEQTYPE) {
        SIM_DIFF = MIN(SIM_PERC,  20);
    } else {
        SIM_DIFF = 0;
    }

    // ATTEMPT_BASE = 10 + SEED_EVALUE * 5;
    COV_SRC_MAX = 1 + (100 - SIM_PERC) / 10;
    SEED_EVALUE = pow(10, 6 - 0.1 * SIM_PERC); // 10; // 1e6 / pow(10, 0.1 * SIM_PERC);
    SEED_MINCNT = 1 + floor(sqrt((100 - SIM_PERC) * 40)); //  20 + SQUARE((100 - SIM_PERC) / 10);

    for (int i = 1; i+1 < argc; i += 2) {
        int is_arg_parsed = 1;
             if (!strcmp("--attempt-base",   argv[i])) { ATTEMPT_BASE          = atoi(argv[i+1]); } 
        else if (!strcmp("--cov-src-max",    argv[i])) { COV_SRC_MAX           = atoi(argv[i+1]); }
        else if (!strcmp("--seed-evalue",    argv[i])) { SEED_EVALUE           = atof(argv[i+1]); } 
        else if (!strcmp("--seed-mincnt",    argv[i])) { SEED_MINCNT           = atoi(argv[i+1]); }
        else if (!strcmp("--sign-length",    argv[i])) { SIGN_LENGTH           = atoi(argv[i+1]); }
        else if (!strcmp("--sign-cntmin",    argv[i])) { SIGN_CNTMIN           = atoi(argv[i+1]); } 
        else if (!strcmp("--sign-zval",      argv[i])) { SIGN_ZVAL             = atoi(argv[i+1]); }
        else if (!strcmp("--sim-diff",       argv[i])) { SIM_DIFF              = atoi(argv[i+1]); }
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }

    IDXENTRY_ITMAX = floor(100 * (10 + SEED_EVALUE));
    
    for (int i = 1; i <= NUM_SIGNATURES; i++) {
        double p = pow((1 - DBL_EPSILON) * SIM_PERC / 100, SIGN_LENGTH);
        double sd = sqrt((i) * p * (1-p)) * NUM_SIGNATURES / i;
        double mean = i * p;
        SIGN_SD_CNTMINS[i] = MAX(ceil(mean - sd * SIGN_ZVAL), 0); // three standard deviations below the normal distribution of true positives.
    }
    
    for (int i = 1; i+1 < argc; i += 2) {
        int is_arg_parsed = 1;
             if (!strcmp("--idxentry-itmax", argv[i])) { IDXENTRY_ITMAX        = atoi(argv[i+1]); } 
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }
    
    for (int i = 1; i < argc; i++) {
        if (0 == are_args_parsed[i]) {
            std::cerr << "Cannot parse command-line parameter " << argv[i] << " at position " << i << std::endl;
            show_usage(argc, argv);
        }
        assert(1 == are_args_parsed[i] 
               || !(std::cerr << "Parsed command-line parameter " << argv[i] << " at position " << i << " " 
                              << are_args_parsed[i] << " times!" << std::endl));
    }
    
    if (1 == SEQTYPE) {
        for (int t = 0; t < 3; t++) {
            if (10 == ALPHA_TYPE_TO_SIZE[t]) {
                alphareduce("DENQ", t);
                alphareduce("FWY", t);
                alphareduce("ILMV", t);
                alphareduce("KR", t);
                alphareduce("ST", t);
            }
            if (15 == ALPHA_TYPE_TO_SIZE[t]) {
                alphareduce("FY", t);
                alphareduce("ILMV", t);
                alphareduce("KR", t);
            }
        }
    }
    showparams();
}

void seq_longword_init(seq_t *const seq_ptr, int idx) {
    for (int seedlen = SEED_LENGTH; seedlen <= SEED_LENGTH + 1; seedlen += 1) {
        if ((int)seq_ptr->seqlen >= seedlen) {
            unsigned int interseed_gap = MAX((unsigned int)1, (unsigned int)((1 + (int)seq_ptr->seqlen - seedlen) / SEED_MINCNT));
            uint64_t hash1 = hash_init(seq_ptr->seq, seedlen);
            if (seedlen == SEED_LENGTH) { seed_add(hash1, idx); }
            for (int i = seedlen; i < (int)seq_ptr->seqlen; i += 1) {
                hash1 = hash_update(hash1, seq_ptr->seq[i-seedlen], seq_ptr->seq[i], seedlen);
                bool is_highlen = ((100 * (1 + i - seedlen) / (1 + (int)seq_ptr->seqlen - seedlen)) >= 100 - SEED_HIGH_LEN_PERC);
                if ((0 == (1 + i - seedlen) % interseed_gap) && ((SEED_LENGTH == seedlen && !is_highlen) || (SEED_LENGTH + 1 == seedlen && is_highlen))) {
                    seed_add(hash1, idx);
                }
            }
        }
    }
}

void seq_signatures_init(seq_t *const seq_ptr) {
#if !SHORTWORD
#if ENTROHASH  
    seq_ptr->compressedsign = 0;
    seq_ptr->compressedsign2 = 0;
    seq_ptr->compressedsign3 = 0;
    seq_ptr->compressedsign4 = 0;
#elif VARSIGNS
    seq_ptr->varsigns = NULL;
#else
    memset(seq_ptr->signatures, 0, NUM_SIGNATURES * sizeof(uint16_t));
#endif
    if ((int)seq_ptr->seqlen >= (int)SIGN_LENGTH) {
        std::vector<uint64_t> signs;
        signs.reserve((int)seq_ptr->seqlen - (int)SIGN_LENGTH + 1);
        uint64_t sign = sign_init(seq_ptr->seq);
        signs.push_back(sign);
        for (unsigned int i = SIGN_LENGTH; i < seq_ptr->seqlen; i += 1) {
            sign = sign_update(sign, seq_ptr->seq[i-SIGN_LENGTH], seq_ptr->seq[i]);
            signs.push_back(sign);
        }
        std::sort(signs.rbegin(), signs.rend());
#if ENTROHASH
        auto it = signs.begin();
        uint64_t stepval = (SIGN_MOD + 256/2) / 256;
        uint64_t currval = SIGN_MOD;
        //std::vector<uint64_t> signvec;
        for (int k = 256 - 1; k >= 0; --k) {
            while (it != signs.end() && (*it) > currval) { ++it; }
            currval -= stepval;
            if (it == signs.end()) {break;}
            if ((k%4) == 0) { seq_ptr->compressedsign  |= ( (__builtin_popcountll (*it) & 0x1UL)        << (k/4)); }
            if ((k%4) == 1) { seq_ptr->compressedsign2 |= ( (                     (*it) & 0x1UL)        << (k/4)); }
            if ((k%4) == 2) { seq_ptr->compressedsign3 |= ( (                    ((*it) & 0x2UL) / 2UL) << (k/4)); }
            if ((k%4) == 3) { seq_ptr->compressedsign4 |= ( (                    ((*it) & 0x4UL) / 4UL) << (k/4)); }
            //signvec.push_back(*it);
            // std::cerr << "I picked sign " << (*it) << " for threshold " << currval << std::endl;
        }
        assert(UINT64_MAX - 256/2 < currval || currval < 256/2 || it == signs.end());
        if (!seq_ptr->compressedsign && signs.size() > 40) { 
            std::cerr << seq_ptr->name << " may have wrong compressed-sign" << std::endl; 
            //for (auto sign : signvec) {
            //    std::cerr << "picked sign = " << sign << std::endl;
            //}
            for (auto sign : signs) {
                std::cerr << "sign = " << sign << std::endl;
            }
            abort();
        }
#elif VARSIGNS
        unsigned int n_varsigns = seqlen_to_n_varsigns(seq_ptr->seqlen);
        std::vector<uint16_t> signatures;
        signatures.reserve(n_varsigns);
        for (auto sign : signs) {
            signatures.push_back((uint16_t)sign);
            if (signatures.size() == n_varsigns) { break; }
        }
        std::sort(signatures.rbegin(), signatures.rend());
        seq_ptr->varsigns = (uint16_t*) malloc(n_varsigns * sizeof(uint16_t));
        unsigned int j = 0;
        for (auto sign : signatures) {
            seq_ptr->varsigns[j] = sign;
            j++;
        }
        assert(n_varsigns == j);
#else    
        int j = 0;
        std::vector<uint16_t> signatures;
        signatures.reserve(NUM_SIGNATURES);
        for (auto sign : signs) {
            signatures.push_back((uint16_t)sign);
            j++;
            if (NUM_SIGNATURES == j) { break; }
        }
        std::sort(signatures.rbegin(), signatures.rend());
        j = 0;
        for (auto sign : signatures) {
            seq_ptr->signatures[j] = sign;
            j++;
        }
        assert(j <= NUM_SIGNATURES);
    
#endif
    }
    /*
    for (int i = 0; i < NUM_SIGNATURES; i++) {
        seq_ptr->compressedsign |= (0x1L << ((seq_ptr->signatures[i] * 23) % 64));
    }
    */
#endif
}

void seed_cov(const seed_t *seed, const uint32_t coveringidx, std::set<uint32_t> & visited) {
    if (IDXENTRY_ITMAX < seed->size && seed->seqidxs[IDXENTRY_ITMAX] <= coveringidx) {
        return;
    }
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    for (unsigned int i = 0; i < seed->size; i++) { 
        unsigned int coveredidx = seed->seqidxs[i];
        if (!fail_len_cutoff(coveringseq, &seq_arrlist.data[coveredidx])) {
            visited.insert(coveredidx);
        }
    }
}

void set_seedsize_histogram(int seedsize_histogram[]) {
    memset(seedsize_histogram, 0, 512 * sizeof(int));
    for (unsigned int i = 0 ; i < DBENTRY_CNT; i++) {
        int histidx = floor(log(seeds[i].size + 1 + DBL_EPSILON) * 10);
        assert(histidx < 512 || !fprintf(stderr, "dbentry %d has seed count %d which is too high!\n", i, seeds[i].size));
        seedsize_histogram[histidx]++;
    }
}

void print_seedsize_histogram(const int hist1[]) {

    std::cerr << "Seedsize histogram with lower/upper delimiters and frequency before/after filter as columns:" << std::endl;
    for (int i = 0; i < 512; i++) {
        if (hist1[i] > 0)
        {
            double lowerlim = (exp(((double)i  ) / 10) - 1 - DBL_EPSILON);
            double upperlim = (exp(((double)i+1) / 10) - 1 - DBL_EPSILON); 
            std::cerr << "\t" << lowerlim << "\t" << upperlim << "\t" << hist1[i] << std::setprecision(15) << std::endl;
        }
    }
    uint64_t n_pairwise_cmp_1 = 0;
    uint64_t n_pairwise_cmp_2 = 0;
    for (unsigned int i = 0 ; i < DBENTRY_CNT; i++) {
        n_pairwise_cmp_1 += (uint64_t)seeds[i].size * (uint64_t)seeds[i].size;
        n_pairwise_cmp_2 += (uint64_t)seeds[i].size * (uint64_t)MIN(seeds[i].size, IDXENTRY_ITMAX);
    }
    std::cerr << "Comparison-counts-totaled: unlimited = " << n_pairwise_cmp_1 << " and limited = " << n_pairwise_cmp_2 << std::endl;
    std::cerr << "Comparison-counts-per-seq: unlimited = " << n_pairwise_cmp_1 / seq_arrlist.size << " and limited = " << n_pairwise_cmp_2 / seq_arrlist.size << std::endl;
}

int main(const int argc, const char *const *const argv) {
    assert(argc > 0);
    if (1 < argc && (!strcmp("--help", argv[1]) || !strcmp("--usage", argv[1]))) {
        if (2 < argc) { std::cerr << "The first command line parameter \"" << argv[1] << "\" must not be followed by any other parameter." << std::endl; }
        show_usage(argc, argv);
    }
    time_t begtime, endtime;

    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;
#ifdef DESC
    std::cerr << "DESC = " << DESC << std::endl;
#endif
    std::cerr << "For usage and help, please enter either one of the following commands:" << std::endl;
    std::cerr << "\t" << argv[0] << " --help" << std::endl;
    std::cerr << "\t" << argv[0] << " --usage" << std::endl;
    std::cerr << "Reading fasta sequences from STDIN" << std::endl; 

    std::set<int> printthresholds;
    for (int i = 0; i < 300; i++) {
        double thres = (double)((i+1)*(BATCHSIZE_INI+1)) * pow(1.05, (double)i);
        if (thres * 1.01 < (double)(INT_MAX)) { printthresholds.insert((int)thres); }
    }

    seq_arrlist_init();
    
    kseq_t *kseq = kseq_init(fileno(stdin));
    unsigned int i = 0;
    time(&begtime);
    while ( kseq_read(kseq) >= 0 ) {
        seq_arrlist_add(kseq);
        if (0 == ((i+1) & i)) {
            time(&endtime);
            fprintf(stderr, "Read %d sequences in %.f seconds.\n", i+1, difftime(endtime, begtime));
        }
        i++;
        if (BATCHSIZE_INI == i) {
            PARAMS_init(argc, argv);
        }
    }
    kseq_destroy(kseq);
    assert(0 < i);
    if (i < BATCHSIZE_INI) { PARAMS_init(argc, argv); }
    
    // reinitialize some vars
    uint64_t num_residues = 0;
    for (unsigned int i = 0 ; i < seq_arrlist.size; i++) {
        num_residues += seq_arrlist.data[i].seqlen;
    }
    DBENTRY_CNT = seq_arrlist.size * SEED_MINCNT / 3; // SEED_N_PER_SEQ;
    std::cerr << "Derived paramters: " << std::endl;
    std::cerr << "\tNUM_SEQS = " << seq_arrlist.size << std::endl;
    std::cerr << "\tNUM_RESIDUES = " << num_residues << std::endl;
    std::cerr << "\tDBENTRY_CNT = " << DBENTRY_CNT << std::endl;
    
    double ch_to_freq_weighted_by_seqlen[256];
    memset(ch_to_freq_weighted_by_seqlen, 0, 256 * sizeof(double));
    // subsample 1000 * 1000 letters in the db
    unsigned int randstate = seq_arrlist.size;
    for (unsigned int i = 0 ; i < 1000 * 1000; i++) {
        unsigned int seqidx = rand_r(&randstate) % seq_arrlist.size;
        unsigned int residx = rand_r(&randstate) % seq_arrlist.data[seqidx].seqlen;
        char ch = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][(int)seq_arrlist.data[seqidx].seq[residx]];
        ch_to_freq_weighted_by_seqlen[(int)ch] += (double)seq_arrlist.data[seqidx].seqlen;
    }
    double ch_totcnt = 0;
    for (int i = 0; i < 256; i++) {
        ch_totcnt += ch_to_freq_weighted_by_seqlen[i];
    }
    double ch_entropy = 0;
    std::cerr << "Letter frequency distribution:" << std::endl;
    for (int i = 0; i < 256; i++) {
        double prob = ch_to_freq_weighted_by_seqlen[i] / ch_totcnt;
        if (0 < prob) { 
            ch_entropy -= prob * log(prob); 
            if (32 <= i && i < 127) {
                std::cerr << "\t" << i << "\t" << prob << "\t" << (char)i << std::endl; // print printable ascii
            } else {
                std::cerr << "\t" << i << "\t" << prob << "\t" << "not-printable" << std::endl; 
            }
        }
    }

    double INFO_PER_LETTER = exp(ch_entropy); 
    // 0.33 is used to normalize whole-protein Shannon-information entropy to about 1 bit per position // https://www.ncbi.nlm.nih.gov/pubmed/8804598
    double SHANNON_INFO_PER_LETTER = pow(INFO_PER_LETTER, 0.33); // estimated
    std::cerr << "INFO_PER_LETTER = " << INFO_PER_LETTER << " , SHANNON_INFO_PER_LETTER = " << SHANNON_INFO_PER_LETTER << std::endl;
    std::cerr << "Recommended similarity-threshold for detecting homology = " << 1 / SHANNON_INFO_PER_LETTER << " , actual similarity-threshold = " << SIM_PERC << std::endl;

    if (SEED_EVALUE > 0) {
        // double seedlen_multi = (double)(MAX(SIM_PERC, 100 / SHANNON_INFO_PER_LETTER)) / (double)(100 / SHANNON_INFO_PER_LETTER);
        double seedlen_fract = log((double)num_residues / SEED_EVALUE + INFO_PER_LETTER) / log(INFO_PER_LETTER);
        // seedlen_fract = MAX(seedlen_fract, seedlen_fract * seedlen_multi - 3);
        seedlen_fract = MIN(MAX(seedlen_fract, 3), 30);
        int    seedlen_floor = (int)floor(seedlen_fract);
        // double seedlen_diff1 = seedlen_fract - seedlen_floor;
        SEED_LENGTH = seedlen_floor;
        SEED_HIGH_LEN_PERC = floor(100 * (seedlen_fract - seedlen_floor));
        //SEED_MINCNT = (int)floor((120 - SIM_PERC) 
                      // / pow(SHANNON_INFO_PER_LETTER, seedlen_diff1 / seedlen_multi) 
                      // / (1 + (seedlen_fract - seedlen_floor) * (100 - SIM_PERC) / SIM_PERC)
        //              * ((1 != SEQTYPE) ? 2 : 1));
        
        std::cerr << "Command-line parameter values after adjustment with SEED_EVALUE = " << SEED_EVALUE << ":" << std::endl;
        std::cerr << "\tadjusted SEED_LENGTH = " << SEED_LENGTH << " from " << seedlen_fract << " with SEED_HIGH_LEN_PERC = " << SEED_HIGH_LEN_PERC <<  std::endl;
        // std::cerr << "\tadjusted SEED_MINCNT = " << SEED_MINCNT << std::endl;
    } else {
        std::cerr << "No adjustment made because SEED_EVALUE = " << SEED_EVALUE << std::endl;
    }

    hash_sign_INIT();

    seeds = (seed_t*) xmalloc(DBENTRY_CNT * sizeof(seed_t));
    memset(seeds, 0, DBENTRY_CNT * sizeof(seed_t));
    
    time(&begtime);
    for (unsigned int i = 0 ; i < seq_arrlist.size; i++) {
        seq_longword_init(&seq_arrlist.data[i], i);
        if (0 == ((i+1) & i)) {
            time(&endtime);
            fprintf(stderr, "Indexed %d sequences in %.f seconds.\n", i+1, difftime(endtime, begtime));
        }
    }
    
    time(&begtime);
    #pragma omp parallel for schedule(dynamic, 9999*10)
    for (unsigned int i = 0 ; i < seq_arrlist.size; i++) {
        seq_signatures_init(&seq_arrlist.data[i]);
        for (unsigned int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
            seq_arrlist.data[i].seq[j] = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][(int)seq_arrlist.data[i].seq[j]];
        }
    }
    time(&endtime);
    std::cerr << "Initialized minhash and reduce alphabet in " << difftime(endtime, begtime) << " seconds." << std::endl;

    int seedsize_hist1[512];
    set_seedsize_histogram(seedsize_hist1); 
    print_seedsize_histogram(seedsize_hist1); 

    printf("%d %d\n", seq_arrlist.size, seq_arrlist.size);
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> coveredarr(BATCHSIZE_INI);
    std::fill(coveredarr.begin(), coveredarr.end(), std::vector<std::pair<uint32_t, uint8_t>>(0));
    unsigned int batchsize = BATCHSIZE_INI;
    unsigned int cov_src_max = COV_SRC_MAX;
    //double clusize_sum = 0;
    //unsigned int clusize_cnt = 0;
    time(&begtime);
    unsigned long nextprint = 0;

    for (unsigned int iter = 0; iter < seq_arrlist.size;) {
        unsigned int itermax = MIN(iter+batchsize, seq_arrlist.size);

        #pragma omp parallel for schedule(dynamic, 1)
        for (unsigned int i = iter; i < itermax; i++)
        {
            std::set<uint32_t> visited;
            int filteredcnt = 0;
            int distcompcnt = 0;
            int max_attempts = 0;
            int max_attempts_arg = 0;
            int nsigns = -1;
            int max_attempts_nsigns = -1;
            if (seq_arrlist.data[i].coveredcnt < cov_src_max && (int)(seq_arrlist.data[i].seqlen) >= (int)SEED_LENGTH) {
                double matchprob = 0;
                double chprobs[256];
                memset(chprobs, 0, 256 * sizeof(double));
                for (unsigned int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
                    chprobs[(int)seq_arrlist.data[i].seq[j]] += 1.0 / seq_arrlist.data[i].seqlen;
                }
                if (SIM_ZVAL) {
                    for (int j = 0; j < 256; j++) {
                        matchprob += chprobs[j] * chprobs[j];
                    }
                }
                std::vector<uint64_t> src_shortwords;
                
                std::vector<uint32_t> srcwords;
                comp_store_words(srcwords, &seq_arrlist.data[i]);
                uint16_t snkhcounts[65536]; 
                uint32_t snkseqidxs[65536];
                memset(snkhcounts, 0, sizeof(*snkhcounts));
                memset(snkseqidxs, 0, sizeof(*snkseqidxs));
                for (int seedlen = SEED_LENGTH; seedlen <= SEED_LENGTH + 1; seedlen += 1) {
                    if ((int)seq_arrlist.data[i].seqlen >= seedlen) {
                        uint64_t hash = hash_init(seq_arrlist.data[i].seq, seedlen);
                        seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                        for (unsigned int j = seedlen; j < seq_arrlist.data[i].seqlen; j++) {
                            hash = hash_update(hash, seq_arrlist.data[i].seq[j-seedlen], seq_arrlist.data[i].seq[j], seedlen);
                            seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                        }
                    }
                }
#if !INDEXWORD
                int NUMSIGNS = NUM_SIGNATURES;
#else
                int NUMSIGNS = seq_arrlist.data[i].seqlen - SIGN_LENGTH + 1;
#endif
                std::vector<std::vector<uint32_t>> nsharedsigns_to_coveredidxs_vec(NUMSIGNS + 1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    if ((i != coveredidx)) {
                        //int n_shared_signatures = comp_n_shared_signatures(src_shortwords, &seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
                        int n_shared_signatures = comp_n_shared_words(srcwords, coveredidx, 
                                                                      snkhcounts, snkseqidxs);
                        assert (n_shared_signatures <= NUMSIGNS || !fprintf(stderr, "n_shared_signatures=%d,NUMSIGNS=%d,seq1=%s,seq2=%s\n", 
                                n_shared_signatures, NUMSIGNS, seq_arrlist.data[i].name, seq_arrlist.data[coveredidx].name));
                        unsigned int norm_signs = ((unsigned int)n_shared_signatures) * seq_arrlist.data[i].seqlen / seq_arrlist.data[coveredidx].seqlen;
                        norm_signs = MIN(norm_signs, (unsigned int) NUMSIGNS);
                        nsharedsigns_to_coveredidxs_vec.at(norm_signs).push_back(coveredidx);
                    }
                }
                for (nsigns = NUMSIGNS; nsigns >= SIGN_CNTMIN; nsigns--) {
                    filteredcnt += nsharedsigns_to_coveredidxs_vec.at(nsigns).size();
                }
                int attempts    = (int)ceil(ATTEMPT_INI * filteredcnt + ATTEMPT_BASE);
                int attempt_inc = (int)ceil(ATTEMPT_INC * filteredcnt + ATTEMPT_BASE);
                int attempt_max = (int)ceil(ATTEMPT_MAX * filteredcnt + ATTEMPT_BASE);
                
                int max_attempts = attempts;
                for (nsigns = NUMSIGNS; nsigns >= SIGN_CNTMIN && attempts > 0 && attempts > max_attempts - attempt_max; nsigns--) {
                    for (auto coveredidx : nsharedsigns_to_coveredidxs_vec.at(nsigns)) {
                        seq_t *coveringseq = &seq_arrlist.data[i];
                        seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                        if (coveredseq->coveredcnt < COV_SNK_MAX) {
                            int sim = comp_perc_seq_sim_editdist(coveringseq, coveredseq, matchprob, chprobs);

                            if (sim >= SIM_PERC) {
                                coveredarr[i-iter].push_back(std::make_pair(coveredidx, sim));
                                attempts += attempt_inc;
                                if (attempts > max_attempts) {
                                    max_attempts = attempts;
                                    max_attempts_arg = distcompcnt + 1;
                                    max_attempts_nsigns = nsigns;
                                }
                            } else {
                                attempts -= 1;
                            }
                            distcompcnt++;
                        }
                        if (!(attempts > 0 && attempts > max_attempts - attempt_max)) { break; }
                    }
                }
                
                std::sort(coveredarr[i-iter].begin(), coveredarr[i-iter].end());
            }
            if (nextprint == i) {
                time(&endtime);
                std::cerr << "In " << difftime(endtime, begtime) << " secs processed " << i+1 << " seqs\t"
                          << "h: " << coveredarr[i-iter].size() << " " << distcompcnt << " " << filteredcnt << " " << visited.size() << "\t" 
                          << "max-att: " << max_attempts << " at " << max_attempts_arg << " sign " << max_attempts_nsigns << " " << nsigns << "\t" 
                          << "coveredcnt: " << seq_arrlist.data[i].coveredcnt << " batchsize " << batchsize << " seqlen " << seq_arrlist.data[i].seqlen << "\n";
                nextprint = MAX(nextprint * 21 / 20, itermax);
            }
        }
#if 0
        if (COV_SRC_ADA > 0) {
            for (unsigned int i = iter; i < MIN(iter + batchsize, itermax); i++) {
                unsigned int clustersize = (coveredarr[i-iter].size() ? coveredarr[i-iter].size() : seq_arrlist.data[i].coveringcnt) + 1;
                clusize_sum += 1 + log(clustersize);
            }
            unsigned int clustercnt = MIN(iter + batchsize, itermax) - iter;
            clusize_cnt += clustercnt;
            cov_src_max = ceil((COV_SRC_MAX + clusize_sum * COV_SRC_ADA) / (1 + clusize_cnt * COV_SRC_ADA));
        }
#endif
        unsigned long maxcov = 1;
        unsigned long extracov = 1;
        for (unsigned int i = iter; i < itermax; i++) {
            maxcov = MAX(maxcov, coveredarr[i-iter].size());
            if (0 < coveredarr[i-iter].size() && COV_SRC_MAX <= seq_arrlist.data[i].coveredcnt) {
                extracov += coveredarr[i-iter].size();
            }
            printf("100 %d", i+1);
            for (auto adj : coveredarr[i-iter]) {
                seq_arrlist.data[adj.first].coveredcnt++;
                //seq_arrlist.data[adj.first].coveringcnt = coveredarr[i-iter].size();
                assert(adj.second <= 100 || ZVAL_AS_SIM);
                printf(" %d %u", (int)adj.second, adj.first+1);
            }
            printf("\n");
            coveredarr[i-iter].clear();
        }
        assert(coveredarr.size() >= batchsize);
        iter += batchsize;
        
        maxcov *= BATCH_COVRAT;
        if (maxcov >= extracov) {
            batchsize += MIN(batchsize, BATCHSIZE_INC * maxcov/extracov);
            while (coveredarr.size() < batchsize) {
                coveredarr.push_back(std::vector<std::pair<uint32_t, uint8_t>>(0));
            }
        } else {
            batchsize -= MIN(batchsize / 2, BATCHSIZE_INC * extracov/maxcov);
        }
    }
}

