#include "edlib.h"
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

#define ENTROHASH 0
#define ORDER_AWARE_FILT 1

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

#define MIN(a, b) (((a) < (b) ? (a) : (b)))
#define MAX(a, b) (((a) > (b) ? (a) : (b)))

const auto SQUARE(const auto v) { return v * v; }
void SWAP(auto &a, auto &b) {auto tmp = a; a = b; b = tmp; }

#ifndef NUM_SIGNATURES
#if ENTROHASH
#define NUM_SIGNATURES (256)
#else
#define NUM_SIGNATURES (32)
#endif
#endif

KSEQ_INIT(int, read)

// constants

const unsigned int BATCHSIZE_INI = 256; // 15999;
const unsigned int BATCHSIZE_INC = 8; // 2; // 15999;

const uint64_t PRIME_BASE = 48271L;
const uint64_t SIGN_BASE  = 48271L; 
const uint64_t PRIME_MOD  = (0x1L << 31L) - 1L; 
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L; 

// variables that are initialized once from constants and other variables

uint64_t SEED_POWER = 0; 
uint64_t SIGN_POWER = 0; 
char ALPHA_TYPE_TO_CHAR_TO_REDUCED[3][256]; // derived from ALPHA_TYPE_TO_SIZE

// variables that are initialized from command line args

int ALPHA_TYPE_TO_SIZE[] = {10, 10, 10};

int ATTEMPT_INI = 50; //50;
int ATTEMPT_INC = 50; //50;
int ATTEMPT_MAX = 50; //50;

uint64_t CHAR_PER_SEED = 40;

unsigned int COV_SRC_MAX = 8; // 5;
unsigned int COV_SNK_MAX = INT_MAX;

unsigned int DBFILT_MINSEED = 1000; // 1000*1000; lower -> more filtering, more time saving later
int DBFILT_SUBSAMP = 50*50; // 800; // lower -> less filtering accuracy, less time
// int DBFILT_ATTEMPT = 50; // 10; // lower -> less filtering accuracy, less time // not useful because same as ATTEMPT_* 
int DBFILT_TRUEHIT = (50*50-1)/100; // 5; // lower -> less filtering accuracy

int IS_INPUT_NUC = -1; // guessed

int LEN_PERC_SRC = -1;
int LEN_PERC_SNK = -1;

int SEED_EVALUE = 1;
int SEED_LENGTH = 10; // can be overriden after determination of db size 
int SEED_MAXGAP = INT_MAX; // can be overriden after determination of db size
int SEED_MINCNT = 10; // can be overriden after determination of db size

int SIGN_LENGTH = -1;
int SIGN_SHARED_CNT_MIN = -1; 

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
    
    std::cerr << "\tCHAR_PER_SEED = "  << CHAR_PER_SEED << std::endl;
    
    std::cerr << "\tCOV_SRC_MAX = " << COV_SRC_MAX  << std::endl;
    std::cerr << "\tCOV_SNK_MAX = " << COV_SNK_MAX << std::endl;
    
    std::cerr << "\tDBFILT_MINSEED = " << DBFILT_MINSEED << std::endl;
    std::cerr << "\tDBFILT_SUBSAMP = " << DBFILT_SUBSAMP << std::endl; 
    std::cerr << "\tDBFILT_TRUEHIT = " << DBFILT_TRUEHIT << std::endl;
    //std::cerr << "\tDBFILT_ATTEMPT = " << DBFILT_ATTEMPT << std::endl;

    std::cerr << "\tIS_INPUT_NUC = " << IS_INPUT_NUC << std::endl;
    
    std::cerr << "\tLEN_PERC_SRC = " << (int)LEN_PERC_SRC << std::endl;
    std::cerr << "\tLEN_PERC_SNK = " << (int)LEN_PERC_SNK << std::endl;

    std::cerr << "\tSEED_EVALUE = " << SEED_EVALUE << std::endl;
    std::cerr << "\tSEED_LENGTH = " << SEED_LENGTH << std::endl;
    std::cerr << "\tSEED_MAXGAP = " << SEED_MAXGAP << std::endl;
    std::cerr << "\tSEED_MINCNT = " << SEED_MINCNT << std::endl;

    std::cerr << "\tSIGN_LENGTH = " << SIGN_LENGTH    << std::endl;
    std::cerr << "\tSIGN_SHARED_CNT_MIN = " << SIGN_SHARED_CNT_MIN      << std::endl;
    
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
    std::cerr << "  --alphasize-sim\t: alphabet size for computation of sequence similarity ["  << ALPHA_TYPE_TO_SIZE[2] << "]" << std::endl;
    
    std::cerr << "  --attempt-ini\t: initial number of attempts. ["                             << ATTEMPT_INI << "]" << std::endl;
    std::cerr << "  --attempt-inc\t: number of attempts incremented per true positive hits. ["  << ATTEMPT_INC << "]" << std::endl;
    std::cerr << "  --attempt-max\t: number of attempts capped at this maximum value. ["        << ATTEMPT_MAX << "]" << std::endl;

    std::cerr << "  --char-per-seed\t: number of residues covered by one seed [" << CHAR_PER_SEED << "]" << std::endl;
    
    std::cerr << "  --cov-snk-max\t: max number of times that the covered sequence can be covered.["     << COV_SNK_MAX   << "]" << std::endl;
    std::cerr << "  --cov-src-max\t: max number of times that the covering sequence can be coverered. [" << COV_SRC_MAX   << "]" << std::endl;
    
    std::cerr << "  --dbfilt-minseed\t: minimum number of times a seed occurs to trigger seed pruning.["           << DBFILT_MINSEED << "]" << std::endl;
    std::cerr << "  --dbfilt-subsamp\t: number of sequence pairs (SP) subsampled for seed pruning .["              << DBFILT_SUBSAMP << "]" << std::endl;
    // std::cerr << "  --dbfilt-attempt\t: number of attempts on most likely (minhash) similar SP for seed pruning.[" << DBFILT_ATTEMPT << "]" << std::endl;
    std::cerr << "  --dbfilt-truepos\t: maximum number of true posive SP hits to trigger pruning.["                << DBFILT_TRUEHIT << "]" << std::endl;
    
    std::cerr << "  --is-input-nuc\t: 0, 1, and 2 mean input is protein, RNA, and DNA, respectively (auto detected). [" << IS_INPUT_NUC << "]" << std::endl;

    std::cerr << "  --len-perc-src\t: A cannot cover B if len(A) * len-perc-src > len(B) * 100. [" << LEN_PERC_SRC << "]" << std::endl; 
    std::cerr << "  --len-perc-snk\t: A cannot cover B if len(B) * len-perc-snk > len(A) * 100. [" << LEN_PERC_SNK << "]" << std::endl; 

    std::cerr << "  --seed-evalue\t: evalue for seed hit. 0 means do not use this parameter ["                   << SEED_EVALUE << "]" << std::endl;
    std::cerr << "  --seed-length\t: length of an indexed seed. Overwritten by nonzero seed-evalue. ["           << SEED_LENGTH << "]" << std::endl;
    std::cerr << "  --seed-maxgap\t: max number of residues between consecutive seeds (not used anymore). ["     << SEED_MAXGAP << "]" << std::endl;
    std::cerr << "  --seed-mincnt\t: minimum number of seeds per sequence. Overwritten by nonzero seed-evalue [" << SEED_MINCNT << "]" << std::endl;
    
    std::cerr << "  --sign-length\t: length of k-mers for computing minhash values. ["                      << SIGN_LENGTH         << "]" << std::endl;
    std::cerr << "  --sign-shared-cnt-min\t: minimum number of minhash values to trigger sequence search [" << SIGN_SHARED_CNT_MIN << "]" << std::endl; 

    std::cerr << "  --sim-zval\t: percent of standard deviations above similarity by chance for a sequence to cover another ["   << SIM_ZVAL << "]" << std::endl;
    std::cerr << "            \t: the effective value is 0.2 to 0.4 times less than this set value due to normal approximation of skewed binomial." << std::endl;
    std::cerr << "  --sim-base\t: A covers B only if (len(B) - edit-distance(A, B)) >= sim_perc * len(B) + sim_base [sim_base="  << SIM_BASE << "]" << std::endl;
    std::cerr << "  --sim-perc\t:   where extra residues at both ends of A are not penalized. [sim_perc="                        << SIM_PERC << "]" << std::endl;
    std::cerr << "  --sim-diff\t: percent of deletions in the covering sequence that are ignored, 0 means skip alignment ["      << SIM_DIFF << "]" << std::endl;

    std::cerr << "  --zval-as-sim\t: used the sim-zval as similarity threshold [" << ZVAL_AS_SIM   << "]" << std::endl;

    std::cerr << "Note: default value of -1 means dependence to other parameters or to the input." << std::endl;
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
    SEED_POWER = 1;
    for (i = 0; i < SEED_LENGTH; i++) {
        SEED_POWER = (SEED_POWER * PRIME_BASE) % PRIME_MOD;
    }
    SIGN_POWER = 1;
    for (i = 0; i < SIGN_LENGTH; i++) { 
        SIGN_POWER = (SIGN_POWER * SIGN_BASE) % SIGN_MOD;
    }
}

const uint64_t hash_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < SEED_LENGTH; i++) {
        ret *= PRIME_BASE;
        ret += (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)beg[i]];
        ret %= PRIME_MOD; 
    }
    return ret;
}

const uint64_t hash_update(uint64_t hash, char prv, char nxt) {
    hash = hash * PRIME_BASE + (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)nxt];
    hash += PRIME_MOD;
    hash -= (SEED_POWER * (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][(int)prv]) % PRIME_MOD;
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
#else
    uint16_t signatures[NUM_SIGNATURES];
#endif
    uint32_t coveredcnt;
    uint32_t seqlen;
}
seq_t; // 16+16*4 bytes +++

const bool fail_len_cutoff(const seq_t *src, const seq_t *snk) {
    return src->seqlen * LEN_PERC_SRC > snk->seqlen * 100 || snk->seqlen * LEN_PERC_SNK > src->seqlen * 100;
}

static const int comp_perc_seq_sim_editdist(const seq_t *src, const seq_t *snk, const double matchprob, const double *chprobs) {
    
    int maxEditDist = snk->seqlen - ceil((1 - DBL_EPSILON) * (sqrt(SQUARE((double)SIM_BASE) + SQUARE((double)(snk->seqlen * (SIM_PERC - SIM_DIFF)) / 100.0))));
    if (maxEditDist < 0) { return 0; }
    
    unsigned int i;
    for (i = 0; i < snk->seqlen; i++) {
        if (chprobs[(int)snk->seq[i]]) { break; }
    }
    if (snk->seqlen == i) { return 0; }

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
        if (ZVAL_AS_SIM) { return MIN(200, 100 + floor(zval)); }
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

int __attribute__((noinline)) comp_n_shared_shortwords(std::vector<uint64_t> &src_shortwords, const seq_t *src, const seq_t *snk) {
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
    int ret2 = ret * 32 / snk_shortwords.size();
    assert (ret2 <= 32);
    return ret2;
}

int comp_n_shared_signatures(std::vector<uint64_t> &src_shortwords, const seq_t *src, const seq_t *snk) {
    // The value 75 is an empirical threshold that works well in practice. O
    // Once a seq A is less than 75% long than another seq B, counting the same number of minhash features does not work.
    if (src->seqlen * 75 > snk->seqlen * 100 || snk->seqlen * 75 > src->seqlen * 100) {
        return comp_n_shared_shortwords(src_shortwords, src, snk); 
    }

#if ENTROHASH
    int ret1 = __builtin_popcountll(src->compressedsign  ^ snk->compressedsign );
    int ret2 = __builtin_popcountll(src->compressedsign2 ^ snk->compressedsign2);
    int ret3 = __builtin_popcountll(src->compressedsign3 ^ snk->compressedsign3);
    int ret4 = __builtin_popcountll(src->compressedsign4 ^ snk->compressedsign4);
    // assert (ret <= NUM_SIGNATURES);
    return MIN(MIN(ret2, ret3), ret4);
    return ret1 + ret2 + ret3 + ret4;
    // return MAX(ret - 32, 0);
#else
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
    LEN_PERC_SRC = LEN_PERC_SNK = 80;
    SIM_BASE = 25;
    if (aa_cnt * 4 > nb_cnt) {
        IS_INPUT_NUC = 0;
        SIM_PERC = 50;
    } else {
        if (baseTcnt < baseUcnt) {
            IS_INPUT_NUC = 1;
            SIM_PERC = 70;
        } else {
            IS_INPUT_NUC = 2;
            LEN_PERC_SRC = LEN_PERC_SNK = 0;
            SIM_PERC = 90;
            SIM_BASE = 0;
        }
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            ALPHA_TYPE_TO_CHAR_TO_REDUCED[j][i] = (char)i;
        }
    }
    std::vector<int> are_args_parsed(argc+1);
    std::fill(are_args_parsed.begin(), are_args_parsed.end(), 0);

    for (int i = 1; i < argc; i += 2) {
        int is_arg_parsed = 1;
        if      (!strcmp("--israndom",       argv[i])) { /* pass */ } 

        else if (!strcmp("--alphasize-seed", argv[i])) { ALPHA_TYPE_TO_SIZE[0] = atoi(argv[i+1]); } 
        else if (!strcmp("--alphasize-sign", argv[i])) { ALPHA_TYPE_TO_SIZE[1] = atoi(argv[i+1]); } 
        else if (!strcmp("--alphasize-sim",  argv[i])) { ALPHA_TYPE_TO_SIZE[2] = atoi(argv[i+1]); } 

        else if (!strcmp("--attempt-ini",    argv[i])) { ATTEMPT_INI           = atoi(argv[i+1]); } 
        else if (!strcmp("--attempt-inc",    argv[i])) { ATTEMPT_INC           = atoi(argv[i+1]); } 
        else if (!strcmp("--attempt-max",    argv[i])) { ATTEMPT_MAX           = atoi(argv[i+1]); } 
       
        else if (!strcmp("--char-per-seed",  argv[i])) { CHAR_PER_SEED         = atoi(argv[i+1]);} 

        else if (!strcmp("--cov-src-max",    argv[i])) { COV_SRC_MAX           = atoi(argv[i+1]); } 
        else if (!strcmp("--cov-snk-max",    argv[i])) { COV_SNK_MAX           = atoi(argv[i+1]); } 
        
        else if (!strcmp("--dbfilt-minseed", argv[i])) { DBFILT_MINSEED        = atoi(argv[i+1]); } 
        else if (!strcmp("--dbfilt-subsamp", argv[i])) { DBFILT_SUBSAMP        = atoi(argv[i+1]); } 
        // else if (!strcmp("--dbfilt-attempt", argv[i])) { DBFILT_ATTEMPT        = atoi(argv[i+1]); } 
        else if (!strcmp("--dbfilt-truehit", argv[i])) { DBFILT_TRUEHIT        = atoi(argv[i+1]); } 

        else if (!strcmp("--is-input-nuc",   argv[i])) { IS_INPUT_NUC          = atoi(argv[i+1]); } 

        else if (!strcmp("--len-perc-src",   argv[i])) { LEN_PERC_SRC          = atoi(argv[i+1]); } 
        else if (!strcmp("--len-perc-snk",   argv[i])) { LEN_PERC_SNK          = atoi(argv[i+1]); } 

        else if (!strcmp("--seed-evalue",    argv[i])) { SEED_EVALUE           = atoi(argv[i+1]); } 
        else if (!strcmp("--seed-length",    argv[i])) { SEED_LENGTH           = atoi(argv[i+1]); } 
        else if (!strcmp("--seed-maxgap",    argv[i])) { SEED_MAXGAP           = atoi(argv[i+1]); }
        else if (!strcmp("--seed-mincnt",    argv[i])) { SEED_MINCNT           = atoi(argv[i+1]); }
        
        else if (!strcmp("--sim-zval",       argv[i])) { SIM_ZVAL              = atoi(argv[i+1]); } 
        else if (!strcmp("--sim-perc",       argv[i])) { SIM_PERC              = atoi(argv[i+1]); } 
        else if (!strcmp("--sim-base",       argv[i])) { SIM_BASE              = atoi(argv[i+1]); } 
        
        else if (!strcmp("--zval-as-sim",    argv[i])) { ZVAL_AS_SIM           = atoi(argv[i+1]); } 
        
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }
    LEN_PERC_SNK = MAX(LEN_PERC_SNK, SIM_PERC);
    
    if (IS_INPUT_NUC) {
        SIGN_LENGTH = (SIM_PERC + 900) / (200 - SIM_PERC); // heuristic values from cd-hit manual page
    } else {
        SIGN_LENGTH = (SIM_PERC + 360) / (150 - SIM_PERC); // heuristic values from cd-hit manual page
    }

    if (2 == IS_INPUT_NUC) {
        SIM_DIFF = MIN(SIM_PERC,  20);
    } else {
        SIM_DIFF = 0;
    }

    for (int i = 1; i < argc; i += 2) {
        int is_arg_parsed = 1;
        if (!strcmp("--sign-length", argv[i]))         { SIGN_LENGTH           = atoi(argv[i+1]); } 
        else if (!strcmp("--sim-diff",       argv[i])) { SIM_DIFF              = atoi(argv[i+1]); }
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }
    
    double p = pow((1 - DBL_EPSILON) * SIM_PERC / 100, SIGN_LENGTH);
    double sd = sqrt(NUM_SIGNATURES * p * (1-p));
    double mean = NUM_SIGNATURES * p;
    SIGN_SHARED_CNT_MIN = MAX(floor(mean - sd * 3), 1); // three standard deviations below the normal distribution of true positives.

    for (int i = 1; i < argc; i += 2) {
        int is_arg_parsed = 1;
        if (!strcmp("--sign-shared-cnt-min", argv[i])) { SIGN_SHARED_CNT_MIN = atoi(argv[i+1]); } 
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
               || !(std::cerr << "Parse command-line parameter " << argv[i] << " at position " << i << " " 
                              << are_args_parsed[i] << " times!" << std::endl));
    }
    
    if (!IS_INPUT_NUC) {
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
    if ((int)seq_ptr->seqlen >= (int)SEED_LENGTH) {
        unsigned int interseed_gap = MAX(1, seq_ptr->seqlen / SEED_MINCNT);
        uint64_t hash = hash_init(seq_ptr->seq);
        seed_add(hash, idx);
        for (int i = SEED_LENGTH; i < (int)seq_ptr->seqlen; i += 1) {
            hash = hash_update(hash, seq_ptr->seq[i-SEED_LENGTH], seq_ptr->seq[i]);
            if (0 == i % interseed_gap) { seed_add(hash, idx); }
        }
    }
}

void seq_signatures_init(seq_t *const seq_ptr) {
#if ENTROHASH  
    seq_ptr->compressedsign = 0;
    seq_ptr->compressedsign2 = 0;
    seq_ptr->compressedsign3 = 0;
    seq_ptr->compressedsign4 = 0;
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
}

void seed_cov(const seed_t *seed, const uint32_t coveringidx, std::set<uint32_t> & visited) {
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

void print_seedsize_histogram(const int hist1[], const int hist2[]) {

    std::cerr << "Seedsize histogram with lower/upper delimiters and frequency before/after filter as columns:" << std::endl;
    for (int i = 0; i < 512; i++) {
        if (hist1[i] > 0 || hist2[i] > 0)
        {
            double lowerlim = (exp(((double)i  ) / 10) - 1 - DBL_EPSILON);
            double upperlim = (exp(((double)i+1) / 10) - 1 - DBL_EPSILON); 
            std::cerr << "\t" << lowerlim << "\t" << upperlim << "\t" << hist1[i] << "\t" << hist2[i] << std::setprecision(15) << std::endl;
        }
    }
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
    DBENTRY_CNT = num_residues / CHAR_PER_SEED + 1;
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

    double INFO_PER_LETTER = exp(ch_entropy); // (IS_INPUT_NUC ? 3.3 : 8.5);
    // 0.33 is used to normalize whole-protein Shannon-information entropy to about 1 bit per position // https://www.ncbi.nlm.nih.gov/pubmed/8804598
    double SHANNON_INFO_PER_LETTER = pow(INFO_PER_LETTER, 0.33); // estimated
    std::cerr << "INFO_PER_LETTER = " << INFO_PER_LETTER << " , SHANNON_INFO_PER_LETTER = " << SHANNON_INFO_PER_LETTER << std::endl;
    std::cerr << "Recommended similarity-threshold for detecting homology = " << 1 / SHANNON_INFO_PER_LETTER << " , actual similarity-threshold = " << SIM_PERC << std::endl;

    if (SEED_EVALUE > 0) {
        double seedlen_fract = log((double)num_residues / SEED_EVALUE + INFO_PER_LETTER) / log(INFO_PER_LETTER);
        int    seedlen_floor = (int)floor(seedlen_fract);
        double seedlen_diff1 = seedlen_fract - seedlen_floor;
        SEED_LENGTH = MIN(MAX(seedlen_floor, 4), 25);
        SEED_MINCNT = (int)floor((100 + 10 - SIM_PERC) / pow(SHANNON_INFO_PER_LETTER, seedlen_diff1) * (IS_INPUT_NUC ? 3 : 1));
        
        std::cerr << "Command-line parameter values after adjustment with SEED_EVALUE = " << SEED_EVALUE << ":" << std::endl;
        std::cerr << "\tSEED_LENGTH = " << SEED_LENGTH << std::endl;
        std::cerr << "\tSEED_MAXGAP = " << SEED_MAXGAP << std::endl;
        std::cerr << "\tSEED_MINCNT = " << SEED_MINCNT << std::endl;
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
    int seedsize_hist2[512];

#if 0
    int maxseedsize = 0;
    for (int i = 0; i < DBENTRY_CNT; i++) {
        maxseedsize = MAX(maxseedsize, seeds[i].size);
    }
    std::cerr << "maxseedsize = " << maxseedsize << std::endl;
#endif

    set_seedsize_histogram(seedsize_hist1);
    
    time(&begtime);
    #pragma omp parallel for schedule(dynamic, 9999*100) 
    for (unsigned int i = 0 ; i < DBENTRY_CNT; i++) {
        if (seeds[i].size > DBFILT_MINSEED) {
            unsigned int randstate = i;
            std::array<std::set<std::pair<uint32_t, uint32_t>>, NUM_SIGNATURES+1> nsigns_to_srcsnk_idx_pairs;
            std::fill(nsigns_to_srcsnk_idx_pairs.begin(), nsigns_to_srcsnk_idx_pairs.end(), std::set<std::pair<uint32_t, uint32_t>>());
            for (int j = 0; j < DBFILT_SUBSAMP; j++) {
                uint32_t srcidx = seeds[i].seqidxs[rand_r(&randstate) % seeds[i].size];
                uint32_t snkidx = seeds[i].seqidxs[rand_r(&randstate) % (seeds[i].size - 1)];
                if (srcidx == snkidx) { snkidx = seeds[i].size-1; }
// Changing seq order affects filter threshold
#if ORDER_AWARE_FILT
                if (seq_arrlist.data[srcidx].seqlen < seq_arrlist.data[snkidx].seqlen) {
                    SWAP(srcidx, snkidx);
                }
#endif
                if (!fail_len_cutoff(&seq_arrlist.data[srcidx], &seq_arrlist.data[snkidx])) {
                    std::vector<uint64_t> src_shortwords;
                    int nsharedsigns = comp_n_shared_signatures(src_shortwords, &seq_arrlist.data[srcidx], &seq_arrlist.data[snkidx]); 
                    nsigns_to_srcsnk_idx_pairs[nsharedsigns].insert(std::make_pair(srcidx, snkidx));
                }
            }
            int nhits = 0;
            int attemptcnt = ATTEMPT_INI;
            for (int nsigns = NUM_SIGNATURES; 
                     nsigns >= SIGN_SHARED_CNT_MIN && attemptcnt > 0 && nhits <= DBFILT_TRUEHIT; 
                     --nsigns) {
                for (auto srcsnk_idx_pair : nsigns_to_srcsnk_idx_pairs[nsigns]) {
                    int srcidx = srcsnk_idx_pair.first;
                    int snkidx = srcsnk_idx_pair.second;
                    double matchprob = 0;
                    double chprobs[256];
                    memset(chprobs, 0, 256 * sizeof(double));
                    for (unsigned int j = 0; j < seq_arrlist.data[srcidx].seqlen; j++) {
                        chprobs[(int)seq_arrlist.data[srcidx].seq[j]] += 1.0 / seq_arrlist.data[srcidx].seqlen;
                    }
                    if (SIM_ZVAL) {
                        for (int j = 0; j < 256; j++) {
                            matchprob += chprobs[j] * chprobs[j];
                        }
                    }
                    // assert(seq_arrlist.data[srcidx].seqlen <= seq_arrlist.data[snkidx].seqlen);
                    int sim = comp_perc_seq_sim_editdist(&seq_arrlist.data[srcidx], &seq_arrlist.data[snkidx], matchprob, chprobs);
                    if (sim >= SIM_PERC) { 
                        attemptcnt = MIN(attemptcnt + ATTEMPT_INC, ATTEMPT_MAX);
                        nhits++;
                        if (!(nhits <= DBFILT_TRUEHIT)) { break; }
                    } else { 
                        attemptcnt--;
                        if (!(attemptcnt > 0)) { break; }
                    }
                }
            }
            if (nhits <= DBFILT_TRUEHIT) {
                seeds[i].size = 0;
                seeds[i].bufsize = 0;
                free(seeds[i].seqidxs);
            }
        }
    }
    time(&endtime);
    
    set_seedsize_histogram(seedsize_hist2);

    print_seedsize_histogram(seedsize_hist1, seedsize_hist2); 
    
    std::cerr << "Filtered " << seedsize_hist2[0] - seedsize_hist1[0] << " seeds in seq index in " << difftime(endtime, begtime) << " seconds." << std::endl;

    printf("%d %d\n", seq_arrlist.size, seq_arrlist.size);
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> coveredarr(BATCHSIZE_INI);
    std::fill(coveredarr.begin(), coveredarr.end(), std::vector<std::pair<uint32_t, uint8_t>>(0));
    unsigned int batchsize = BATCHSIZE_INI;
    time(&begtime);

    for (unsigned int iter = 0; iter < seq_arrlist.size;) {
        unsigned int itermax = MIN(iter+batchsize, seq_arrlist.size);

        #pragma omp parallel for schedule(dynamic, 1)
        for (unsigned int i = iter; i < itermax; i++)
        {
            std::set<uint32_t> visited;
            int filteredcnt = 0;
            int distcompcnt = 0;
            int max_attempts = ATTEMPT_INI;
            int max_attempts_arg = 0;
            if (seq_arrlist.data[i].coveredcnt < COV_SRC_MAX && (int)(seq_arrlist.data[i].seqlen) >= (int)SEED_LENGTH) {
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
                uint64_t hash = hash_init(seq_arrlist.data[i].seq);
                seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                for (unsigned int j = SEED_LENGTH; j < seq_arrlist.data[i].seqlen; j++) {
                    hash = hash_update(hash, seq_arrlist.data[i].seq[j-SEED_LENGTH], seq_arrlist.data[i].seq[j]);
                    seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                }
                std::vector<std::vector<uint32_t>> nsharedsigns_to_coveredidxs_vec(NUM_SIGNATURES + 1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    if ((i != coveredidx)) {
                        int n_shared_signatures = comp_n_shared_signatures(src_shortwords, &seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
                        nsharedsigns_to_coveredidxs_vec.at(n_shared_signatures).push_back(coveredidx);
                    }
                }
                for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_SHARED_CNT_MIN; nsigns--) {
                    filteredcnt += nsharedsigns_to_coveredidxs_vec.at(nsigns).size();
                } 
                int attempts = ATTEMPT_INI;
                for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_SHARED_CNT_MIN && attempts > 0 && attempts > max_attempts - ATTEMPT_MAX; nsigns--) {
                    for (auto coveredidx : nsharedsigns_to_coveredidxs_vec.at(nsigns)) {
                        seq_t *coveringseq = &seq_arrlist.data[i];
                        seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                        if (coveredseq->coveredcnt < COV_SNK_MAX) {
                            int sim = comp_perc_seq_sim_editdist(coveringseq, coveredseq, matchprob, chprobs);

                            if (sim >= SIM_PERC) {
                                coveredarr[i-iter].push_back(std::make_pair(coveredidx, sim));
                                attempts += ATTEMPT_INC;
                                if (attempts > max_attempts) {
                                    max_attempts = attempts;
                                    max_attempts_arg = distcompcnt + 1;
                                }
                            } else {
                                attempts -= 1;
                            }
                            distcompcnt++;
                        }
                        if (!(attempts > 0 && attempts > max_attempts - ATTEMPT_MAX)) { break; }
                    }
                }
                
                std::sort(coveredarr[i-iter].begin(), coveredarr[i-iter].end());
            }
            if (printthresholds.find(i) != printthresholds.end()) {
                time(&endtime);
                fprintf(stderr, "In %.f seconds processed %u sequences.\th1to4=%lu/%i/%i/%lu.\tmax_attempts=%i at %i\tATTEMPT_INI=%i\tATTEMPT_INC=%i\tseqlen=%u\tcoveredcnt=%u\n", 
                        difftime(endtime, begtime), i+1, coveredarr[i-iter].size(), distcompcnt, filteredcnt, visited.size(), max_attempts, max_attempts_arg, 
                        ATTEMPT_INI, ATTEMPT_INC, seq_arrlist.data[i].seqlen, seq_arrlist.data[i].coveredcnt); 
            }
        }
        for (unsigned int i = iter; i < itermax; i++) {
            printf("100 %d", i+1);
            for (auto adj : coveredarr[i-iter]) {
                seq_arrlist.data[adj.first].coveredcnt++;
                assert(adj.second <= 100 || ZVAL_AS_SIM);
                printf(" %d %u", (int)adj.second, adj.first+1);
            }
            printf("\n");
            coveredarr[i-iter].clear();
        }
        assert(coveredarr.size() == batchsize);
        iter += batchsize;
        for (unsigned int i = 0; i < BATCHSIZE_INC; i++) {
            batchsize++;
            coveredarr.push_back(std::vector<std::pair<uint32_t, uint8_t>>(0));
        }
    }
}

