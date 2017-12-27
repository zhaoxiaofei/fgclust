#include "kseq.h"
#include "edlib.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <unordered_map>
#include <set>

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <time.h> 
#include <unistd.h>

#define ENTROHASH 0

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

#define SWAP(a, b) { auto tmp = (a); (a) = (b); (b) = (tmp); }
#define MIN(a, b) (((a) < (b) ? (a) : (b)))
#define MAX(a, b) (((a) > (b) ? (a) : (b)))
// #define SQUARE(a) ((a)*(a))

const auto SQUARE(const auto v) { return v * v; }

#if ENTROHASH
#define NUM_SIGNATURES (256)
#else
#define NUM_SIGNATURES (32)
#endif

KSEQ_INIT(int, read)

// constants
const int PARAM_INIT_SEQCNT = 100;

const int BATCHSIZE_INI = 256; // 15999;
const int BATCHSIZE_INC = 2; // 15999;

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

int COV_SRC_MAX = 5; // 5;
int COV_SNK_MAX = 1000*1000;

int DBFILT_MINSEED = 1000; // 1000*1000; lower -> more filtering, more time saving later
int DBFILT_SUBSAMP = 1000; // 800; // lower -> less filtering accuracy, less time
int DBFILT_ATTEMPT = 50; // 10; // lower -> less filtering accuracy, less time
int DBFILT_TRUEHIT = 5; // lower -> less filtering accuracy

int IS_INPUT_NUC = -1; // guessed

int LEN_PERC = 80;

int SEED_EVALUE = 10;
int SEED_LENGTH = 10; // can be overriden after determination of db size 
int SEED_MAXGAP = 10; // can be overriden after determination of db size
int SEED_MINCNT = 10; // can be overriden after determination of db size

int SIGN_LENGTH = 0;
int SIGN_SHARED_CNT_MIN = 0; 

int SIM_PERC = 0;
int SIM_BASE = 25;
int SIM_ZVAL = 10 * 100; // 11; // real z-score threshold is about 75% of this value due to innacuracy of normal approximation of binomial

bool ZVAL_AS_SIM = false;

// derived vars
uint64_t DBENTRY_CNT = 0; // can be overriden after determination of db size

void showparams() {
    
    std::cerr << " BATCHSIZE_INI = " << BATCHSIZE_INI << std::endl;
    std::cerr << " NUM_SIGNATURES = " << NUM_SIGNATURES << std::endl;

    std::cerr << " ALPHA_TYPE_TO_SIZE[0] = ALPHASIZE_SEED = " << ALPHA_TYPE_TO_SIZE[0] << std::endl;
    std::cerr << " ALPHA_TYPE_TO_SIZE[1] = ALPHASIZE_SIGN = " << ALPHA_TYPE_TO_SIZE[1] << std::endl;
    std::cerr << " ALPHA_TYPE_TO_SIZE[2] = ALPHASIZE_SIM  = " << ALPHA_TYPE_TO_SIZE[2] << std::endl;

    std::cerr << " ATTEMPT_INI = " << ATTEMPT_INI << std::endl;
    std::cerr << " ATTEMPT_INC = " << ATTEMPT_INC << std::endl;
    std::cerr << " ATTEMPT_MAX = " << ATTEMPT_MAX << std::endl;
    
    std::cerr << " CHAR_PER_SEED = "  << CHAR_PER_SEED << std::endl;
    
    std::cerr << " COV_SRC_MAX = " << COV_SRC_MAX  << std::endl;
    std::cerr << " COV_SNK_MAX = " << COV_SNK_MAX << std::endl;
    
    std::cerr << " DBFILT_MINSEED = " << DBFILT_MINSEED << std::endl;
    std::cerr << " DBFILT_SUBSAMP = " << DBFILT_SUBSAMP << std::endl; 
    std::cerr << " DBFILT_TRUEHIT = " << DBFILT_TRUEHIT << std::endl;
    std::cerr << " DBFILT_ATTEMPT = " << DBFILT_ATTEMPT << std::endl;

    std::cerr << " IS_INPUT_NUC = " << IS_INPUT_NUC         << std::endl;
    
    std::cerr << " LEN_PERC = " << (int)LEN_PERC << std::endl;
        
    std::cerr << " SEED_EVALUE = " << SEED_EVALUE << std::endl;
    std::cerr << " SEED_LENGTH = " << SEED_LENGTH << std::endl;
    std::cerr << " SEED_MAXGAP = " << SEED_MAXGAP << std::endl;
    std::cerr << " SEED_MINCNT = " << SEED_MINCNT << std::endl;

    std::cerr << " SIGN_LENGTH = " << SIGN_LENGTH    << std::endl;
    std::cerr << " SIGN_SHARED_CNT_MIN = " << SIGN_SHARED_CNT_MIN      << std::endl;
    
    std::cerr << " SIM_PERC = " << (int)SIM_PERC << std::endl;
    std::cerr << " SIM_BASE = " << (int)SIM_BASE << std::endl;
    std::cerr << " SIM_ZVAL = " << (int)SIM_ZVAL << std::endl;

    std::cerr << " ZVAL_AS_SIM = " << ZVAL_AS_SIM << std::endl;    
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
    std::cerr << "  --dbfilt-attempt\t: number of attempts on most likely (minhash) similar SP for seed pruning.[" << DBFILT_ATTEMPT << "]" << std::endl;
    std::cerr << "  --dbfilt-truepos\t: maximum number of true posive SP hits to trigger pruning.["                << DBFILT_TRUEHIT << "]" << std::endl;
    
    std::cerr << "  --is-input-nuc\t: 0, 1, and 2 mean input is protein, RNA, and DNA, respectively (auto detected). [" << IS_INPUT_NUC << "]" << std::endl;

    std::cerr << "  --len-perc\t: A covers B only if min(len(A) / len(B), len(B) / len(A)) >= len_perc. [" << LEN_PERC << "]" << std::endl; 
    
    std::cerr << "  --seed-length\t: length of an indexed seed. overwrite --seed-evalue. ["     << SEED_LENGTH << "]" << std::endl;
    std::cerr << "  --seed-maxgap\t: max number of residues between consecutive seeds. ["       << SEED_MAXGAP << "]" << std::endl;
    std::cerr << "  --seed-mincnt\t: minimum number of seeds per sequence.["                    << SEED_MINCNT << "]" << std::endl;
    std::cerr << "  --seed-evalue\t: evalue for seed hit (higher leads to slower runtime). ["   << SEED_EVALUE << "]" << std::endl;

    std::cerr << "  --sign-length\t: length of k-mers for computing minhash values. ["                      << SIGN_LENGTH         << "]" << std::endl;
    std::cerr << "  --sign-shared-cnt-min\t: minimum number of minhash values to trigger sequence search [" << SIGN_SHARED_CNT_MIN << "]" << std::endl; 

    std::cerr << "  --sim-zval\t: percent of standard deviations above similarity by chance for a sequence to cover another ["  << SIM_ZVAL << "]" << std::endl;
    std::cerr << "            \t: the effective value is 0.2 to 0.4 times less than this set value due to normal approximation of skewed binomial" << std::endl;
    std::cerr << "  --sim-base\t: A covers B only if (len(B) - edit-distance(A, B)) >= sim_perc * len(B) + sim_base [sim_base=" << SIM_BASE << "]" << std::endl;
    std::cerr << "  --sim-perc\t:   where gaps at ends of B are not penalized. [sim_perc="                                      << SIM_PERC << "]" << std::endl;
               
    std::cerr << "  --zval-as-sim\t: used the sim-zval as similarity threshold [" << ZVAL_AS_SIM   << "]" << std::endl;

    std::cerr << "Note: default value of 0 means dependence to other parameters or to the input." << std::endl;
    exit(-1);
}

void alphareduce(const char *const strarg, const int reducetype) {
    const char *str = strarg;
    for (; *str; str++) {
        ALPHA_TYPE_TO_CHAR_TO_REDUCED[reducetype][*str] = *strarg;
    }
}

int calc_vecnorm(int a, int b) {
    return (int)ceil(sqrt(a * a + b * b));
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
        ret += (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][beg[i]];
        ret %= PRIME_MOD; 
    }
    return ret;
}

const uint64_t hash_update(uint64_t hash, char prv, char nxt) {
    hash = hash * PRIME_BASE + (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][nxt];
    hash += PRIME_MOD;
    hash -= (SEED_POWER * (uint64_t)ALPHA_TYPE_TO_CHAR_TO_REDUCED[0][prv]) % PRIME_MOD;
    return hash % PRIME_MOD;
}

const uint64_t sign_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < SIGN_LENGTH; i++) {
        ret *= SIGN_BASE;
        ret += (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][beg[i]];
        ret %= SIGN_MOD; 
    }
    return ret;
}

const uint64_t sign_update(uint64_t hash, char prv, char nxt) {
    hash = hash * SIGN_BASE + (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][nxt];
    hash += SIGN_MOD;
    hash -= (SIGN_POWER * (uint64_t) ALPHA_TYPE_TO_CHAR_TO_REDUCED[1][prv]) % SIGN_MOD;
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

const bool fail_len_cutoff(const seq_t *seq1, const seq_t *seq2) {
    return seq1->seqlen * LEN_PERC > seq2->seqlen * 100 || seq2->seqlen * LEN_PERC > seq1->seqlen * 100;
}

static const int calc_perc_seq_sim_editdist(const seq_t *seq1, const seq_t *seq2, const double matchprob) {
    
    int maxEditDist = seq1->seqlen - ceil((1 - DBL_EPSILON) * (sqrt(SQUARE((double)SIM_BASE) + SQUARE((double)(seq1->seqlen * SIM_PERC) / 100.0))));
    if (maxEditDist < 0) { return 0; }
    
    EdlibAlignResult result;
    result = edlibAlign(seq1->seq, seq1->seqlen, seq2->seq, seq2->seqlen,
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
    int editdist = result.editDistance;
    edlibFreeAlignResult(result);

    assert(editdist >= -1);
    assert(editdist <= (int)seq1->seqlen || fprintf(stderr, ">%s\n%s\n>%s\n%s\neditdist=%d>seq1len\n", seq1->name, seq1->seq, seq2->name, seq2->seq, editdist));
    // assert(editdist <= (int)seq2->seqlen || fprintf(stderr, ">%s\n%s\n>%s\n%s\neditdist=%d>seq2len\n", seq1->name, seq1>seq, seq2->name, seq2->seq, editdist));
    
    if (-1 == editdist) { return 0; }
    if (SIM_ZVAL) {
        
        double mean = seq1->seqlen * matchprob;
        double var  = seq1->seqlen * matchprob * (1 - matchprob);

        double zval = ((double)(seq1->seqlen - editdist) - mean) /  sqrt(var);
        //fprintf(stderr, "seq1=%s,seq2=%s,mean=%f,var=%f,zval=%f,seqlen1=%d,editdist=%d\n", seq1->name, seq2->name, mean, var, zval, seq1->seqlen, editdist);
        
        if (zval < SIM_ZVAL / 100.0) { return 0; }
        if (ZVAL_AS_SIM) { return MIN(200, 100 + floor(zval)); }
    }
    int ret = 100 * (seq1->seqlen - editdist) / seq1->seqlen;
    return ret;
}

std::vector<uint64_t> compsigns(const seq_t *seq_ptr) {
    std::vector<uint64_t> signs;
    if ((int)seq_ptr->seqlen >= (int)SIGN_LENGTH) {
        signs.reserve((int)seq_ptr->seqlen - (int)SIGN_LENGTH + 1);
        uint64_t sign = sign_init(seq_ptr->seq);
        signs.push_back(sign);
        for (int i = SIGN_LENGTH; i < seq_ptr->seqlen; i += 1) {
            sign = sign_update(sign, seq_ptr->seq[i-SIGN_LENGTH], seq_ptr->seq[i]);
            signs.push_back(sign);
        }
        std::sort(signs.rbegin(), signs.rend());
    }
    return signs;
}

int __attribute__((noinline)) calc_n_shared_kmers(const seq_t *seq1, const seq_t *seq2) {
    std::vector<uint64_t> signs1 = compsigns(seq1);
    std::vector<uint64_t> signs2 = compsigns(seq2);
    auto it1 = signs1.begin();
    auto it2 = signs2.begin();
    int ret = 0;
    while (it1 != signs2.end() && it2 != signs2.end()) {
        if ((*it1) == (*it2)) {
            ret++;
            it1++;
        } else if ((*it1) < (*it2)) {
            it2++;
        } else if ((*it1) > (*it2)) {
            it1++;
        }
    }
    int ret2 = ret * 32 / signs2.size();
    assert (ret2 <= 32);
    return ret2;
}

int calc_n_shared_signatures(const seq_t *seq1, const seq_t *seq2) {
    // The value 75 is an empirical threshold that works well in practice. O
    // Once a seq A is less than 75% long than another seq B, counting the same number of minhash features does not work.
    if (seq1->seqlen * 75 > seq2->seqlen * 100 || seq2->seqlen * 75 > seq1->seqlen * 100) {
        return calc_n_shared_kmers(seq1, seq2); 
    }

#if ENTROHASH
    int ret1 = __builtin_popcountll(seq1->compressedsign  ^ seq2->compressedsign );
    int ret2 = __builtin_popcountll(seq1->compressedsign2 ^ seq2->compressedsign2);
    int ret3 = __builtin_popcountll(seq1->compressedsign3 ^ seq2->compressedsign3);
    int ret4 = __builtin_popcountll(seq1->compressedsign4 ^ seq2->compressedsign4);
    // assert (ret <= NUM_SIGNATURES);
    return MIN(MIN(ret2, ret3), ret4);
    return ret1 + ret2 + ret3 + ret4;
    // return MAX(ret - 32, 0);
#else
    int i = 0;
    int j = 0;
    int ret = 0;
    while (i != NUM_SIGNATURES && j != NUM_SIGNATURES) {
        if (seq1->signatures[i] == seq2->signatures[j]) {
            ret++; 
            i++; 
            j++;
        } else if (seq1->signatures[i] < seq2->signatures[j]) {
            j++;
        } else if (seq1->signatures[i] > seq2->signatures[j]) {
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
    for (int i = 0; i < MIN(seq_arrlist.size, PARAM_INIT_SEQCNT); i++) {
        for (int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
            if (NULL != strchr("ACGTUacgtu", seq_arrlist.data[i].seq[j])) {
                nb_cnt++;
            } else {
                aa_cnt++;
            }
        }
    }
    if (aa_cnt * 4 > nb_cnt) {
        IS_INPUT_NUC = 0;
        SIM_PERC = 50;
    } else {
        IS_INPUT_NUC = 1;
        SIM_PERC = 90;
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            ALPHA_TYPE_TO_CHAR_TO_REDUCED[j][i] = (char)i;
        }
    }
    std::vector<int> are_args_parsed(argc+1);
    std:fill(are_args_parsed.begin(), are_args_parsed.end(), 0);

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
        else if (!strcmp("--dbfilt-attempt", argv[i])) { DBFILT_ATTEMPT        = atoi(argv[i+1]); } 
        else if (!strcmp("--dbfilt-truehit", argv[i])) { DBFILT_TRUEHIT        = atoi(argv[i+1]); } 

        else if (!strcmp("--is-input-nuc",   argv[i])) { IS_INPUT_NUC          = atoi(argv[i+1]); } 

        else if (!strcmp("--len-perc",       argv[i])) { LEN_PERC              = atoi(argv[i+1]); } 
                
        else if (!strcmp("--seed-evalue",    argv[i])) { SEED_EVALUE           = atoi(argv[i+1]); } 
        else if (!strcmp("--seed-length",    argv[i])) { SEED_LENGTH           = atoi(argv[i+1]); } 
        else if (!strcmp("--seed-maxgap",    argv[i])) { SEED_MAXGAP           = atoi(argv[i+1]); }
        else if (!strcmp("--seed-mincnt",    argv[i])) { SEED_MINCNT           = atoi(argv[i+1]); }
        
        else if (!strcmp("--sim-zval",       argv[i])) { SIM_ZVAL              = atoi(argv[i+1]);} 
        else if (!strcmp("--sim-perc",       argv[i])) { SIM_PERC              = atoi(argv[i+1]);} 
        else if (!strcmp("--sim-base",       argv[i])) { SIM_BASE              = atoi(argv[i+1]); } 

        else if (!strcmp("--zval-as-sim",    argv[i])) { ZVAL_AS_SIM           = atoi(argv[i+1]);} 
        
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }

    if (IS_INPUT_NUC) {
        SIGN_LENGTH = (SIM_PERC + 900) / (200 - SIM_PERC); // heuristic values from cd-hit manual page
    } else {
        SIGN_LENGTH = (SIM_PERC + 360) / (150 - SIM_PERC); // heuristic values from cd-hit manual page
    }
    
    for (int i = 1; i < argc; i += 2) {
        int is_arg_parsed = 1;
        if (!strcmp("--sign-length", argv[i])) { SIGN_LENGTH = atoi(argv[i+1]); } 
        else { is_arg_parsed = 0; }
        are_args_parsed[i]   += is_arg_parsed;
        are_args_parsed[i+1] += is_arg_parsed;
    }
    
    double p = pow(SIM_PERC / 100.001, SIGN_LENGTH);
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
        int kmerspace = MAX(1, MIN(SEED_MAXGAP, seq_ptr->seqlen / SEED_MINCNT));
        uint64_t hash = hash_init(seq_ptr->seq);
        seed_add(hash, idx);
        for (int i = SEED_LENGTH; i < (int)seq_ptr->seqlen; i += 1) {
            hash = hash_update(hash, seq_ptr->seq[i-SEED_LENGTH], seq_ptr->seq[i]);
            if (0 == i % kmerspace) { seed_add(hash, idx); }
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
        for (int i = SIGN_LENGTH; i < seq_ptr->seqlen; i += 1) {
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
    for (int i = 0; i < seed->size; i++) { 
        int coveredidx = seed->seqidxs[i];
        if (!fail_len_cutoff(coveringseq, &seq_arrlist.data[coveredidx])) {
            visited.insert(coveredidx);
        }
    }
}

void print_seedsize_histogram(int seedsize_histogram[], const char *name) {
    memset(seedsize_histogram, 0, (1000+1) * sizeof(int));
    for (int i = 0 ; i < DBENTRY_CNT; i++) {
        uint32_t seedsize = (seeds[i].size > 1000 ? 1000 : seeds[i].size);
        seedsize_histogram[seedsize]++;
    }
    std::cerr << "Start of " << name << std::endl;
    for (int i = 0; i < 1000+1; i++) {
        std::cerr << i << "\t";
    }
    std::cerr << std::endl;
    for (int i = 0; i < 1000+1; i++) {
        std::cerr << seedsize_histogram[i] << "\t";
    }
    std::cerr << std::endl;
    std::cerr << "End of " << name << std::endl;
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
    std::cerr << "For usage and help, please enter either one of the following commands:" << std::endl;
    std::cerr << "  " << argv[0] << " --help" << std::endl;
    std::cerr << "  " << argv[0] << " --usage" << std::endl;
    std::cerr << "Reading fasta sequences from STDIN" << std::endl; 

    std::set<int> printthresholds;
    for (int i = 0; i < 300; i++) {
        double thres = (double)((i+1)*(BATCHSIZE_INI+1)) * pow(1.05, (double)i);
        if (thres * 1.01 < (double)(INT_MAX)) { printthresholds.insert((int)thres); }
    }

    seq_arrlist_init();
    
    kseq_t *kseq = kseq_init(fileno(stdin));
    int i = 0;
    time(&begtime);
    while ( kseq_read(kseq) >= 0 ) {
        seq_arrlist_add(kseq);
        if (printthresholds.find(++i) != printthresholds.end()) {
            time(&endtime);
            fprintf(stderr, "Read %d sequences in %.f seconds.\n", i, difftime(endtime, begtime));
        }
        if (PARAM_INIT_SEQCNT == i) {
            PARAMS_init(argc, argv);
        }
    }
    kseq_destroy(kseq);
    
    if (i < PARAM_INIT_SEQCNT) { PARAMS_init(argc, argv); }
    
    // reinitialize some vars
    uint64_t num_residues = 0;
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        num_residues += seq_arrlist.data[i].seqlen;
    }
    DBENTRY_CNT = num_residues / CHAR_PER_SEED + 1;
    std::cerr << "Derived paramters: " << std::endl;
    std::cerr << "  NUM_SEQS = " << seq_arrlist.size << std::endl;
    std::cerr << "  NUM_RESIDUES = " << num_residues << std::endl;
    std::cerr << "  DBENTRY_CNT = " << DBENTRY_CNT << std::endl;
    
    uint64_t ch_to_cnt[256];
    memset(ch_to_cnt, 0, 256 * sizeof(uint64_t));
    // subsample 1000 sequences in the database and 1000 letters in each sequence
    for (int i = 0 ; i < seq_arrlist.size; i += MAX(1, seq_arrlist.size/1000)) {
        for (int j = 0; j < seq_arrlist.data[i].seqlen; j += MAX(1, seq_arrlist.data[i].seqlen/1000)) {
            char c = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][seq_arrlist.data[i].seq[j]];
            ch_to_cnt[c]++;
        }
    }
    uint64_t ch_totcnt = 0;
    for (int i = 0; i < 256; i++) {
        ch_totcnt += ch_to_cnt[i];
    }
    double ch_to_prob[256];
    double ch_entropy = 0;
    std::cerr << "Letter frequency distribution:" << std::endl;
    for (int i = 0; i < 256; i++) {
        ch_to_prob[i] = ((double)ch_to_cnt[i] + 1.0 / 256) / ((double)ch_totcnt + 1);
        ch_entropy -= ch_to_prob[i] * log(ch_to_prob[i]);
        std::cerr << i << " : " << ch_to_cnt[i] << " , "; // print all ascii characters
    }
    std::cerr << std::endl;
    for (int i = 32; i < 127; i++) {
        std::cerr << (char)i << " : " << ch_to_cnt[i] << " , "; // print printable ascii characters
    }
    std::cerr << std::endl;

    double INFO_PER_LETTER = exp(ch_entropy); // (IS_INPUT_NUC ? 3.3 : 8.5);
    std::cerr << "INFO_PER_LETTER = " << INFO_PER_LETTER << std::endl;
    if (SEED_EVALUE > 0) {
        double adjusted_kmer_size = log((double)num_residues / SEED_EVALUE + INFO_PER_LETTER) / log(INFO_PER_LETTER);
        int inf_kmer_size = (int)floor(adjusted_kmer_size);
        int sim_kmer_size = 150 / (100 - MIN((int)SIM_PERC, 99));
        if (inf_kmer_size > sim_kmer_size) { 
            int adjusted_kmer_space_ratio = (int)(100.0 * (adjusted_kmer_size - floor(adjusted_kmer_size)));
            SEED_LENGTH = inf_kmer_size;
            SEED_MAXGAP = (SIM_PERC / 5) * (100 + adjusted_kmer_space_ratio) / 100; 
            SEED_MINCNT = (120 - SIM_PERC) * 100 / (100 + adjusted_kmer_space_ratio);
        } else {
            SEED_LENGTH = sim_kmer_size;
            SEED_MAXGAP = (SIM_PERC / 5);
            SEED_MINCNT = (120 - SIM_PERC); 
        }
        SEED_LENGTH = MIN(MAX(SEED_LENGTH, 7), 25);
        SEED_MAXGAP = MAX(SEED_MAXGAP, 1);
        SEED_MINCNT = MAX(SEED_MINCNT, 20);
        
        std::cerr << "Command-line parameter values after adjustment with SEED_EVALUE = " << SEED_EVALUE << ":" << std::endl;
        std::cerr << "  SEED_LENGTH = " << SEED_LENGTH << std::endl;
        std::cerr << "  SEED_MAXGAP = " << SEED_MAXGAP << std::endl;
        std::cerr << "  SEED_MINCNT = " << SEED_MINCNT << std::endl;
    } else {
        std::cerr << "No adjustment made because SEED_EVALUE = " << SEED_EVALUE << std::endl;
    }

    hash_sign_INIT();

    seeds = (seed_t*) xmalloc(DBENTRY_CNT * sizeof(seed_t));
    memset(seeds, 0, DBENTRY_CNT * sizeof(seed_t));
    
    time(&begtime);
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        seq_longword_init(&seq_arrlist.data[i], i);
        if (printthresholds.find(i) != printthresholds.end()) {
            time(&endtime);
            fprintf(stderr, "Indexed %d sequences in %.f seconds.\n", i, difftime(endtime, begtime));
        }
    }

    #pragma omp parallel for schedule(dynamic, 9999*10)
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        seq_signatures_init(&seq_arrlist.data[i]);
        for (int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
            seq_arrlist.data[i].seq[j] = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][seq_arrlist.data[i].seq[j]];
        }
    }
    
    int seedsize_histogram[1000+1];
    print_seedsize_histogram(seedsize_histogram, "seedsize_histogram 1"); 
    
    int maxseedsize = 0;
    for (int i = 0; i < DBENTRY_CNT; i++) {
        maxseedsize = MAX(maxseedsize, seeds[i].size);
    }
    
    std::cerr << "maxseedsize = " << maxseedsize << std::endl;

    #pragma omp parallel for schedule(dynamic, 9999*100) 
    for (int i = 0 ; i < DBENTRY_CNT; i++) {
        if (seeds[i].size > DBFILT_MINSEED) {
            unsigned int randstate = 1 + (unsigned int)i;
            int probsimcnt = 0;
            std::array<std::set<std::pair<uint32_t, uint32_t>>, NUM_SIGNATURES+1> nsigns_to_seqidx_pairs;
            std::fill(nsigns_to_seqidx_pairs.begin(), nsigns_to_seqidx_pairs.end(), std::set<std::pair<uint32_t, uint32_t>>());
            for (int j = 0; j < DBFILT_SUBSAMP; j++) {
                uint32_t seqidx1 = seeds[i].seqidxs[rand_r(&randstate) % seeds[i].size];
                uint32_t seqidx2 = seeds[i].seqidxs[rand_r(&randstate) % (seeds[i].size - 1)];
                if (seqidx1 == seqidx2) { seqidx2 = seeds[i].size-1; }
                if (seq_arrlist.data[seqidx1].seqlen > seq_arrlist.data[seqidx2].seqlen) {
                    SWAP(seqidx1, seqidx2);
                }
                if (!fail_len_cutoff(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2])) {
                    int nsharedsigns = calc_n_shared_signatures(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]); 
                    nsigns_to_seqidx_pairs[nsharedsigns].insert(std::make_pair(seqidx1, seqidx2));
                }
            }
            int nhits = 0;
            int attemptcnt = 0;
            for (int nsigns = NUM_SIGNATURES; 
                    nsigns >= SIGN_SHARED_CNT_MIN && attemptcnt < DBFILT_ATTEMPT && nhits <= DBFILT_TRUEHIT; 
                    nsigns--, attemptcnt++) {
                for (auto seqidxpair : nsigns_to_seqidx_pairs[nsigns]) {
                    int seqidx1 = seqidxpair.first;
                    int seqidx2 = seqidxpair.second;
                    double chprobs[256];
                    memset(chprobs, 0, 256 * sizeof(double));
                    for (int j = 0; j < seq_arrlist.data[seqidx2].seqlen; j++) {
                        chprobs[seq_arrlist.data[seqidx2].seq[j]] += 1.0 / seq_arrlist.data[seqidx2].seqlen;
                    }
                    double matchprob = 0;
                    for (int j = 0; j < 256; j++) {
                        matchprob += chprobs[j] * chprobs[j];
                    }
                    int sim = calc_perc_seq_sim_editdist(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2], matchprob);
                    if (sim >= SIM_PERC) { nhits++; }
                }
            }
            if (nhits <= DBFILT_TRUEHIT) {
                seeds[i].size = 0;
                seeds[i].bufsize = 0;
                free(seeds[i].seqidxs);
            }
        }
    }
    print_seedsize_histogram(seedsize_histogram, "seedsize_histogram 2"); 
    
    printf("%d %d\n", seq_arrlist.size, seq_arrlist.size);
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> coveredarr(BATCHSIZE_INI);
    std::fill(coveredarr.begin(), coveredarr.end(), std::vector<std::pair<uint32_t, uint8_t>>(0));
    int batchsize = BATCHSIZE_INI;
    time(&begtime);

    for (int64_t iter = 0; iter < seq_arrlist.size;) {
        int itermax = MIN(iter+batchsize, seq_arrlist.size);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = iter; i < itermax; i++)
        {
            std::set<uint32_t> visited;
            int filteredcnt = 0;
            int distcompcnt = 0;
            int max_attempts = ATTEMPT_INI;
            int max_attempts_arg = 0;
            if (seq_arrlist.data[i].coveredcnt <= COV_SRC_MAX && (int)(seq_arrlist.data[i].seqlen) >= (int)SEED_LENGTH) {
                double chprobs[256];
                memset(chprobs, 0, 256 * sizeof(double));
                for (int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
                    chprobs[seq_arrlist.data[i].seq[j]] += 1.0 / seq_arrlist.data[i].seqlen;
                }
                double matchprob = 0;
                for (int j = 0; j < 256; j++) {
                    matchprob += chprobs[j] * chprobs[j];
                }

                uint64_t hash = hash_init(seq_arrlist.data[i].seq);
                seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                for (int j = SEED_LENGTH; j < seq_arrlist.data[i].seqlen; j++) {
                    hash = hash_update(hash, seq_arrlist.data[i].seq[j-SEED_LENGTH], seq_arrlist.data[i].seq[j]);
                    seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                }
                std::vector<std::vector<uint32_t>> nsharedsigns_to_coveredidxs_vec(NUM_SIGNATURES + 1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                    // fprintf(stderr, "%s has link to %s\n", seq_arrlist.data[i].name, coveredseq->name);
                    if ((i != coveredidx)) {
                        int n_shared_signatures = calc_n_shared_signatures(&seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
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
                        if (coveredseq->coveredcnt <= COV_SNK_MAX) {
                            int sim = calc_perc_seq_sim_editdist(coveredseq, coveringseq, matchprob);

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
                fprintf(stderr, "In %.f seconds processed %d sequences.\th1to4=%u/%i/%i/%u.\tmax_attempts=%i at %i\tATTEMPT_INI=%i\tATTEMPT_INC=%i\tseqlen=%u\tcoveredcnt=%u\n", 
                        difftime(endtime, begtime), i+1, coveredarr[i-iter].size(), distcompcnt, filteredcnt, visited.size(), max_attempts, max_attempts_arg, 
                        ATTEMPT_INI, ATTEMPT_INC, seq_arrlist.data[i].seqlen, seq_arrlist.data[i].coveredcnt); 
            }
        }
        for (int i = iter; i < itermax; i++) {
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
        for (int i = 0; i < BATCHSIZE_INC; i++) {
            batchsize++;
            coveredarr.push_back(std::vector<std::pair<uint32_t, uint8_t>>(0));
        }
    }
}

