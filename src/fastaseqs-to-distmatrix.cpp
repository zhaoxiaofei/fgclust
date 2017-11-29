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
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <unistd.h>

void *xmalloc(size_t size) {
    void *ret = malloc(size);
    if (NULL == ret) {
        fprintf(stderr, "malloc failed!\n");
        abort();
    }
    return ret;
}

void *xcalloc(size_t num, size_t size) {
    void *ret = calloc(num, size);
    if (NULL == ret) {
        fprintf(stderr, "calloc failed!\n");
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
#define UPDATE_MIN(a, b) ((a = MIN(a, b)))
#define UPDATE_MAX(a, b) ((a = MAX(a, b)))
#define SQUARE(a) ((a)*(a))

KSEQ_INIT(int, read)

// constants

const uint64_t PRIME_BASE = 48271L;
const uint64_t SIGN_BASE  = 48271L; 
const uint64_t PRIME_MOD  = (0x1L << 31L) - 1L; 
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L; 

// variables that are initialized once from constants and other variables

uint64_t SEED_POWER = 0; 
uint64_t SIGN_POWER = 0; 
char ALPHA_TYPE_TO_CHAR_TO_REDUCED[3][256]; // derived from ALPHA_TYPE_TO_SIZE

// variables that are initialized from command line args

int LEN_PERC = 0;

uint64_t BATCH_SIZE = 256; // 15999;
uint64_t RESIDUE_CNT_TO_SEED_CNT_RATIO = 40;

int IS_INPUT_NUC = 0; // guessed
int ALPHA_TYPE_TO_SIZE[] = {10, 10, 10};
int SIM_MODE = 0;
int SIM_PERC = 0;
int SIM_BASE = 0;

int DBENTRY_FILT_OCC_MIN = 1000;
int DBENTRY_FILT_PAIR_SUBSAMP_CNT = 800;
int DBENTRY_FILT_PAIR_ATTEMPT_CNT = 10;
int DBENTRY_FILT_PAIR_TRUEHIT_MAX = 4;

// dbsize-related variables

uint64_t DBENTRY_CNT = 0; // can be overriden after determination of db size

// variables that are derived from SIM_PERC but can be overriden in command line

int COV_SRC_MAX = 2;
int COV_SNK_MAX = 16;

int SEED_LENGTH = 0; // can be overriden after determination of db size 
int SEED_MAXGAP = 0; // can be overriden after determination of db size
int SEED_MINCNT = 0; // can be overriden after determination of db size

int SIGN_LENGTH = 0;
int SIGN_SHARED_PERC_MIN = 0; 
int SIGN_SHARED_PERC_ZSCORE = 300;
int HASH_MIN_RATIO = 8;

int ATTEMPT_INI = 50;
int ATTEMPT_INC = 75;
int ATTEMPT_MAX = 100;
int ATTEMPT_MIN = -100;

void showparams() {
    
    std::cerr << " BATCH_SIZE = " << BATCH_SIZE << std::endl;

    std::cerr << " RESIDUE_CNT_TO_SEED_CNT_RATIO = " << RESIDUE_CNT_TO_SEED_CNT_RATIO << std::endl;
    std::cerr << " HASH_MIN_RATIO = " << HASH_MIN_RATIO << std::endl;

    std::cerr << " IS_INPUT_NUC = " << IS_INPUT_NUC         << std::endl;
    std::cerr << " SIM_MODE = " << (int)SIM_MODE << std::endl;
    std::cerr << " SIM_PERC = " << (int)SIM_PERC << std::endl;
    std::cerr << " SIM_BASE = " << (int)SIM_BASE << std::endl;
    std::cerr << " LEN_PERC = " << (int)LEN_PERC << std::endl;

    std::cerr << " SEED_LENGTH = " << SEED_LENGTH    << std::endl;
    std::cerr << " SEED_MAXGAP = " << SEED_MAXGAP    << std::endl;
    std::cerr << " SIGN_LENGTH = " << SIGN_LENGTH    << std::endl;
    std::cerr << " SIGN_SHARED_PERC_MIN = " << SIGN_SHARED_PERC_MIN      << std::endl;
    std::cerr << " SIGN_SHARED_PERC_ZSCORE = " << SIGN_SHARED_PERC_ZSCORE      << std::endl;
    std::cerr << " COV_SRC_MAX = " << COV_SRC_MAX  << std::endl;
    std::cerr << " COV_SNK_MAX = " << COV_SNK_MAX << std::endl;

    std::cerr << " ATTEMPT_INI = " << ATTEMPT_INI << std::endl;
    std::cerr << " ATTEMPT_INC = " << ATTEMPT_INC << std::endl;
    std::cerr << " ATTEMPT_MAX = " << ATTEMPT_MAX << std::endl;
    std::cerr << " ATTEMPT_MIN = " << ATTEMPT_MIN << std::endl;

    std::cerr << " DBENTRY_FILT_OCC_MIN = " << DBENTRY_FILT_OCC_MIN << std::endl;
    std::cerr << " DBENTRY_FILT_PAIR_SUBSAMP_CNT = " << DBENTRY_FILT_PAIR_SUBSAMP_CNT << std::endl;
    
    std::cerr << " DBENTRY_FILT_PAIR_TRUEHIT_MAX = " << DBENTRY_FILT_PAIR_TRUEHIT_MAX << std::endl;
    std::cerr << " DBENTRY_FILT_PAIR_ATTEMPT_CNT = " << DBENTRY_FILT_PAIR_ATTEMPT_CNT << std::endl;
}

void show_usage(const int argc, const char *const *const argv) {
    std::cerr << "Program : " << argv[0] << std::endl;
    std::cerr << "Command-line arguments with [default-values]:" << std::endl;
    std::cerr << "  --sim-mode\t: type of similarity ["                                                    << SIM_MODE << "]" << std::endl;
    std::cerr << "            \t:   1 : use (len(A) - edit-distance(A, B)) as sim" << std::endl;
    std::cerr << "            \t:   2 : use (alignment-len - edit-distance(A, B)) as sim" << std::endl;
    std::cerr << "  --sim-base\t: A covers B only if sim^2 >= (sim_perc*len(B))^2 + sim-base^2 [sim-base=" << SIM_BASE << "]" << std::endl;
    std::cerr << "  --sim-perc\t:   where gaps at ends of B are not penalized. [sim-perc="                 << SIM_PERC << "]" << std::endl;
    std::cerr << "  --len-perc\t: A covers B only if min(len(A) / len(B), len(B) / len(A)) >= len_perc. [" << LEN_PERC << "]" << std::endl;
    std::cerr << "  --cov-snk-max\t: max number of times that the covered sequence can be covered.["       << COV_SNK_MAX << "]" << std::endl;
    std::cerr << "  --cov-src-max\t: max number of times that the covering sequence can be coverered. ["   << COV_SRC_MAX << "]" << std::endl;
    std::cerr << "  --attempt-ini\t: initial number of attempts. ["                            << ATTEMPT_INI << "]" << std::endl;
    std::cerr << "  --attempt-inc\t: number of attempts incremented per true positive hits. [" << ATTEMPT_INC << "]" << std::endl;
    std::cerr << "  --attempt-max\t: number of attempts capped at this maximum value. ["       << ATTEMPT_MAX << "]" << std::endl;
    std::cerr << "  --attempt-min\t: number of attempts capped at this minimum value. ["       << ATTEMPT_MIN << "]" << std::endl;
    std::cerr << "  --sign-length\t: length of k-mers for computing minhash values. ["         << SIGN_LENGTH << "]" << std::endl;
    std::cerr << "  --seed-length\t: length of an indexed seed. ["                             << SEED_LENGTH << "]" << std::endl;
    std::cerr << "  --seed-maxgap\t: max number of residues between consecutive seeds. ["      << SEED_MAXGAP << "]" << std::endl;
    std::cerr << "  --seed-mincnt\t: minimum number of seeds per sequence.["                   << SEED_MINCNT << "]" << std::endl;
    std::cerr << "  --sign-shared-perc-min\t: minimum number of signatures shared between query and target. ["       << SIGN_SHARED_PERC_MIN << "]" << std::endl;
    std::cerr << "  --sign-shared-perc-zscore\t: z-score of shared signatures below which further sequence search is skipped [" << SIGN_SHARED_PERC_ZSCORE << "]" << std::endl;
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

int seqlen_to_nsigns(int seqlen) {
    int ret = (seqlen - SIGN_LENGTH + 1) / HASH_MIN_RATIO;
    return MAX(ret, 1);
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
    uint32_t coveredcnt;
    uint32_t seqlen;
    uint16_t *signatures; // 11x less than seqlen
}
seq_t; // 16+16*4 bytes +++

static const int calc_perc_seq_id_editdist(const seq_t *seq1, const seq_t *seq2) {
    EdlibAlignResult result;
    result = edlibAlign(seq1->seq, seq1->seqlen, seq2->seq, seq2->seqlen,
                        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH));
    int editdist = result.editDistance;
    int alnlen = result.alignmentLength;
    edlibFreeAlignResult(result);
    int identities = alnlen - editdist;
    assert (identities <= seq1->seqlen);
    if (identities <= SIM_BASE) { return 0; }
    return floor((1 + DBL_EPSILON) * 100 * sqrt(SQUARE(identities) - SQUARE(SIM_BASE)) / seq1->seqlen);
}

static const int calc_perc_seq_my_editdist(const seq_t *seq1, const seq_t *seq2) {
    int maxEditDist = seq1->seqlen - ceil((1 - DBL_EPSILON) * (sqrt(SQUARE(SIM_BASE) + SQUARE(seq1->seqlen * SIM_PERC / 100.0))));
    if (maxEditDist < 0) { return 0; }
    
    EdlibAlignResult result;
    result = edlibAlign(seq1->seq, seq1->seqlen, seq2->seq, seq2->seqlen,
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
    int editdist = result.editDistance;
    edlibFreeAlignResult(result);

    assert(editdist >= -1);
    assert(editdist <= (int)seq1->seqlen);
    assert(editdist <= (int)seq2->seqlen);
    
    if (-1 == editdist) { return 0; }
    int ret = 100 * (seq1->seqlen - editdist) / seq1->seqlen;
    return ret;
}

static const int calc_perc_seq_sim_editdist(const seq_t *seq1, const seq_t *seq2) {
    if (1 == SIM_MODE) {
        return calc_perc_seq_my_editdist(seq1, seq2);
    }
    if (2 == SIM_MODE) {
        return calc_perc_seq_id_editdist(seq1, seq2);
    }
    abort();
}

int calc_n_shared_signatures(const seq_t *seq1, const seq_t *seq2) {
    int i = 0;
    int j = 0;
    int ret = 0;
    int nsignatures1 = seqlen_to_nsigns(seq1->seqlen);
    int nsignatures2 = seqlen_to_nsigns(seq2->seqlen);
    while (i != nsignatures1 && j != nsignatures2) {
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
    int dna_cnt = 0;
    int rna_cnt = 0;
    for (int i = 0; i < MIN(seq_arrlist.size, 100); i++) {
        for (int j = 0; j < seq_arrlist.data[i].seqlen; j++) {
            if (NULL != strchr("ACGTUacgtu", seq_arrlist.data[i].seq[j])) {
                nb_cnt++;
                if (NULL != strchr("Tt", seq_arrlist.data[i].seq[j])) {
                    dna_cnt++;
                } else if (NULL != strchr("Uu", seq_arrlist.data[i].seq[j])) {
                    rna_cnt++;
                }
            } else if (NULL == strchr("Nn-?*", seq_arrlist.data[i].seq[j])) {
                aa_cnt++;
            }
        }
    }
    if (aa_cnt * 4 > nb_cnt) {
        IS_INPUT_NUC = 0; // protein
        LEN_PERC = 80;
        SIM_BASE = 25;
        SIM_PERC = 50;
        SIM_MODE = 1;
    } else if (dna_cnt < rna_cnt) {
        IS_INPUT_NUC = 1; // rna
        LEN_PERC = 80;
        SIM_BASE = 25;
        SIM_PERC = 70;
        SIM_MODE = 1;
    } else {
        IS_INPUT_NUC = 2; // dna
        LEN_PERC = 0;
        SIM_BASE = 0;
        SIM_PERC = 90;
        SIM_MODE = 1;
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            ALPHA_TYPE_TO_CHAR_TO_REDUCED[j][i] = (char)i;
        }
    }
    std::vector<bool> are_args_parsed(argc+1);
    std:fill(are_args_parsed.begin(), are_args_parsed.end(), true);

    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--sim-mode", argv[i])) {
            SIM_MODE = atoi(argv[i+1]);
        } else if (!strcmp("--sim-perc", argv[i])) {
            SIM_PERC = atoi(argv[i+1]);
        } else if (!strcmp("--len-perc", argv[i])) {
            LEN_PERC = atoi(argv[i+1]);
        } else if (!strcmp("--is-input-nuc", argv[i])) {
            IS_INPUT_NUC = atoi(argv[i+1]);
        } else if (!strcmp("--dbentry-filt-occ-min", argv[i])) {
            DBENTRY_FILT_OCC_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--batch-size", argv[i])) {
            BATCH_SIZE = atoi(argv[i+1]);
        } else {
            are_args_parsed[i] = false;
            are_args_parsed[i+1] = false;
        }
    }

    //COV_SRC_MAX = (150 - SIM_PERC) / 20;
    //COV_SNK_MAX = (150 - SIM_PERC) / 5;

    SIGN_LENGTH = (SIM_PERC + 360) / (150 - SIM_PERC);
    SIGN_SHARED_PERC_MIN = 0; // ceil(MAX(0, SIM_PERC - sqrt(SIM_PERC * (100 - SIM_PERC) / 100.0) * 3))
     
    //ATTEMPT_INI = (120 - SIM_PERC); //* 16; // / 2;
    //ATTEMPT_INC = ATTEMPT_INI;
    //ATTEMPT_MAX = ATTEMPT_INI;
    //ATTEMPT_MIN = -2 * ATTEMPT_MAX / ATTEMPT_INC;

    if (IS_INPUT_NUC) {
        SIGN_LENGTH = (SIM_PERC + 900) / (200 - SIM_PERC);
    }
    
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--sim-base", argv[i])) {
            SIM_BASE = atoi(argv[i+1]);
        } else if (!strcmp("--cov-src-max", argv[i])) {
            COV_SRC_MAX = atoi(argv[i+1]);
        } else if (!strcmp("--cov-snk-max", argv[i])) {
            COV_SNK_MAX = atoi(argv[i+1]);
        } else if (!strcmp("--attempt-ini", argv[i])) {
            ATTEMPT_INI = atoi(argv[i+1]);
        } else if (!strcmp("--attempt-inc", argv[i])) {
            ATTEMPT_INC = atoi(argv[i+1]);
        } else if (!strcmp("--attempt-max", argv[i])) {
            ATTEMPT_MAX = atoi(argv[i+1]);
        }  else if (!strcmp("--attempt-min", argv[i])) {
            ATTEMPT_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--sign-shared-perc-min", argv[i])) {
            SIGN_SHARED_PERC_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--sign-shared-perc-zscore", argv[i])) {
            SIGN_SHARED_PERC_ZSCORE = atoi(argv[i+1]);
        } else if (!strcmp("--sign-length", argv[i])) {
            SIGN_LENGTH = atoi(argv[i+1]);
        } else if (!strcmp("--israndom", argv[i])) {
            // pass
        } else if (!are_args_parsed[i]) {
            if (strcmp("--seed-length", argv[i]) && strcmp("--seed-maxgap", argv[i]) && strcmp("--seed-mincnt", argv[i])) {
                show_usage(argc, argv);
            }
        } 
    }
    ATTEMPT_MAX = MAX(ATTEMPT_INI, ATTEMPT_MAX);

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
    int nsigns = seqlen_to_nsigns(seq_ptr->seqlen);
    seq_ptr->signatures = (uint16_t*) xcalloc(nsigns, sizeof(uint16_t));
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
        int j = 0;
        std::vector<uint16_t> signatures;
        signatures.reserve(nsigns);
        for (auto sign : signs) {
            signatures.push_back((uint16_t)sign);
            j++;
            if (nsigns == j) { break; }
        }
        std::sort(signatures.rbegin(), signatures.rend());
        j = 0;
        for (auto sign : signatures) {
            seq_ptr->signatures[j] = sign;
            j++;
        }
        assert(j <= nsigns);
    }
}

void seed_cov(const seed_t *seed, const uint32_t coveringidx, std::set<uint32_t> & visited) {
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    for (int i = 0; i < seed->size; i++) { 
        int coveredidx = seed->seqidxs[i];
        visited.insert(coveredidx);
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
    
    for (int i = 1; i < argc; i++) {
        if (!strcmp("--help", argv[i])) {
            show_usage(argc, argv);
        }
    }
    
    time_t begtime, endtime;

    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;
        
    std::set<int> printthresholds;
    for (int i = 0; i < 300; i++) {
        double thres = (double)((i+1)*(BATCH_SIZE+1)) * pow(1.05, (double)i);
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
    }
    kseq_destroy(kseq);
    
    PARAMS_init(argc, argv);
    showparams();

    // reinitialize some vars
    uint64_t num_residues = 0;
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        num_residues += seq_arrlist.data[i].seqlen;
    }
    DBENTRY_CNT = num_residues / RESIDUE_CNT_TO_SEED_CNT_RATIO + 1;
    std::cerr << "NUM_SEQS = " << seq_arrlist.size << std::endl;
    std::cerr << "NUM_RESIDUES = " << num_residues << std::endl;
    std::cerr << "DBENTRY_CNT = " << DBENTRY_CNT << std::endl;
    
    double INFO_PER_LETTER = (IS_INPUT_NUC ? 3.3 : 8.5);
    double adjusted_kmer_size = log((double)num_residues+INFO_PER_LETTER) / log(INFO_PER_LETTER);
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

    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--seed-length", argv[i])) { SEED_LENGTH = atoi(argv[i+1]); } 
        else if (!strcmp("--seed-maxgap", argv[i])) { SEED_MAXGAP = atoi(argv[i+1]); }
        else if (!strcmp("--seed-mincnt", argv[i])) { SEED_MINCNT = atoi(argv[i+1]); }
        
    }

    std::cerr << "After adjustment, (SEED_LENGTH, SEED_MAXGAP, SEED_MINCNT) = " << SEED_LENGTH << ", " << SEED_MAXGAP << "," << SEED_MINCNT << std::endl;

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

    #pragma omp parallel for schedule(dynamic, 9999*100) 
    for (int i = 0 ; i < DBENTRY_CNT; i++) {
        if (seeds[i].size > DBENTRY_FILT_OCC_MIN) {
            unsigned int randstate = 1 + (unsigned int)i;
            int probsimcnt = 0;
            std::array<std::set<std::pair<uint32_t, uint32_t>>, 101> psigns_to_seqidx_pairs;
            std::fill(psigns_to_seqidx_pairs.begin(), psigns_to_seqidx_pairs.end(), std::set<std::pair<uint32_t, uint32_t>>());
            for (int j = 0; j < DBENTRY_FILT_PAIR_SUBSAMP_CNT; j++) {
                uint32_t seqidx1 = seeds[i].seqidxs[rand_r(&randstate) % seeds[i].size];
                uint32_t seqidx2 = seeds[i].seqidxs[rand_r(&randstate) % (seeds[i].size - 1)];
                if (seqidx1 == seqidx2) { seqidx2 = seeds[i].size-1; }
                int nsharedsigns = calc_n_shared_signatures(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]);
                int percsharedsigns1 = 100 * nsharedsigns / seqlen_to_nsigns(seq_arrlist.data[seqidx1].seqlen);
                int percsharedsigns2 = 100 * nsharedsigns / seqlen_to_nsigns(seq_arrlist.data[seqidx2].seqlen);
                psigns_to_seqidx_pairs[percsharedsigns1].insert(std::make_pair(seqidx1, seqidx2));
                psigns_to_seqidx_pairs[percsharedsigns2].insert(std::make_pair(seqidx2, seqidx1));
            }
            int nhits = 0;
            int attemptcnt = 0;
            for (int psigns = 100; 
                    psigns >= SIGN_SHARED_PERC_MIN && attemptcnt < DBENTRY_FILT_PAIR_ATTEMPT_CNT && nhits <= DBENTRY_FILT_PAIR_TRUEHIT_MAX; 
                    psigns--, attemptcnt++) {
                for (auto seqidxpair : psigns_to_seqidx_pairs[psigns]) {
                    const seq_t *seq1 = &seq_arrlist.data[seqidxpair.first];
                    const seq_t *seq2 = &seq_arrlist.data[seqidxpair.second];
                    if (seq1->seqlen * LEN_PERC <= seq2->seqlen * 100 && seq2->seqlen * LEN_PERC <= seq1->seqlen * 100) { 
                        uint8_t sim = calc_perc_seq_sim_editdist(seq1, seq2);
                        if (sim >= SIM_PERC) { nhits++; }
                    }
                }
            }
            if (nhits <= DBENTRY_FILT_PAIR_TRUEHIT_MAX) {
                seeds[i].size = 0;
                seeds[i].bufsize = 0;
                free(seeds[i].seqidxs);
            }
        }
    }
    print_seedsize_histogram(seedsize_histogram, "seedsize_histogram 2"); 
    
    printf("%d %d\n", seq_arrlist.size, seq_arrlist.size);
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> coveredarr(BATCH_SIZE);
    std::fill(coveredarr.begin(), coveredarr.end(), std::vector<std::pair<uint32_t, uint8_t>>(0));

    unsigned int randstate = -1;
    time(&begtime);

    for (int64_t iter = 0; iter < seq_arrlist.size;) {
        int itermax = MIN(iter+BATCH_SIZE, seq_arrlist.size);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = iter; i < itermax; i++)
        {
            std::set<uint32_t> visited;
            int hash2size = 0;
            int filteredcnt = 0;
            int iteratedseqcnt = 0;
            int tp_distcompcnt = 0;
            int fp_distcompcnt = 0;
            int max_attempts = ATTEMPT_INI;
            
            int iteratedseqcnt_atmax = 0;
            int tp_distcompcnt_atmax = 0;
            int fp_distcompcnt_atmax = 0;
            if (seq_arrlist.data[i].coveredcnt <= COV_SRC_MAX && (int)(seq_arrlist.data[i].seqlen) >= (int)SEED_LENGTH) {
                uint64_t hash = hash_init(seq_arrlist.data[i].seq);
                seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                for (int j = SEED_LENGTH; j < strlen(seq_arrlist.data[i].seq); j++) {
                    hash = hash_update(hash, seq_arrlist.data[i].seq[j-SEED_LENGTH], seq_arrlist.data[i].seq[j]);
                    seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                }
                std::vector<std::vector<uint32_t>> psharedsigns_to_coveredidxs_vec(100+1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                    if ((i != coveredidx)) {
                        int n_shared_signatures = calc_n_shared_signatures(&seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
                        int psharedsigns = 100 * n_shared_signatures / seqlen_to_nsigns(seq_arrlist.data[coveredidx].seqlen);
                        psharedsigns_to_coveredidxs_vec.at(psharedsigns).push_back(coveredidx);
                    }
                }
                for (int psigns = 100; psigns >= SIGN_SHARED_PERC_MIN; psigns--) {
                    filteredcnt += psharedsigns_to_coveredidxs_vec.at(psigns).size();
                }
                int attempts = ATTEMPT_INI;
#define HASH2 0
#if HASH2
                // hash2size = 0;
                // for (auto seqidx : visited) { hash2size += seqlen_to_nsigns(seq_arrlist.data[seqidx].seqlen); }
                // hash2size = MIN((int)UINT16_MAX, hash2size);
                // hash2size = MIN(seqlen_to_nsigns(seq_arrlist.data[i].seqlen) * 4, hash2size);
                // hash2size = 1; // FIXME: remove this perf test
                hash2size = seqlen_to_nsigns(seq_arrlist.data[i].seqlen) * 4;

                int *hash2_to_attcnt = (int*) xmalloc(hash2size * sizeof(int));
                for (int i = 0; i < hash2size; i++) { hash2_to_attcnt[i] = ATTEMPT_INI; }
#else
                hash2size = -1;
#endif
                for (int psigns = 100; psigns >= SIGN_SHARED_PERC_MIN; psigns--) {
                    for (auto coveredidx : psharedsigns_to_coveredidxs_vec.at(psigns)) {
                        iteratedseqcnt++;
                        seq_t *coveringseq = &seq_arrlist.data[i];
                        seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                        double p = pow(SIM_PERC / 100.0, SIGN_LENGTH);
                        double n = seqlen_to_nsigns(coveredseq->seqlen);
                        int sd = ceil(sqrt(p * (1 - p) * n) / n * SIGN_SHARED_PERC_ZSCORE);
                        if (coveredseq->coveredcnt <= COV_SNK_MAX 
                                && coveringseq->seqlen * LEN_PERC <= coveredseq->seqlen * 100 
                                && coveredseq->seqlen * LEN_PERC <= coveringseq->seqlen * 100
                                && p <= psigns + sd) {
                            
                            int tot_attcnt = attempts;
#if HASH2
                            tot_attcnt = 0;
                            for (int i = 0; i < seqlen_to_nsigns(coveredseq->seqlen); i++) {
                                tot_attcnt += hash2_to_attcnt[coveredseq->signatures[i] % hash2size];
                            }
#endif

                            if (tot_attcnt > 0 
                                    // || seqlen_to_nsigns(seq_arrlist.data[coveredidx].seqlen) < 4
                                    ) {
                                uint8_t sim = calc_perc_seq_sim_editdist(coveredseq, coveringseq);
                                if (sim >= SIM_PERC) {
                                    coveredarr[i-iter].push_back(std::make_pair(coveredidx, sim));
                                    attempts += ATTEMPT_INC;
                                    UPDATE_MIN(attempts, ATTEMPT_MAX);
#if HASH2
                                    for (int i = 0; i < seqlen_to_nsigns(coveredseq->seqlen); i++) {
                                        hash2_to_attcnt[coveredseq->signatures[i] % hash2size] += ATTEMPT_INC;
                                        UPDATE_MIN(hash2_to_attcnt[coveredseq->signatures[i] % hash2size], ATTEMPT_MAX);
                                    }
#endif
                                    tp_distcompcnt++;
                                    if (attempts >= max_attempts) {
                                        max_attempts = attempts;
                                        tp_distcompcnt_atmax = tp_distcompcnt;
                                        fp_distcompcnt_atmax = fp_distcompcnt;
                                        iteratedseqcnt_atmax = iteratedseqcnt;
                                    }

                                } else {
                                    attempts -= 1;
#if HASH2
                                    for (int i = 0; i < seqlen_to_nsigns(coveredseq->seqlen); i++) {
                                        hash2_to_attcnt[coveredseq->signatures[i] % hash2size] -= 1;
                                        UPDATE_MAX(hash2_to_attcnt[coveredseq->signatures[i] % hash2size], ATTEMPT_MIN);
                                    }
#endif
                                    fp_distcompcnt++;
                                }
                            }
                        }
                    }
                }
#if HASH2
                free(hash2_to_attcnt);
#endif
                std::sort(coveredarr[i-iter].begin(), coveredarr[i-iter].end());
            }
            if (printthresholds.find(i) != printthresholds.end()) {
                time(&endtime);
                fprintf(stderr, "In %.f secs processed %d seqs\t"
                        "h1to4=%u/%i/%i/%u.\t"
                        "max_atts=%i at ic=%i,tp=%i,fp=%i\t"
                        "ATT_INI=%i,ATT_INC=%i\tseqlen=%u,coveredcnt=%u,hash2size=%d\n", 
                        difftime(endtime, begtime), i+1, 
                        coveredarr[i-iter].size(), tp_distcompcnt + fp_distcompcnt, filteredcnt, visited.size(), 
                        max_attempts, iteratedseqcnt_atmax, tp_distcompcnt_atmax, fp_distcompcnt_atmax,
                        ATTEMPT_INI, ATTEMPT_INC, seq_arrlist.data[i].seqlen, seq_arrlist.data[i].coveredcnt, hash2size); 
                assert(tp_distcompcnt == coveredarr[i-iter].size());
            }
        }
        for (int i = iter; i < itermax; i++) {
            printf("100 %d", i+1);
            int randmax = ceil(1000 * log(2) / log(1.0001 + coveredarr[i-iter].size()));
            for (auto adj : coveredarr[i-iter]) {
                int randval = rand_r(&randstate) % 1000;
                if (randval < randmax) { seq_arrlist.data[adj.first].coveredcnt++; }
                assert(adj.second <= 100);
                printf(" %d %u", (int)adj.second, adj.first+1);
            }
            printf("\n");
            coveredarr[i-iter].clear();
        }
        assert(coveredarr.size() == BATCH_SIZE);
        iter += BATCH_SIZE;
        for (int i = 0; i < 2; i++) {
            BATCH_SIZE++;
            coveredarr.push_back(std::vector<std::pair<uint32_t, uint8_t>>(0));
        }
    }
}

