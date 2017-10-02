#include "kseq.h"
#include "edlib.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <unordered_map>
#include <set>

#include <assert.h>
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

#define MIN(a, b) (((a) < (b) ? (a) : (b)))
#define MAX(a, b) (((a) > (b) ? (a) : (b)))
#define SQUARE(a) ((a)*(a))

#define NUM_SIGNATURES (32)

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

// variables that are initialzed from command line args

uint64_t BATCH_SIZE = 15999;
uint64_t RESIDUE_CNT_TO_SEED_CNT_RATIO = 40;

bool IS_INPUT_NUC = false;
int ALPHA_TYPE_TO_SIZE[] = {10, 10, 10};
int SIM_PERC = 90;
int SIM_BASE = 25;

int DBENTRY_FILT_OCC_MIN = 1000;
int DBENTRY_FILT_PAIR_SUBSAMP_CNT = 800;
int DBENTRY_FILT_PAIR_ATTEMPT_CNT = 10;
int DBENTRY_FILT_PAIR_TRUEHIT_MAX = 4;

// dbsize-related variables

uint64_t DBENTRY_CNT = 0; // can be overriden after determination of db size

// variables that are derived from SIM_PERC but can be overriden in command line

int COV_SRC_MAX = 0;
int COV_SNK_MAX = 0;

int SEED_LENGTH = 0; // can be overriden after determination of db size 
int SEED_MAXGAP = 0; // can be overriden after determination of db size
int SEED_MINCNT = 0; // can be overriden after determination of db size

int SIGN_LENGTH = 0;
int SIGN_SHARED_CNT_MIN = 0; 

int ATTEMPT_INI = 0;
int ATTEMPT_INC = 0;
int ATTEMPT_MAX = 0;


void showparams() {
    std::cerr << " BATCH_SIZE = " << BATCH_SIZE << std::endl;

    std::cerr << " RESIDUE_CNT_TO_SEED_CNT_RATIO = " << RESIDUE_CNT_TO_SEED_CNT_RATIO << std::endl;
    std::cerr << " NUM_SIGNATURES = " << NUM_SIGNATURES << std::endl;

    std::cerr << " IS_INPUT_NUC = " << IS_INPUT_NUC         << std::endl;
    std::cerr << " SIM_PERC = " << (int)SIM_PERC << std::endl;
    std::cerr << " SIM_BASE = " << (int)SIM_BASE << std::endl;
    
    std::cerr << " SEED_LENGTH = " << SEED_LENGTH     << std::endl;
    std::cerr << " SEED_MAXGAP = " << SEED_MAXGAP    << std::endl;
    std::cerr << " SIGN_LENGTH = " << SIGN_LENGTH     << std::endl;
    std::cerr << " SIGN_SHARED_CNT_MIN = " << SIGN_SHARED_CNT_MIN      << std::endl;
    std::cerr << " COV_SRC_MAX = " << COV_SRC_MAX  << std::endl;
    std::cerr << " COV_SNK_MAX = " << COV_SNK_MAX << std::endl;

    std::cerr << " ATTEMPT_INI = " << ATTEMPT_INI << std::endl;
    std::cerr << " ATTEMPT_INC = " << ATTEMPT_INC << std::endl;
    std::cerr << " ATTEMPT_MAX = " << ATTEMPT_MAX << std::endl;
    std::cerr << " DBENTRY_FILT_OCC_MIN = " << DBENTRY_FILT_OCC_MIN << std::endl;
    std::cerr << " DBENTRY_FILT_PAIR_SUBSAMP_CNT = " << DBENTRY_FILT_PAIR_SUBSAMP_CNT << std::endl;
    
    std::cerr << " DBENTRY_FILT_PAIR_TRUEHIT_MAX = " << DBENTRY_FILT_PAIR_TRUEHIT_MAX << std::endl;
    std::cerr << " DBENTRY_FILT_PAIR_ATTEMPT_CNT = " << DBENTRY_FILT_PAIR_ATTEMPT_CNT << std::endl;
}

void show_usage(const int argc, const char *const *const argv) {
    std::cerr << "Program : " << argv[0] << std::endl;
    std::cerr << "Command-line arguments:" << std::endl;
    std::cerr << "--edsim\tA covers B if and only if (len(B) - edit-distance(A, B)) / len(B) >= edsim, where gaps at ends of B are not penalized." << SIM_PERC << std::endl;
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

void PARAMS_init(const int argc, const char *const *const argv) {
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            ALPHA_TYPE_TO_CHAR_TO_REDUCED[j][i] = (char)i;
        }
    }
    
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--sim-perc", argv[i])) {
            SIM_PERC = atoi(argv[i+1]);
        } else if (!strcmp("--sim-base", argv[i])) {
            SIM_BASE = atoi(argv[i+1]);
        } else if (!strcmp("--is-input-nuc", argv[i])) {
            IS_INPUT_NUC = atoi(argv[i+1]);
        } 
        
        else if (!strcmp("--dbentry-filt-occ-min", argv[i])) {
            DBENTRY_FILT_OCC_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--sign-shared-cnt-min", argv[i])) {
            SIGN_SHARED_CNT_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--batch-size", argv[i])) {
            BATCH_SIZE = atoi(argv[i+1]);
        } else {
            show_usage(argc, argv);
        }
    }
    
    COV_SRC_MAX = (150 - SIM_PERC) / 20;
    COV_SNK_MAX = (150 - SIM_PERC) / 5;

    SIGN_LENGTH = (SIM_PERC + 360) / (150 - SIM_PERC);
    SIGN_SHARED_CNT_MIN = MAX(1, SIM_PERC / 10 - 4); 
    
    ATTEMPT_INI = 100 - SIM_PERC;
    ATTEMPT_INC = 100 - SIM_PERC;
    ATTEMPT_MAX = 110 - SIM_PERC;

    if (IS_INPUT_NUC) {
        SIGN_LENGTH = (SIM_PERC + 900) / (200 - SIM_PERC);
    }
    
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--cov-src-max", argv[i])) {
            COV_SRC_MAX = atoi(argv[i+1]);
        } else if (!strcmp("--cov-snk-max", argv[i])) {
            COV_SNK_MAX = atoi(argv[i+1]);
        } else {
            show_usage(argc, argv);
        }
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
    uint16_t signatures[NUM_SIGNATURES];
    uint32_t coveredcnt;
    uint32_t seqlen;
}
seq_t; // 16+16*4 bytes +++

static const int calc_perc_seq_sim_editdist(const seq_t *seq1, const seq_t *seq2) {
    
    if (seq1->seqlen * 80 > seq2->seqlen * 100 || seq2->seqlen * 80 > seq1->seqlen * 100) { return 0; }
    int maxEditDist = seq1->seqlen - ceil(sqrt(SQUARE((double)SIM_BASE) + SQUARE((double)(seq1->seqlen * SIM_PERC) / 100.0)));
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

int calc_n_shared_signatures(const seq_t *seq1, const seq_t *seq2) {
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
        seed->seqidxs = (uint32_t*) realloc(seed->seqidxs, seed->bufsize * sizeof(uint32_t));
    }
    seed->seqidxs[seed->size] = seqidx;
    seed->size++;
}

void seq_arrlist_init() {
    seq_arrlist.data = (seq_t*) malloc(4 * sizeof(seq_t));
    seq_arrlist.bufsize = 4;
    seq_arrlist.size = 0;
}

void seq_arrlist_add(const kseq_t *kseq) {
    if (seq_arrlist.size == seq_arrlist.bufsize) {
        seq_arrlist.bufsize *= 2;
        seq_arrlist.data = (seq_t*)realloc(seq_arrlist.data, seq_arrlist.bufsize * sizeof(seq_t));
    }
    seq_arrlist.data[seq_arrlist.size].name = (char*)malloc(kseq->name.l + 1);
    seq_arrlist.data[seq_arrlist.size].seq = (char*)malloc(kseq->seq.l + 1);
    strcpy(seq_arrlist.data[seq_arrlist.size].name, kseq->name.s);
    for (size_t i = 0; i <= kseq->seq.l; i++) { seq_arrlist.data[seq_arrlist.size].seq[i] = ALPHA_TYPE_TO_CHAR_TO_REDUCED[2][kseq->seq.s[i]]; }
    assert( seq_arrlist.data[seq_arrlist.size].seq[kseq->seq.l] == '\0' );
    seq_arrlist.data[seq_arrlist.size].seqlen = kseq->seq.l;
    seq_arrlist.data[seq_arrlist.size].coveredcnt = 0;
    seq_arrlist.size++;
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
    memset(seq_ptr->signatures, 0, NUM_SIGNATURES * sizeof(uint16_t));
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
    time_t begtime, endtime;

    std::cerr << "GITCOMMIT = " << GITCOMMIT << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl;
    PARAMS_init(argc, argv);
    std::set<int> printthresholds;
    for (int i = 0; i < 190; i++) {
        double thres = (double)((i+1)*(BATCH_SIZE+1)) * pow(1.05, (double)i);
        if (thres * 1.01 < (double)(INT_MAX)) { printthresholds.insert((int)thres); }
    }

    seq_arrlist_init();
    
    showparams();
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
        if (!strcmp("--kmer_size",     argv[i])) { SEED_LENGTH     = atoi(argv[i+1]); } 
        if (!strcmp("--kmer_space",    argv[i])) { SEED_MAXGAP    = atoi(argv[i+1]); }
        if (!strcmp("--min_num_kmers", argv[i])) { SEED_MINCNT = atoi(argv[i+1]); }
    }

    std::cerr << "After adjustment, (SEED_LENGTH, SEED_MAXGAP, SEED_MINCNT) = " << SEED_LENGTH << ", " << SEED_MAXGAP << "," << SEED_MINCNT << std::endl;

    hash_sign_INIT();

    seeds = (seed_t*) malloc(DBENTRY_CNT * sizeof(seed_t));
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
    }
    
    int seedsize_histogram[1000+1];
    print_seedsize_histogram(seedsize_histogram, "seedsize_histogram 1"); 

    #pragma omp parallel for schedule(dynamic, 9999*100) 
    for (int i = 0 ; i < DBENTRY_CNT; i++) {
        if (seeds[i].size > DBENTRY_FILT_OCC_MIN) {
            unsigned int randstate = 1 + (unsigned int)i;
            int probsimcnt = 0;
            std::array<std::set<std::pair<uint32_t, uint32_t>>, NUM_SIGNATURES+1> nsigns_to_seqidx_pairs;
            std::fill(nsigns_to_seqidx_pairs.begin(), nsigns_to_seqidx_pairs.end(), std::set<std::pair<uint32_t, uint32_t>>());
            for (int j = 0; j < DBENTRY_FILT_PAIR_SUBSAMP_CNT; j++) {
                uint32_t seqidx1 = seeds[i].seqidxs[rand_r(&randstate) % seeds[i].size];
                uint32_t seqidx2 = seeds[i].seqidxs[rand_r(&randstate) % (seeds[i].size - 1)];
                if (seqidx1 == seqidx2) { seqidx2 = seeds[i].size-1; }
                int nsharedsigns = calc_n_shared_signatures(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]); 
                nsigns_to_seqidx_pairs[nsharedsigns].insert(std::make_pair(seqidx1, seqidx2));
            }
            int nhits = 0;
            int attemptcnt = 0;
            for (int nsigns = NUM_SIGNATURES; 
                    nsigns >= SIGN_SHARED_CNT_MIN && attemptcnt < DBENTRY_FILT_PAIR_ATTEMPT_CNT && nhits <= DBENTRY_FILT_PAIR_TRUEHIT_MAX; 
                    nsigns--, attemptcnt++) {
                for (auto seqidxpair : nsigns_to_seqidx_pairs[nsigns]) {
                    int seqidx1 = seqidxpair.first;
                    int seqidx2 = seqidxpair.second;
                    uint8_t sim = calc_perc_seq_sim_editdist(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]);
                    if (sim >= SIM_PERC) { nhits++; }
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

    time(&begtime);

    for (int64_t iter = 0; iter < seq_arrlist.size; iter += BATCH_SIZE) {
        int itermax = MIN(iter+BATCH_SIZE, seq_arrlist.size);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = iter; i < itermax; i++)
        {
            std::set<uint32_t> visited;
            int filteredcnt = 0;
            int distcompcnt = 0;
            int max_attempts = ATTEMPT_INI;
            int max_attempts_arg = 0;
            if (seq_arrlist.data[i].coveredcnt <= COV_SRC_MAX && (int)(seq_arrlist.data[i].seqlen) >= (int)SEED_LENGTH) {
                uint64_t hash = hash_init(seq_arrlist.data[i].seq);
                seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                for (int j = SEED_LENGTH; j < strlen(seq_arrlist.data[i].seq); j++) {
                    hash = hash_update(hash, seq_arrlist.data[i].seq[j-SEED_LENGTH], seq_arrlist.data[i].seq[j]);
                    seed_cov(&seeds[hash % DBENTRY_CNT], i, visited);
                }
                std::vector<std::vector<uint32_t>> nsharedsigns_to_coveredidxs_vec(NUM_SIGNATURES + 1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    seq_t *coveredseq = &seq_arrlist.data[coveredidx];
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
                            uint8_t sim = calc_perc_seq_sim_editdist(coveredseq, coveringseq);

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
                assert(adj.second <= 100);
                printf(" %d %u", (int)adj.second, adj.first+1);
            }
            printf("\n");
            coveredarr[i-iter].clear();
        }
    }
}

