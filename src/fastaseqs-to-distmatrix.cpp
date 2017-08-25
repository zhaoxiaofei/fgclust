#include "kseq.h"
#include "edlib.h"

#include "parasail/matrices/blosum62.h"
#include "parasail.h"

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


KSEQ_INIT(int, read)

#define MIN(a, b) (((a) < (b) ? (a) : (b)))
#define MAX(a, b) (((a) > (b) ? (a) : (b)))
#define UPDATE_MIN(a, b) { if ((a) > (b)) { (a) = (b); } } 
#define UPDATE_MAX(a, b) { if ((a) < (b)) { (a) = (b); } } 
#define SWAP(a, b) { (a)^=(b); (b)^=(a) ; (a)^=(b); }
#define SQUARE(a) ((a)*(a))

#define NUM_SIGNATURES (32)

const uint64_t PRIME_BASE = 48271L; // 727L;
const uint64_t SIGN_BASE  = 48271L; // 10007L; // 17001L; // 1009L;
const uint64_t PRIME_MOD  = (0x1L << 31L) - 1L; //1000*1000*1000+7;
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L; //1000*1000*1000+7;

int BATCH_SIZE = 15999;

const uint64_t NUM_RESIDUES_TO_NUM_SEEDS_RATIO = 40; // 33;

uint64_t NUM_SEEDS = 0;

bool ISNUC = false;

int SIM_ALPHASIZE = 20;

uint8_t PERC_SIM = 90; // 33;
uint8_t CANOPY_PERC_SIM = 95;

int BASE_SIM = 25; //30; // 26;
int KMER_SIZE = 7; // 12; //11; // 6; // 15;
int KMER_SPACE = 8; // 8; // 5
int MIN_NUM_KMERS = 70;

int SIGN_SIZE = 7;
int SIGN_MIN = 6; // 8-2; // 4;

int MAX_COV_AS_QUERY = 5; // 5; //10;
int MAX_COV_AS_TARGET = 20; // 1000*1000*1000;

int ATTEMPT_INI = 10;
int ATTEMPT_INC = 10;
int ATTEMPT_MAX = 20;
// int ATTEMPT_RES_PER_DEC = 900; // 1200;
// int ATTEMPT_LEARNING_RATE = 5;

int SEED_SIZE_MAX = 1000; // 3*300; // 30; // 1000;
int SEED_SEQIDX_PAIR_MAX = 800; //400;
// int SEED_SEQIDX_PAIR_PASS_SIGN_MIN_PERC = 25;

int SEED_COV_MAX = 1000*1000; // 1000*10;

int SEED_TRUE_HIT_MIN = 4;
int SEED_ATTEMPT_MIN = 10;

//int PRIOR_EDIT_DIST = 2;

char RED_ALPHA[3][256]; // derived

uint64_t PRIME_POWER = 0; // derived constant
uint64_t SIGN_POWER = 0; // derived constant

void showparams() {
    std::cerr << " BATCH_SIZE = " << BATCH_SIZE << std::endl;

    std::cerr << " NUM_RESIDUES_TO_NUM_SEEDS_RATIO = " << NUM_RESIDUES_TO_NUM_SEEDS_RATIO << std::endl;
    std::cerr << " NUM_SIGNATURES  = " << NUM_SIGNATURES << std::endl;

    std::cerr << " ISNUC       = " << ISNUC         << std::endl;
    std::cerr << " PERC_SIM    = " << (int)PERC_SIM << std::endl;
    std::cerr << " BASE_SIM    = " << (int)BASE_SIM << std::endl;
    
    std::cerr << " KMER_SIZE   = " << KMER_SIZE     << std::endl;
    std::cerr << " KMER_SPACE  = " << KMER_SPACE    << std::endl;
    std::cerr << " SIGN_SIZE   = " << SIGN_SIZE     << std::endl;
    std::cerr << " SIGN_MIN    = " << SIGN_MIN      << std::endl;
    std::cerr << " MAX_COV_AS_QUERY  = " << MAX_COV_AS_QUERY  << std::endl;
    std::cerr << " MAX_COV_AS_TARGET = " << MAX_COV_AS_TARGET << std::endl;

    std::cerr << " ATTEMPT_INI = " << ATTEMPT_INI   << std::endl;
    std::cerr << " ATTEMPT_INC = " << ATTEMPT_INC   << std::endl;
    std::cerr << " ATTEMPT_MAX = " << ATTEMPT_MAX   << std::endl;
    // std::cerr << " ATTEMPT_RES_PER_DEC = "   << ATTEMPT_RES_PER_DEC << std::endl;
    // std::cerr << " ATTEMPT_LEARNING_RATE = " << ATTEMPT_LEARNING_RATE << std::endl;
    std::cerr << " SEED_SIZE_MAX                       = " << SEED_SIZE_MAX                       << std::endl;
    std::cerr << " SEED_SEQIDX_PAIR_MAX                = " << SEED_SEQIDX_PAIR_MAX                << std::endl;
    // std::cerr << " SEED_SEQIDX_PAIR_PASS_SIGN_MIN_PERC = " << SEED_SEQIDX_PAIR_PASS_SIGN_MIN_PERC << std::endl;
    
    std::cerr << " SEED_TRUE_HIT_MIN                   = " << SEED_TRUE_HIT_MIN                 << std::endl;
    std::cerr << " SEED_ATTEMPT_MIN                    = " << SEED_ATTEMPT_MIN                  << std::endl;

    std::cerr << " SEED_COV_MAX = " << SEED_COV_MAX << std::endl;
    // std::cerr << " PRIOR_EDIT_DIST                     = " << PRIOR_EDIT_DIST                     << std::endl;
}

void show_usage(const int argc, const char *const *const argv) {
    std::cerr << "Program : " << argv[0] << std::endl;
    std::cerr << "Command-line arguments:" << std::endl;
    std::cerr << "--edsim\tA covers B if and only if (len(B) - edit-distance(A, B)) / len(B) >= edsim, where gaps at ends of B are not penalized." << PERC_SIM << std::endl;
    exit(-1);
}


void alphareduce(const char *const strarg, const int reducetype) {
    const char *str = strarg;
    for (; *str; str++) {
        RED_ALPHA[reducetype][*str] = *strarg;
    }
}

int calc_vecnorm(int a, int b) {
    return (int)ceil(sqrt(a * a + b * b));
}

void PARAMS_init(const int argc, const char *const *const argv) {
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--edsim", argv[i])) {
            PERC_SIM = atoi(argv[i+1]);
        } else if (!strcmp("--flatsim", argv[i])) {
            BASE_SIM = atoi(argv[i+1]);
        } else if (!strcmp("--kmersize", argv[i])) {
            KMER_SIZE = atoi(argv[i+1]);
        } else if (!strcmp("--isnuc", argv[i])) {
            ISNUC = atoi(argv[i+1]);
        } else if (!strcmp("--alphasize_sim", argv[i])) {
            SIM_ALPHASIZE = atoi(argv[i+1]);
        } else if (!strcmp("--seed_size_max", argv[i])) {
            SEED_SIZE_MAX = atoi(argv[i+1]);
        } else if (!strcmp("--sign_min", argv[i])) {
            SIGN_MIN = atoi(argv[i+1]);
        } else if (!strcmp("--batch_size", argv[i])) {
            BATCH_SIZE = atoi(argv[i+1]);
        } else {
            show_usage(argc, argv);
        }
    }
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {
            RED_ALPHA[j][i] = (char)i;
        }
    }
    
    CANOPY_PERC_SIM = (PERC_SIM + 100) / 2;

    if (ISNUC) {
        KMER_SIZE = calc_vecnorm(3 * (PERC_SIM + 10) / (110 - PERC_SIM), 9);
        SIGN_SIZE = calc_vecnorm(2 * (PERC_SIM + 10) / (110 - PERC_SIM), 6);
    } else {
        for (int t = 0; t < 3; t++) {
            alphareduce("DENQ", t);
            alphareduce("FWY", t);
            alphareduce("ILMV", t);
            alphareduce("KR", t);
            alphareduce("ST", t);
        }
        if (62 <= PERC_SIM && PERC_SIM < 81) {
            // KMER_SIZE = 12;
            // KMER_SPACE = 7;
            SIGN_SIZE = 6;
            SIGN_MIN = 5;
            // MAX_COV_AS_QUERY = 10;
            ATTEMPT_INI = 30; //50;
            ATTEMPT_INC = 30; //50;
            ATTEMPT_MAX = 40; //50;
            #if 0
            for (int t = 0; t < 3; t++) {
                alphareduce("FY", t);
                alphareduce("ILMV", t);
                //alphareduce("LVIM");
                alphareduce("KR", t);
            }
            #endif
        }
        if (10 <= PERC_SIM && PERC_SIM < 62) {
            // KMER_SIZE = 11;
            // KMER_SPACE = 6;
            SIGN_SIZE = 4; // 5
            SIGN_MIN = 1; // 4;
            // if (PERC_SIM < 45 || SIM_ALPHASIZE < 20) { SIGN_MIN = 1; }
            // MAX_COV_AS_QUERY = 5;
            ATTEMPT_INI = 50; // 125;
            ATTEMPT_INC = 50; // 125;
            ATTEMPT_MAX = 60; // 125;
            #if 0
            for (int t = 0; t < 3; t++) {
                alphareduce("DENQ", t);
                // alphareduce("EDNQ");
                alphareduce("FWY", t);
                alphareduce("ILMV", t);
                // alphareduce("LVIM");
                alphareduce("KR", t);
                alphareduce("ST", t);
            }
            #endif
        }
        assert(20 == SIM_ALPHASIZE || 15 == SIM_ALPHASIZE || 10 == SIM_ALPHASIZE || 0 == SIM_ALPHASIZE);
        if (15 == SIM_ALPHASIZE) {
            alphareduce("FY", 2);
            alphareduce("ILMV", 2);
            alphareduce("KR", 2);
        }
        if (10 == SIM_ALPHASIZE) {
            alphareduce("DENQ", 2);
            alphareduce("FWY", 2);
            alphareduce("ILMV", 2);
            alphareduce("KR", 2);
            alphareduce("ST", 2);
        }
    }
}

void hash_sign_INIT() {
    int i;
    PRIME_POWER = 1;
    for (i = 0; i < KMER_SIZE; i++) {
        PRIME_POWER = (PRIME_POWER * PRIME_BASE) % PRIME_MOD;
    }
    SIGN_POWER = 1;
    for (i = 0; i < SIGN_SIZE; i++) { 
        SIGN_POWER = (SIGN_POWER * SIGN_BASE) % SIGN_MOD;
    }
}

const uint64_t hash_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < KMER_SIZE; i++) {
        ret *= PRIME_BASE;
        ret += (uint64_t) RED_ALPHA[0][beg[i]];
        ret %= PRIME_MOD; 
    }
    return ret;
}

const uint64_t hash_update(uint64_t hash, char prv, char nxt) {
    hash = hash * PRIME_BASE + (uint64_t)RED_ALPHA[0][nxt];
    hash += PRIME_MOD;
    hash -= (PRIME_POWER * (uint64_t)RED_ALPHA[0][prv]) % PRIME_MOD;
    return hash % PRIME_MOD;
}

const uint64_t sign_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < SIGN_SIZE; i++) {
        ret *= SIGN_BASE;
        ret += (uint64_t) RED_ALPHA[1][beg[i]];
        ret %= SIGN_MOD; 
    }
    return ret;
}

const uint64_t sign_update(uint64_t hash, char prv, char nxt) {
    hash = hash * SIGN_BASE + (uint64_t) RED_ALPHA[1][nxt];
    hash += SIGN_MOD;
    hash -= (SIGN_POWER * (uint64_t) RED_ALPHA[1][prv]) % SIGN_MOD;
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

static const int calc_perc_seq_sim_editdist(const seq_t *seq1, const seq_t *seq2 /*, char redseqs[2][1024*64] */) {
    if (seq1->seqlen * 80 > seq2->seqlen * 100 || seq2->seqlen * 80 > seq1->seqlen * 100) { return 0; }
    int maxEditDist = seq1->seqlen - ceil(sqrt(SQUARE((double)BASE_SIM) + SQUARE((double)(seq1->seqlen * PERC_SIM) / 100.0)));
    if (maxEditDist < 0) { return 0; }
    // if (seq1->seqlen < BASE_SIM || seq2->seqlen < BASE_SIM) { return 0; }
    // int maxEditDist = MIN((int)(seq1->seqlen * (100 - PERC_SIM) / 100), (int)(seq1->seqlen - BASE_SIM)); // PRIOR_EDIT_DIST;
    // assert (maxEditDist >= 0);
    // if (maxEditDist < 0) { return 0; }
    EdlibAlignResult result;
    /*
    for (int i = 0; i < seq1->seqlen; i++) {
        redseqs[0][i] = RED_ALPHA[2][seq1->seq[i]];
    }
    for (int i = 0; i < seq2->seqlen; i++) {
        redseqs[1][i] = RED_ALPHA[2][seq2->seq[i]];
    }
    */
    result = edlibAlign(seq1->seq, seq1->seqlen, seq2->seq, seq2->seqlen,
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
    int editdist = result.editDistance;
    edlibFreeAlignResult(result);

    assert(editdist >= -1);
    assert(editdist <= (int)seq1->seqlen);
    assert(editdist <= (int)seq2->seqlen);
    
    if (-1 == editdist) { return 0; }
    int ret = 100 * (seq1->seqlen - editdist) / ( seq1->seqlen /*strlen(s1) + 18 */); 
    return ret;
}

int calc_n_shared_signatures(const seq_t *seq1, const seq_t *seq2) {
    //fprintf(stderr, "Comparing the sim btw %s and %s\n", seq1->name, seq2->name);
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

seed_t *seeds; //[NUM_SEEDS];

void seed_add(uint64_t hash, uint32_t seqidx) {
    seed_t *seed = &seeds[hash % NUM_SEEDS];
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

    size_t i;
    
    if (seq_arrlist.size == seq_arrlist.bufsize) {
        seq_arrlist.bufsize *= 2;
        seq_arrlist.data = (seq_t*)realloc(seq_arrlist.data, seq_arrlist.bufsize * sizeof(seq_t));
    }
    seq_arrlist.data[seq_arrlist.size].name = (char*)malloc(kseq->name.l + 1);
    seq_arrlist.data[seq_arrlist.size].seq = (char*)malloc(kseq->seq.l + 1);
    strcpy(seq_arrlist.data[seq_arrlist.size].name, kseq->name.s);
    for (i = 0; i <= kseq->seq.l; i++) { seq_arrlist.data[seq_arrlist.size].seq[i] = RED_ALPHA[2][kseq->seq.s[i]]; }
    assert( seq_arrlist.data[seq_arrlist.size].seq[kseq->seq.l] == '\0' );
    // strcpy(seq_arrlist.data[seq_arrlist.size].seq, kseq->seq.s);
    
    seq_arrlist.data[seq_arrlist.size].seqlen = kseq->seq.l;
    seq_arrlist.data[seq_arrlist.size].coveredcnt = 0;
    #if 0 
    memset(seq_arrlist.data[seq_arrlist.size].signatures, 0, NUM_SIGNATURES * sizeof(uint16_t));
    if ((int)kseq->seq.l >= (int)SIGN_SIZE) {
        std::vector<uint16_t> signs(kseq->seq.l * 8 * 0);
        //fprintf(stderr, "Start %s\n", kseq->seq.s);
        uint64_t sign = sign_init(kseq->seq.s);
        signs.push_back((uint16_t)sign);
        for (i = SIGN_SIZE; i < kseq->seq.l; i += 1) {
            sign = sign_update(sign, kseq->seq.s[i-SIGN_SIZE], kseq->seq.s[i]);
            // if ( kseq->seq.l - 1 == i ) fprintf(stderr, "Pushing %lu %u th time / %u\n", sign, i, kseq->seq.l);
            signs.push_back((uint16_t)sign);
        }
        #if 1
        std::sort(signs.begin(), signs.end());
        int j = 0;
        for (auto sign : signs) {
            seq_arrlist.data[seq_arrlist.size].signatures[j] = (uint16_t)sign;
            j++;
            if (NUM_SIGNATURES == j) { break; }
        }
        #endif
    }
    #endif
    seq_arrlist.size++;
}

void seq_longword_init(seq_t *const seq_ptr, int idx) {
    if ((int)seq_ptr->seqlen >= (int)KMER_SIZE) {
        int kmerspace = MAX(1, MIN(KMER_SPACE, seq_ptr->seqlen / MIN_NUM_KMERS));
        uint64_t hash = hash_init(seq_ptr->seq);
        seed_add(hash, idx);
        for (int i = KMER_SIZE; i < (int)seq_ptr->seqlen; i += 1) {
            hash = hash_update(hash, seq_ptr->seq[i-KMER_SIZE], seq_ptr->seq[i]);
            if (0 == i % kmerspace) { seed_add(hash, idx); }
        }
    }
}

void seq_signatures_init(seq_t *const seq_ptr) {
    memset(seq_ptr->signatures, 0, NUM_SIGNATURES * sizeof(uint16_t));
    if ((int)seq_ptr->seqlen >= (int)SIGN_SIZE) {
        std::vector<uint64_t> signs;
        signs.reserve((int)seq_ptr->seqlen - (int)SIGN_SIZE + 1);
        uint64_t sign = sign_init(seq_ptr->seq);
        signs.push_back(sign);
        for (int i = SIGN_SIZE; i < seq_ptr->seqlen; i += 1) {
            sign = sign_update(sign, seq_ptr->seq[i-SIGN_SIZE], seq_ptr->seq[i]);
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

#if 0
void seed_cover(const seed_t *seed, const uint32_t coveringidx, std::set<uint32_t> & visited, std::vector<std::pair<uint32_t, uint8_t>> & covered, int & filteredcnt) {
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    for (int i = 0; i < seed->size; i++) {
        int coveredidx = seed->seqidxs[i];
        if (visited.find(coveredidx) == visited.end()) { 
            seq_t *coveredseq = &seq_arrlist.data[coveredidx];
            if ((coveringidx != coveredidx)) {
                int probsim = calc_n_shared_signatures(&seq_arrlist.data[coveringidx], &seq_arrlist.data[coveredidx]);
                if (probsim >= SIGN_MIN) {
                    uint8_t sim = calc_perc_seq_sim_editdist(coveredseq->seq, coveringseq->seq); 
                    if (sim >= PERC_SIM) {
                        covered.push_back(std::make_pair(coveredidx, sim));
                    }
                    filteredcnt++;
                }
            }
            visited.insert(coveredidx);
        }
    }
}
#endif

void seed_cov(const seed_t *seed, const uint32_t coveringidx, std::set<uint32_t> & visited) {
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    if (seed->size <= SEED_COV_MAX) {
        for (int i = 0; i < seed->size; i++) { 
            int coveredidx = seed->seqidxs[i];
            visited.insert(coveredidx);
        }
    } else {
        unsigned int randstate = 1 + (unsigned int)coveringidx;
        for (int i = 0; i < seed->size; i++) {
            int coveredidx = seed->seqidxs[rand_r(&randstate)%seed->size];
            visited.insert(coveredidx);
        }
    }
}

void print_seedsize_histogram(int seedsize_histogram[], const char *name) {
    memset(seedsize_histogram, 0, (1000+1) * sizeof(int));
    for (int i = 0 ; i < NUM_SEEDS; i++) {
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
    PARAMS_init(argc, argv);
    // hash_sign_INIT(); 
    std::set<int> printthresholds;
    for (int i = 0; i < 190; i++) {
        double thres = (double)((i+1)*(BATCH_SIZE+1)) * pow(1.05, (double)i);
        if (thres * 1.01 < (double)(INT_MAX)) { printthresholds.insert((int)thres); }
    }

    seq_arrlist_init();
    
    showparams();
    kseq_t *kseq = kseq_init(fileno(stdin));
    int i = 0;
    while ( kseq_read(kseq) >= 0 ) {
        seq_arrlist_add(kseq);
        i++;
        if (printthresholds.find(i) != printthresholds.end()) {
            fprintf(stderr, "Read %d sequences.\n", i);
        }
    }
    kseq_destroy(kseq);
    
    // reinitialize some vars
    uint64_t num_residues = 0;
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        num_residues += seq_arrlist.data[i].seqlen;
    }
    NUM_SEEDS = num_residues / NUM_RESIDUES_TO_NUM_SEEDS_RATIO + 1; 
    std::cerr << "NUM_SEEDS = " << NUM_SEEDS << std::endl;

    double adjusted_kmer_size = MAX(log((double)num_residues+1.0) / log(8.5), 7.0); /// about 3 bits per reduced amino acid
    KMER_SIZE = MAX(KMER_SIZE, (int)floor(adjusted_kmer_size)); // heuristic estimate of kmer size
    int adjusted_kmer_space_ratio = (int)(100.0 * (adjusted_kmer_size - floor(adjusted_kmer_size)));
    
    KMER_SPACE = KMER_SPACE * (100 + adjusted_kmer_space_ratio) / 100; 
    MIN_NUM_KMERS = MIN_NUM_KMERS * 100 / (100 + adjusted_kmer_space_ratio);
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--kmersize", argv[i])) { KMER_SIZE = atoi(argv[i+1]); } 
    }

    std::cerr << "After adjustment, (KMER_SIZE, KMER_SPACE, MIN_NUM_KMERS) = " << KMER_SIZE << ", " << KMER_SPACE << "," << MIN_NUM_KMERS << std::endl;

    hash_sign_INIT();

    seeds = (seed_t*) malloc(NUM_SEEDS * sizeof(seed_t));
    memset(seeds, 0, NUM_SEEDS * sizeof(seed_t));
    
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        seq_longword_init(&seq_arrlist.data[i], i);
        if (printthresholds.find(i) != printthresholds.end()) {
            fprintf(stderr, "Indexed %d sequences.\n", i);
        }
    }

    #pragma omp parallel for schedule(dynamic, 9999*10)
    for (int i = 0 ; i < seq_arrlist.size; i++) {
        seq_signatures_init(&seq_arrlist.data[i]);
        /*
        if (printthresholds.find(i) != printthresholds.end()) {
            fprintf(stderr, "Minhashed %d sequences.\n", i);
        }
        */
    }
    
    int seedsize_histogram[1000+1];
    print_seedsize_histogram(seedsize_histogram, "seedsize_histogram 1"); 

    #pragma omp parallel for schedule(dynamic, 9999*100) 
    for (int i = 0 ; i < NUM_SEEDS; i++) {
        if (seeds[i].size > SEED_SIZE_MAX) {
            unsigned int randstate = 1 + (unsigned int)i;
            int probsimcnt = 0;
            std::array<std::set<std::pair<uint32_t, uint32_t>>, NUM_SIGNATURES+1> nsigns_to_seqidx_pairs;
            std::fill(nsigns_to_seqidx_pairs.begin(), nsigns_to_seqidx_pairs.end(), std::set<std::pair<uint32_t, uint32_t>>());
            for (int j = 0; j < SEED_SEQIDX_PAIR_MAX; j++) {
                uint32_t seqidx1 = seeds[i].seqidxs[rand_r(&randstate) % seeds[i].size];
                uint32_t seqidx2 = seeds[i].seqidxs[rand_r(&randstate) % (seeds[i].size - 1)];
                if (seqidx1 == seqidx2) { seqidx2 = seeds[i].size-1; }
                int nsharedsigns = calc_n_shared_signatures(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]); 
                nsigns_to_seqidx_pairs[nsharedsigns].insert(std::make_pair(seqidx1, seqidx2));
                /*
                if (nsharedsigns >= SIGN_MIN) {
                    probsimcnt += 1;
                }*/
            }
            int nhits = 0;
            int attemptcnt = 0;
            for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_MIN && attemptcnt < SEED_ATTEMPT_MIN && nhits <= 4; nsigns--, attemptcnt++) {
                for (auto seqidxpair : nsigns_to_seqidx_pairs[nsigns]) {
                    int seqidx1 = seqidxpair.first;
                    int seqidx2 = seqidxpair.second;
                    uint8_t sim = calc_perc_seq_sim_editdist(&seq_arrlist.data[seqidx1], &seq_arrlist.data[seqidx2]);
                    if (sim >= PERC_SIM) { nhits++; }
                }
            }
            if ( /*probsimcnt * 100 / SEED_SEQIDX_PAIR_MAX < SEED_SEQIDX_PAIR_PASS_SIGN_MIN_PERC*/ 
                    nhits <= SEED_TRUE_HIT_MIN // 4
                ) {
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
    
    // this constant is estimated empirically, the intuition is that clusters do not grow as fast as db
    uint64_t attempt_norm_fact = (uint64_t)(1000 * pow((double)(100*1000*1000) / (double)(seq_arrlist.size + 1), 0.2));
    int64_t edgecnt = 0;

    time_t begtime, endtime;
    time(&begtime);

    for (int64_t iter = 0; iter < seq_arrlist.size; iter += BATCH_SIZE) {
        int itermax = MIN(iter+BATCH_SIZE, seq_arrlist.size);

        // ensure that runtime is linear in size of output as well 
        int attempt_ini = ATTEMPT_INI; //MAX(5, ATTEMPT_INI * attempt_norm_fact * (1000*200 + 0 * iter) / (1000*200 + 0 * (iter + edgecnt)) / 2000);
        int attempt_inc = ATTEMPT_INC; //MAX(5, ATTEMPT_INC * attempt_norm_fact * (1000*200 + 0 * iter) / (1000*200 + 0 * (iter + edgecnt)) / 2000);
        int attempt_max = ATTEMPT_MAX; //MAX(5, ATTEMPT_MAX * attempt_norm_fact * (1000*200 + 0 * iter) / (1000*200 + 0 * (iter + edgecnt)) / 2000); 
        
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = iter; i < itermax; i++)
        {
            // char redseqs[2][1024*64];
            std::set<uint32_t> visited;
            int filteredcnt = 0;
            int distcompcnt = 0;
            int max_attempts = attempt_ini;
            int max_attempts_arg = 0;
            // int is_covered = i % MAX(seq_arrlist.data[i].coveredcnt, 1);
            if (seq_arrlist.data[i].coveredcnt <= MAX_COV_AS_QUERY && (int)(seq_arrlist.data[i].seqlen) >= (int)KMER_SIZE) {
                uint64_t hash = hash_init(seq_arrlist.data[i].seq);
                seed_cov(&seeds[hash % NUM_SEEDS], i, visited);
                // seed_cover(&seeds[hash % NUM_SEEDS], i, visited, coveredarr[i-iter], filteredcnt);
                for (int j = KMER_SIZE; j < strlen(seq_arrlist.data[i].seq); j++) {
                    hash = hash_update(hash, seq_arrlist.data[i].seq[j-KMER_SIZE], seq_arrlist.data[i].seq[j]);
                    seed_cov(&seeds[hash % NUM_SEEDS], i, visited);
                    // seed_cover(&seeds[hash % NUM_SEEDS], i, visited, coveredarr[i-iter], filteredcnt);
                }
                std::vector<std::vector<uint32_t>> nsharedsigns_to_coveredidxs_vec(NUM_SIGNATURES + 1, std::vector<uint32_t>());
                for (auto coveredidx: visited) {
                    seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                    if ((i != coveredidx)) {
                        int n_shared_signatures = calc_n_shared_signatures(&seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
                        nsharedsigns_to_coveredidxs_vec.at(n_shared_signatures).push_back(coveredidx);
                    }
                }
                for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_MIN; nsigns--) {
                    filteredcnt += nsharedsigns_to_coveredidxs_vec.at(nsigns).size();
                } 
                int attempts = attempt_ini;
                for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_MIN && attempts > 0 && attempts > max_attempts - ATTEMPT_MAX; nsigns--) {
                    for (auto coveredidx : nsharedsigns_to_coveredidxs_vec.at(nsigns)) {
                        seq_t *coveringseq = &seq_arrlist.data[i];
                        seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                        if (coveredseq->coveredcnt <= MAX_COV_AS_TARGET) {
                            uint8_t sim = calc_perc_seq_sim_editdist(coveredseq, coveringseq);
                            // uint8_t sim = calc_perc_seq_sim_editdist(coveredseq, coveringseq, redseqs);

                            if (sim >= PERC_SIM) {
                                coveredarr[i-iter].push_back(std::make_pair(coveredidx, sim));
                                attempts += attempt_inc;
                                if (attempts > max_attempts) {
                                    max_attempts = attempts;
                                    max_attempts_arg = distcompcnt + 1;
                                }
                            } else {
                                attempts -= 1 ; // + coveringseq->seqlen / ATTEMPT_RES_PER_DEC;
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
                fprintf(stderr, "In %.f seconds processed %d sequences.\th1to4=%u/%i/%i/%u.\tmax_attempts=%i at %i\tattempt_ini=%i\tattempt_inc=%i\tseqlen=%u\tcoveredcnt=%u\n", 
                        difftime(endtime, begtime), i+1, coveredarr[i-iter].size(), distcompcnt, filteredcnt, visited.size(), max_attempts, max_attempts_arg, 
                        attempt_ini, attempt_inc, seq_arrlist.data[i].seqlen, seq_arrlist.data[i].coveredcnt); 
            }
        }
        for (int i = iter; i < itermax; i++) {
            printf("100 %d", i+1);
            for (auto adj : coveredarr[i-iter]) {
                // if (adj.second >= CANOPY_PERC_SIM) {
                seq_arrlist.data[adj.first].coveredcnt++;
                // }
                assert(adj.second <= 100);
                printf(" %d %u", (int)adj.second, adj.first+1);
            }
            printf("\n");
            edgecnt += coveredarr[i-iter].size();
            coveredarr[i-iter].clear();
        }
        #if 0
        attempt_ini = (attempt_ini * (100 - ATTEMPT_LEARNING_RATE) + 400 * (itermax - iter) * ATTEMPT_LEARNING_RATE / coveredtotal) / 100;
        attempt_ini = MAX(1, MIN(ATTEMPT_INI, attempt_ini));
        attempt_inc = (attempt_inc * (100 - ATTEMPT_LEARNING_RATE) + 400 * (itermax - iter) * ATTEMPT_LEARNING_RATE / coveredtotal) / 100;
        attempt_inc = MAX(1, MIN(ATTEMPT_INC, attempt_inc));
        #endif
    }
}

