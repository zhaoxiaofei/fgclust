#include "kseq.h"
#include "edlib.h"

#include "parasail/matrices/blosum62.h"
#include "parasail.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

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

#define NUM_SEEDS (1000*1000*1000)
#define NUM_SIGNATURES 31

const uint64_t PRIME_BASE = 727L;
const uint64_t SIGN_BASE  = 10007L; // 1009L;
const uint64_t PRIME_MOD  = (0x1L << 31L) - 1L; //1000*1000*1000+7;
const uint64_t SIGN_MOD   = (0x1L << 31L) - 1L; //1000*1000*1000+7;

uint8_t PERC_SIM = 90; // 33;
int KMER_SIZE = 13; // 6; // 15;
int KMER_SPACE = 8; // 8; // 5
int SIGN_SIZE = 7;
int SIGN_MIN = 7; // 8-2; // 4;

int MAX_COV = 15;

int ATTEMPT_INI = 100;
int ATTEMPT_INC = 50;

char RED_ALPHA[256]; // derived

uint64_t PRIME_POWER = 0; // derived constant
uint64_t SIGN_POWER = 0; // derived constant

void showparams() {
    std::cerr << " NUM_SEEDS       = " << NUM_SEEDS      << std::endl;
    std::cerr << " NUM_SIGNATURES  = " << NUM_SIGNATURES << std::endl;

    std::cerr << " PERC_SIM    = " << (int)PERC_SIM << std::endl;
    std::cerr << " KMER_SIZE   = " << KMER_SIZE     << std::endl;
    std::cerr << " KMER_SPACE  = " << KMER_SPACE    << std::endl;
    std::cerr << " SIGN_SIZE   = " << SIGN_SIZE     << std::endl;
    std::cerr << " SIGN_MIN    = " << SIGN_MIN      << std::endl;
    std::cerr << " MAX_COV     = " << MAX_COV       << std::endl;
    std::cerr << " ATTEMPT_INI = " << ATTEMPT_INI   << std::endl;
    std::cerr << " ATTEMPT_INC = " << ATTEMPT_INC   << std::endl;
}

void show_usage(const int argc, const char *const *const argv) {
    std::cerr << "Program : " << argv[0] << std::endl;
    std::cerr << "Command-line arguments:" << std::endl;
    std::cerr << "--edsim\tA covers B if and only if (len(B) - edit-distance(A, B)) / len(B) >= edsim, where gaps at ends of B are not penalized." << PERC_SIM << std::endl;
    exit(-1);
}


void alphareduce(const char *const strarg) {
    const char *str = strarg;
    for (; *str; str++) {
        RED_ALPHA[*str] = *strarg;
    }
}

void PARAMS_init(const int argc, const char *const *const argv) {
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp("--edsim", argv[i])) {
            PERC_SIM = atoi(argv[i+1]);
        } else {
            show_usage(argc, argv);
        }
    }
    
    for (int i = 0; i < 256; i++) {
        RED_ALPHA[i] = (char)i;
    }

    if (PERC_SIM < 85) {
        KMER_SIZE = 12;
        KMER_SPACE = 7;
        SIGN_SIZE = 6;
        SIGN_MIN = 6;
        MAX_COV = 10;
        ATTEMPT_INI = 150;
        ATTEMPT_INC = 75;
        alphareduce("FY");
        alphareduce("ILMV");
        alphareduce("KR");
    } else if (PERC_SIM < 62) {
        KMER_SIZE = 11;
        KMER_SPACE = 6;
        SIGN_SIZE = 5;
        SIGN_MIN = 5;
        MAX_COV = 5;
        ATTEMPT_INI = 200;
        ATTEMPT_INC = 100;
        alphareduce("DENQ");
        alphareduce("FWY");
        alphareduce("ST");
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
        ret += (uint64_t) RED_ALPHA[beg[i]];
        ret %= PRIME_MOD; 
    }
    return ret;
}

const uint64_t hash_update(uint64_t hash, char prv, char nxt) {
    hash = hash * PRIME_BASE + (uint64_t)RED_ALPHA[nxt];
    hash += PRIME_MOD;
    hash -= (PRIME_POWER * (uint64_t)RED_ALPHA[prv]) % PRIME_MOD;
    return hash % PRIME_MOD;
}

const uint64_t sign_init(const char *beg) {
    uint64_t ret = 0;
    int i;
    for (i = 0; i < SIGN_SIZE; i++) {
        ret *= SIGN_BASE;
        ret += (uint64_t) RED_ALPHA[beg[i]];
        ret %= SIGN_MOD; 
    }
    return ret;
}

const uint64_t sign_update(uint64_t hash, char prv, char nxt) {
    hash = hash * SIGN_BASE + (uint64_t) RED_ALPHA[nxt];
    hash += SIGN_MOD;
    hash -= (SIGN_POWER * (uint64_t) RED_ALPHA[prv]) % SIGN_MOD;
    return hash % SIGN_MOD;
}

static const int calc_perc_seq_sim_editdist(const char *s1, const char *s2) {
    if (strlen(s1) * 80 > strlen(s2) * 100 || strlen(s2) * 80 > strlen(s1) * 100) { return 0; }
    int maxEditDist = strlen(s1) * (100 - PERC_SIM) / 100;
    EdlibAlignResult result;
    result = edlibAlign(s1, strlen(s1), s2, strlen(s2),
                        edlibNewAlignConfig(maxEditDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
    int editdist = result.editDistance;
    edlibFreeAlignResult(result);

    assert(editdist >= -1);
    assert(editdist <= (int)strlen(s1));
    assert(editdist <= (int)strlen(s2));
    
    if (-1 == editdist) { return 0; }
    int ret = 100 * (strlen(s1) - editdist) / ( strlen(s1) /*strlen(s1) + 18 */); 
    return ret;
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
    uint16_t coveredcnt;
}
seq_t; // 16+16*4 bytes +++

int are_probably_sim(const seq_t *seq1, const seq_t *seq2) {
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
            i++;
        } else if (seq1->signatures[i] > seq2->signatures[j]) {
            j++;
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
    strcpy(seq_arrlist.data[seq_arrlist.size].seq, kseq->seq.s);
    
    seq_arrlist.data[seq_arrlist.size].coveredcnt = 0;
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

    seq_arrlist.size++; 
    
    if ((int)kseq->seq.l >= (int)KMER_SIZE) {
        uint64_t hash = hash_init(kseq->seq.s);
        seed_add(hash, seq_arrlist.size - 1);
        for (i = KMER_SIZE; i < kseq->seq.l; i += 1) {
            hash = hash_update(hash, kseq->seq.s[i-KMER_SIZE], kseq->seq.s[i]);
            if (0 == i % KMER_SPACE) { seed_add(hash, seq_arrlist.size - 1); }
        }
    }
}

void seed_cover(const seed_t *seed, const uint32_t coveringidx, std::unordered_set<uint32_t> & visited, std::vector<std::pair<uint32_t, uint8_t>> & covered, int & filteredcnt) {
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    for (int i = 0; i < seed->size; i++) {
        int coveredidx = seed->seqidxs[i];
        if (visited.find(coveredidx) == visited.end()) { 
            seq_t *coveredseq = &seq_arrlist.data[coveredidx];
            if ((coveringidx != coveredidx)) {
                int probsim = are_probably_sim(&seq_arrlist.data[coveringidx], &seq_arrlist.data[coveredidx]);
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

void seed_cov(const seed_t *seed, const uint32_t coveringidx, std::unordered_set<uint32_t> & visited) {
    seq_t *coveringseq = &seq_arrlist.data[coveringidx];
    for (int i = 0; i < seed->size; i++) { 
        int coveredidx = seed->seqidxs[i];
        visited.insert(coveredidx);
    }
}

int main(const int argc, const char *const *const argv) {
    PARAMS_init(argc, argv);
    hash_sign_INIT(); 
    std::unordered_set<int> printthresholds;
    for (int i = 0; i < 400; i++) {
        printthresholds.insert((i+1)*1000*10 + (int)pow(1.05, (double)i));
    }

    seeds = (seed_t*) malloc(NUM_SEEDS * sizeof(seed_t));
    memset(seeds, 0, NUM_SEEDS * sizeof(seed_t));
    seq_arrlist_init();
    
    showparams();
    kseq_t *kseq = kseq_init(fileno(stdin));
    int i = 0;
    while ( kseq_read(kseq) >= 0 ) {
        seq_arrlist_add(kseq);
        i++;
        if (printthresholds.find(i) != printthresholds.end()) {
            fprintf(stderr, "Read and indexed %d sequences.\n", i);
        }
    }
    kseq_destroy(kseq);
    
    time_t begtime, endtime;
    time(&begtime);
    printf("%d %d\n", seq_arrlist.size, seq_arrlist.size);
    std::array<std::vector<std::pair<uint32_t, uint8_t>>, 9999> coveredarr;
    std::fill(coveredarr.begin(), coveredarr.end(), std::vector<std::pair<uint32_t, uint8_t>>(0));
    for (int iter = 0; iter < seq_arrlist.size; iter += 9999) {
        int itermax = MIN(iter+9999, seq_arrlist.size);
        
        #pragma omp parallel for schedule(dynamic, 1)
        for (int i = iter; i < itermax; i++)
        {
            std::unordered_set<uint32_t> visited;
            int filteredcnt = 0;
            if (seq_arrlist.data[i].coveredcnt <= MAX_COV && (int)strlen(seq_arrlist.data[i].seq) >= (int)KMER_SIZE) {
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
                        int probsim = are_probably_sim(&seq_arrlist.data[i], &seq_arrlist.data[coveredidx]);
                        nsharedsigns_to_coveredidxs_vec.at(probsim).push_back(coveredidx);
                    } 
                }
                
                int attempts = ATTEMPT_INI;
                for (int nsigns = NUM_SIGNATURES; nsigns >= SIGN_MIN && attempts > 0; nsigns--) {
                    for (auto coveredidx : nsharedsigns_to_coveredidxs_vec.at(nsigns)) {
                        seq_t *coveringseq = &seq_arrlist.data[i];
                        seq_t *coveredseq = &seq_arrlist.data[coveredidx];
                        uint8_t sim = calc_perc_seq_sim_editdist(coveredseq->seq, coveringseq->seq);
                        if (sim >= PERC_SIM) {
                            coveredarr[i-iter].push_back(std::make_pair(coveredidx, sim));
                            attempts += ATTEMPT_INC;
                        } else {
                            attempts--;
                        }
                        filteredcnt++;
                        if (attempts <= 0) { break; }
                    }
                }

                std::sort(coveredarr[i-iter].begin(), coveredarr[i-iter].end());
            }
            if (printthresholds.find(i) != printthresholds.end()) {
                time(&endtime);
                fprintf(stderr, "In %.f seconds processed %d sequences. h1/h2/h3 = %u/%u/%u.\n", difftime(endtime, begtime), i+1, coveredarr[i-iter].size(), filteredcnt, visited.size()); 
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
