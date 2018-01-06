#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <omp.h>
#include <string.h>

// breadth first search
void bfs(uint32_t node, const std::vector<std::vector<uint32_t>> &node_to_adjs, int depth, std::vector<std::pair<uint32_t, uint8_t>> &elemsims) {
    std::set<uint32_t> visited;
    std::vector<uint32_t> currbatch;
    for (auto adj : node_to_adjs[node]) {
        visited.insert(adj);
        currbatch.push_back(adj);
    }
    while (currbatch.size() != 0 && depth > 0) {
        std::vector<uint32_t> nextbatch;
        for (uint32_t v : currbatch) {
            for (uint32_t adj : node_to_adjs[v]) {
                if (visited.find(adj) == visited.end()) {
                    elemsims.push_back(std::make_pair(adj, 1+depth));
                    visited.insert(adj);
                    nextbatch.push_back(adj);
                }
            }
        }
        currbatch.clear();
        for (auto v : nextbatch) {
            currbatch.push_back(v);
        }
        depth--;
    }
}

int main(int argc, char** argv) {
        
    int SETCOVER_DEPTH = 0;
    for (int i = 1; i < argc; i+=2) {
        if (i+1 < argc && !strcmp("--setcover-depth", argv[i])) {
            SETCOVER_DEPTH = atoi(argv[i+1]);
        } else {
            std::cerr << "Program : " << argv[0] << std::endl;
            std::cerr << "  version " << GITCOMMIT << " compiled by " << CXXVERSION << std::endl;
            std::cerr << "Command-line arguments with [default-values]:" << std::endl;
            std::cerr << "  --setcover-depth\t: A covers B if a path of at most this length from A to B exists. [" << SETCOVER_DEPTH << "]" << std::endl;
            exit(-1);
 
        }
    }
    std::cerr << "GITCOMMIT = "  << GITCOMMIT  << std::endl;
    std::cerr << "CXXVERSION = " << CXXVERSION << std::endl; 

    uint32_t nsets, nelems;
    
    std::string line;
    std::getline(std::cin, line);
    std::stringstream ss(line);
    ss >> nelems;
    ss >> nsets;
    
    std::vector<std::vector<uint32_t>>                     set_to_elems(nsets, std::vector<uint32_t>());
    std::vector<std::vector<std::pair<uint32_t, uint8_t>>> elem_to_setsims(nelems, std::vector<std::pair<uint32_t, uint8_t>>());
    std::vector<uint32_t>                                  set_to_nelems(nsets, 0);
    std::vector<std::set<uint32_t>>                        nelems_to_sets(nelems+1, std::set<uint32_t>());
    std::vector<bool>                                      set_to_iscovering(nsets, false);
    std::vector<bool>                                      elem_to_iscovered(nelems, false);
    
    for (uint32_t set = 0; set < nsets; set++) {
        std::getline(std::cin, line);
        std::stringstream ss(line);
        uint32_t tokval;
        uint32_t elem;
        uint8_t sim = 0;
        uint32_t tokcnt = 0;
        while (ss >> tokval) {
            tokcnt += 1;
            if (tokcnt % 2 == 0) {
                elem = tokval - 1;
                assert(0 <= elem);
                set_to_elems[set].push_back(elem);
                elem_to_setsims[elem].push_back(std::make_pair(set, sim));
            } else {
                sim = tokval;
            }
        }
        assert(tokcnt > 1 || !(std::cerr << "The line " << line << " is invalid!" << std::endl));
    }
    
    if (SETCOVER_DEPTH > 0) { 
        // transitive closure of sequence similarity somewhat work in practice, 
        // meaning if (A is similar to B and B is similar to C) then probably (A is similar to C).
        std::cerr << "Computing transitive closure in parallel.\n";
        // compute transitive closure
        std::vector<std::vector<std::pair<uint32_t, uint8_t>>> set_to_elemsims(nelems, std::vector<std::pair<uint32_t, uint8_t>>());
#pragma omp parallel for schedule(dynamic, 9999)
        for (uint32_t set = 0; set < nsets; set++) {
            std::vector<std::pair<uint32_t, uint8_t>> elemsims;
            bfs(set, set_to_elems, SETCOVER_DEPTH, elemsims);
            set_to_elemsims[set].insert(set_to_elemsims[set].end(), elemsims.begin(), elemsims.end());
        }
        std::cerr << "Converting transitive closure to apply to input.\n";
        for (uint32_t set = 0; set < nsets; set++) {
            for (auto elemsim : set_to_elemsims[set]) {
                auto elem = elemsim.first;
                auto sim = elemsim.second;
                assert(sim > 0);
                set_to_elems[set].push_back(elem);
                elem_to_setsims[elem].push_back(std::make_pair(set, sim));
            }
        }
    }
    
    for (uint32_t set = 0; set < nsets; set++) {
        uint32_t nelems = set_to_elems[set].size();
        assert (nelems <= elem_to_setsims.size());
        set_to_nelems[set] = nelems;
        nelems_to_sets[nelems].insert(set);
    }

    for (uint32_t tmpnelems = nelems; tmpnelems > 0; tmpnelems--) {
        std::unordered_set<uint32_t> nextsets;
        for (auto set : nelems_to_sets[tmpnelems]) {
            if (set_to_nelems[set] >= tmpnelems) {
                assert(set_to_nelems[set] == tmpnelems);
                uint32_t nuncovelems = 0;
                for (auto elem : set_to_elems[set]) { 
                    if (!elem_to_iscovered[elem]) {
                        nuncovelems++;
                    }
                }
                assert(tmpnelems == nuncovelems);
                for (auto elem : set_to_elems[set]) {
                    if (!elem_to_iscovered[elem]) {
                        for (auto setsim : elem_to_setsims[elem]) {
                            set_to_nelems[setsim.first]--;
                        }
                        elem_to_iscovered[elem] = true;
                    }
                }
                set_to_iscovering[set] = true;
            } else {
                assert(nextsets.find(set) == nextsets.end());
                nextsets.insert(set);
            }
        }

        for (auto set : nextsets) {
            assert(set_to_nelems[set] < tmpnelems);
            nelems_to_sets[set_to_nelems[set]].insert(set);
        }
    }
    
    std::unordered_set<uint32_t> coverings;
    for (uint32_t elem = 0; elem < nelems; elem++) {
        std::pair<uint32_t, uint8_t> coveringsetsim = std::make_pair(0, 0);
        for (auto setsim : elem_to_setsims[elem]) {
            if (set_to_iscovering[setsim.first] && setsim.second > coveringsetsim.second) {
                coveringsetsim = setsim;
            }
        }
        assert(coveringsetsim.second > 0);
        coverings.insert(coveringsetsim.first);
        std::cout << (coveringsetsim.first+1) << "\t" << (elem+1) << "\t" << (int)coveringsetsim.second << "\n";
    }
    std::cerr << "linsetcover is done : " << nelems <<  " elements are clustered into " << coverings.size() << " sets." << std::endl;
}

