#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

int main(int argc, char** argv) {
    uint32_t nsets, nelems;
    
    std::string line;
    std::getline(std::cin, line);
    std::stringstream ss(line);
    ss >> nelems;
    ss >> nsets;
    // assert(nsets == nobjs);
    
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
        uint8_t sim;
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
        uint32_t nelems = set_to_elems[set].size();
        assert(0 < nelems);
        set_to_nelems[set] = nelems;
        nelems_to_sets[nelems].insert(set);
    }
    //assert(!std::getline(std::cin, line));

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
    
    for (uint32_t elem = 0; elem < nelems; elem++) {
        std::pair<uint32_t, uint8_t> coveringsetsim = std::make_pair(0, 0);
        for (auto setsim : elem_to_setsims[elem]) {
            if (set_to_iscovering[setsim.first] && setsim.second > coveringsetsim.second) {
                coveringsetsim = setsim;
            }
        }
        assert(coveringsetsim.second > 0);
        std::cout << (coveringsetsim.first+1) << "\t" << (elem+1) << "\t" << coveringsetsim.second << std::endl;
    }
}
