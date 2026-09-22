/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#ifndef __FAITH_PROPMAP_H
#define __FAITH_PROPMAP_H 1

#include <vector>
#include <stack>
#include <unordered_map>

#include "tree.h"
#include "assay.h"

namespace su {

class PropMap {
    public:
        PropMap(uint32_t vecsize);
        virtual ~PropMap();
        void clear(uint32_t i);
        void update(uint32_t i, std::vector<double> vec);
        std::vector<double> get(uint32_t i);
        
    private:
        std::unordered_map<uint32_t, std::vector<double>> prop_map;
        uint32_t defaultsize;
};

// Helper class that splits the full proportions vector into smaller chunks of
// pre-defined size
class PropMapMulti {
    public:
        PropMapMulti(uint32_t _vecsize);
        ~PropMapMulti();
        
        uint32_t get_num_stacks() const;
        uint32_t get_start(uint32_t idx) const;
        uint32_t get_end(uint32_t idx) const;
        PropMap & get_prop_map(uint32_t idx);
        
    private:
        const uint32_t vecsize; // Size of the full vector, equal to n_samples
        static const uint32_t DEF_VEC_SIZE = 1024; // size of the sub-vectors
        std::vector<PropMap> multi;
};

std::vector<double> set_proportions(const BPTree & tree, uint32_t node,
                                    const Assay & table,
                                    PropMap & pm,
                                    bool normalize = true);

std::vector<double> set_proportions_range(const su::BPTree & tree,
                                          uint32_t node,
                                          const su::Assay & table,
                                          unsigned int start,
                                          unsigned int end,
                                          PropMap & pm,
                                          bool normalize = true);

}

#endif /* __FAITH_PROPMAP_H */
