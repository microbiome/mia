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

// Helper class
// To allow chunked processing, stores PropMap with vecsize-sized vectors
class PropMapMulti {
public:
    PropMapMulti(uint32_t _vecsize)
        : vecsize(_vecsize)
    , multi(get_num_stacks(), PropMap(DEF_VEC_SIZE)) // round up
    {}
    ~PropMapMulti(){}
    
    // Number of stacks = number of def_sizes that go in vecsize
    // Rounding up ensures that there are always enough stacks for full vecsize
    uint32_t get_num_stacks() const {return (vecsize + (DEF_VEC_SIZE-1)) / DEF_VEC_SIZE;}
    // These are used only for passing the value to set_prop_range and embed_prop_range
    uint32_t get_start(uint32_t idx) const {return idx*DEF_VEC_SIZE;}
    uint32_t get_end(uint32_t idx) const   {return std::min((idx+1)*DEF_VEC_SIZE, vecsize);}
    PropMap & get_prop_map(uint32_t idx) {return multi[idx];}
    
protected:
    const uint32_t vecsize; // equal to number of samples
    static const uint32_t DEF_VEC_SIZE = 1024; // size of the sub-vectors, small enough to fit in L1 cache
    std::vector<PropMap> multi; // Holds a StripeMap for each chunk
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
