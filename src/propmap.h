/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#ifndef __FAITH_PROPMAP
#define __FAITH_PROPMAP 1

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

std::vector<double> set_proportions(const BPTree &tree, uint32_t node,
                                    const Assay &table,
                                    PropMap &ps,
                                    bool normalize = true);
}

#endif /* __FAITH_PROPMAP */
