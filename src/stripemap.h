/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#ifndef __UNIFRAC_STRIPEMAP_H
#define __UNIFRAC_STRIPEMAP_H 1

#include <vector>
#include <stack>
#include <unordered_map>
#include <cstdint>

namespace su {

class StripeMap {
    public:
        StripeMap(uint32_t n_samples);
        virtual ~StripeMap();
        
        void clear(uint32_t i);
        void update(uint32_t i, std::vector<double> vec);
        std::vector<double> get(uint32_t i);
        bool is_empty(uint32_t i);
        
    private:
        std::unordered_map<uint32_t, std::vector<double>> stripe_map;
        uint32_t n_samples; 
        uint32_t n_stripes;
};

}

#endif /* __UNIFRAC_STRIPEMAP_H */
