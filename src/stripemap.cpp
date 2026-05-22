/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#include "tree.h"
#include "assay.h"
#include "stripemap.h"

using namespace su;

StripeMap::StripeMap(uint32_t n_samples)
    : stripe_map()
    , vecsize(n_samples)
{
    n_stripes = (n_samples + 1) / 2;
    for( unsigned int i = 0; i < n_stripes; i++ ){
        this->update(i, std::vector<double>(vecsize, 0.0));
    }
}

StripeMap::~StripeMap(){
}

std::vector<double> StripeMap::get(uint32_t i){
    if( stripe_map.count(i) > 0 ){
        return stripe_map.at(i);
    } else {
        return(std::vector<double>());  
    }
}

void StripeMap::clear(uint32_t i){
    stripe_map[i] = std::vector<double>();
}

void StripeMap::update(uint32_t node, std::vector<double> vec){
    stripe_map[node] = vec;
}

bool StripeMap::is_empty(uint32_t i){
    return get(i).empty();
}