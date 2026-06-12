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
#include "propmap.h"

#include <Rcpp.h>

using namespace su;

PropMap::PropMap(uint32_t vecsize)
    : prop_map()
    , defaultsize(vecsize)
{
    prop_map.reserve(1000);
}

PropMap::~PropMap(){}

std::vector<double> PropMap::get(uint32_t i){
    if( prop_map.count(i) > 0 ){
        return prop_map.at(i);
    } else {
        return(std::vector<double>());
    }
}

void PropMap::clear(uint32_t i){
    prop_map[i] = std::vector<double>();
}

void PropMap::update(uint32_t node, std::vector<double> vec){
    prop_map[node] = vec;
}

PropMapMulti::PropMapMulti(uint32_t _vecsize)
    : vecsize(_vecsize)
    , multi(get_num_stacks(), PropMap(DEF_VEC_SIZE)) {}

PropMapMulti::~PropMapMulti(){}

// Number of stacks = number of def_sizes that go in vecsize
// Rounding up ensures that there are always enough stacks for full vecsize
uint32_t PropMapMulti::get_num_stacks() const {
    return (vecsize + (DEF_VEC_SIZE-1)) / DEF_VEC_SIZE;
}

// get_start and get_end are used only for passing the values to set_prop_range
// and embed_prop_range
uint32_t PropMapMulti::get_start(uint32_t idx) const {
    return idx*DEF_VEC_SIZE;
}

uint32_t PropMapMulti::get_end(uint32_t idx) const {
    return std::min((idx+1)*DEF_VEC_SIZE, vecsize);
}

PropMap & PropMapMulti::get_prop_map(uint32_t idx){
    return multi[idx];
}

std::vector<double> su::set_proportions(const BPTree & tree,
                                        uint32_t node,
                                        const Assay & table,
                                        PropMap & pm,
                                        bool normalize){
    std::vector<double> props = std::vector<double>(table.n_samples, 0.0);
    if( tree.isleaf(node) ){
        std::string leaf = tree.names[node];
        props = table.get_obs_data(leaf); // get row for the specified node
        if( normalize ){
            for( unsigned int i = 0; i < table.n_samples; i++ ){
                props[i] /= table.sample_counts[i];
            }
        }
    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

        while( current <= right && current != 0 ){
            std::vector<double> vec = pm.get(current);  // Pull from prop map
            pm.clear(current);  // Remove from prop map

            for( unsigned int i = 0; i < table.n_samples; i++ ){
                props[i] = props[i] + vec[i];
            }

            current = tree.rightsibling(current);
        }
    }

    pm.update(node, props);
    return(props);
}

std::vector<double> su::set_proportions_range(const su::BPTree & tree,
                                              uint32_t node,
                                              const su::Assay & table,
                                              unsigned int start,
                                              unsigned int end,
                                              PropMap & pm,
                                              bool normalize){
    const unsigned int els = end-start;
    std::vector<double> props = std::vector<double>(els, 0.0);
    if(tree.isleaf(node)) {
        std::string leaf = tree.names[node];
        props = table.get_obs_data_range(leaf, start, end, normalize);
    } else {
        unsigned int current = tree.leftchild(node);
        unsigned int right = tree.rightchild(node);

        while(current <= right && current != 0) {
            std::vector<double> vec = pm.get(current);  // pull from prop map
            pm.clear(current);  // remove from prop map, place back on stack

            for(unsigned int i = 0; i < els; i++){
                props[i] = props[i] + vec[i];
            }

            current = tree.rightsibling(current);
        }
    }

    pm.update(node, props);
    return props;
}
