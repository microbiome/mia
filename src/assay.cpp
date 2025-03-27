/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "assay.h"

#include <Rcpp.h>

using namespace su;

Assay::Assay(const Rcpp::NumericMatrix & assay){
    table = assay;
    
    sample_ids = std::vector<std::string>(); 
    obs_ids = std::vector<std::string>();
    
    Rcpp::StringVector rownames = Rcpp::rownames(table);
    obs_ids = Rcpp::as<std::vector<std::string>>(rownames);
    
    n_samples = table.ncol();
    n_obs = obs_ids.size();
    
    /* Define a mapping between an ID and its corresponding offset */
    obs_id_index = std::unordered_map<std::string, uint32_t>();
    
    create_id_index(obs_ids, obs_id_index);
    
    sample_counts = get_sample_counts();
}

Assay::~Assay(){
}

void Assay::create_id_index(std::vector<std::string> &ids, 
                            std::unordered_map<std::string, uint32_t> &map){
    uint32_t count = 0;
    map.reserve(ids.size());
    for( auto i = ids.begin(); i != ids.end(); i++, count++ ){
        map[*i] = count;
    }
}

std::vector<double> Assay::get_obs_data(const std::string &id) const {
    std::vector<double> out = std::vector<double>();
    uint32_t idx = obs_id_index.at(id);
    for( unsigned int i = 0; i < n_samples; i++ ){
        out.push_back(table(idx, i));
    }
    return out;
}

std::vector<double> Assay::get_sample_counts(){
    std::vector<double> sample_counts = std::vector<double>();
    
    for( unsigned int i = 0; i < n_samples; i++ ){
        unsigned int sum = 0;
        for( unsigned int j = 0; j < n_obs; j++ ){
            sum += table(j, i);
        }
        sample_counts.push_back(sum);
    }
    
    return(sample_counts);
}
