/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree.h"

#include <stack>
#include <algorithm>

#include <Rcpp.h>

using namespace su;

BPTree::BPTree(std::vector<bool> input_structure, std::vector<double> input_lengths, std::vector<std::string> input_names) {
    
    structure = input_structure;
    lengths = input_lengths;
    names = input_names;
    
    nparens = structure.size();

    openclose = std::vector<uint32_t>();
    select_0_index = std::vector<uint32_t>();
    select_1_index = std::vector<uint32_t>();
    openclose.resize(nparens);
    select_0_index.resize(nparens / 2);
    select_1_index.resize(nparens / 2);
    excess.resize(nparens);

    structure_to_openclose();
    index_and_cache();
}

BPTree::BPTree(const Rcpp::List & rowTree) {
    
    //Initialize vectors
    openclose = std::vector<uint32_t>();
    lengths = std::vector<double>();
    names = std::vector<std::string>();
    excess = std::vector<uint32_t>();
    
    select_0_index = std::vector<uint32_t>();
    select_1_index = std::vector<uint32_t>();
    
    //Load the tree structure
    structure = std::vector<bool>();
    structure.reserve(500000);  // a fair sized tree... avoid reallocs, and its not _that_ much waste if this is wrong
    rowTree_to_bp(rowTree); //Also sets the size of nparens
    
    //Resize vectors
    // resize is correct here as we are not performing a push_back
    openclose.resize(nparens);
    lengths.resize(nparens);
    names.resize(nparens);
    excess.resize(nparens);
    
    select_0_index.resize(nparens / 2);
    select_1_index.resize(nparens / 2);
    
    //Builds a vector that lets us find the corresponding indices for each true/false pair
    structure_to_openclose();
    //Get metadata
    rowTree_to_metadata(rowTree);
    
    //Finalize
    index_and_cache();
}


BPTree BPTree::mask(std::vector<bool> topology_mask, std::vector<double> in_lengths) {
    
    std::vector<bool> new_structure = std::vector<bool>();
    std::vector<double> new_lengths = std::vector<double>();
    std::vector<std::string> new_names = std::vector<std::string>();

    uint32_t count = 0;
    for(auto i = topology_mask.begin(); i != topology_mask.end(); i++) {
        if(*i)
            count++;
    }

    new_structure.resize(count);
    new_lengths.resize(count);
    new_names.resize(count);

    auto mask_it = topology_mask.begin();
    auto base_it = this->structure.begin();
    uint32_t new_idx = 0;
    uint32_t old_idx = 0;
    for(; mask_it != topology_mask.end(); mask_it++, base_it++, old_idx++) {
        if(*mask_it) {
            new_structure[new_idx] = this->structure[old_idx];
            new_lengths[new_idx] = in_lengths[old_idx];
            new_names[new_idx] = this->names[old_idx];
            new_idx++;
        }
    }
    
    return BPTree(new_structure, new_lengths, new_names);
}

std::unordered_set<std::string> BPTree::get_tip_names() {
    std::unordered_set<std::string> observed;
	
    for(unsigned int i = 0; i < this->nparens; i++) {
        if(this->isleaf(i)) {
            observed.insert(this->names[i]);
        }
    }

    return observed;
}

BPTree BPTree::shear(std::unordered_set<std::string> to_keep) {
    std::vector<bool> shearmask = std::vector<bool>(this->nparens);
    int32_t p;

	for(unsigned int i = 0; i < this->nparens; i++) {
        if(this->isleaf(i) && to_keep.count(this->names[i]) > 0) {
            shearmask[i] = true;
            shearmask[i+1] = true;

            p = this->parent(i);
            while(p != -1 && !shearmask[p]) {
                shearmask[p] = true;
                shearmask[this->close(p)] = true;
                p = this->parent(p);
            }
        }
    }
    return this->mask(shearmask, this->lengths);
}

BPTree BPTree::collapse() {
    std::vector<bool> collapsemask = std::vector<bool>(this->nparens);
    std::vector<double> new_lengths = std::vector<double>(this->lengths);

    uint32_t current, first, last;

    for(uint32_t i = 0; i < this->nparens / 2; i++) {
        current = this->preorderselect(i);

        if(this->isleaf(current) or (current == 0)) {  // 0 == root
            collapsemask[current] = true;
            collapsemask[this->close(current)] = true;
        } else {
            first = this->leftchild(current);
            last = this->rightchild(current);

            if(first == last) {
                new_lengths[first] = new_lengths[first] + new_lengths[current];
            } else {
                collapsemask[current] = true;
                collapsemask[this->close(current)] = true;
            }
        }
    }

    return this->mask(collapsemask, new_lengths);
}

BPTree::~BPTree() {
}

void BPTree::index_and_cache() {
    // should probably do the open/close in here too
    unsigned int idx = 0;
    auto i = structure.begin();
    auto k0 = select_0_index.begin();
    auto k1 = select_1_index.begin();
    auto e_it = excess.begin();
    unsigned int e = 0;  
    
    for(; i != structure.end(); i++, idx++ ) {
        if(*i) {
            *(k1++) = idx;
            *(e_it++) = ++e;
        }
        else {
            *(k0++) = idx;
            *(e_it++) = --e;
        }
    }
}

uint32_t BPTree::postorderselect(uint32_t k) const { 
    return open(select_0_index[k]);
}

uint32_t BPTree::preorderselect(uint32_t k) const {
    return select_1_index[k];
}

inline uint32_t BPTree::open(uint32_t i) const {
    return structure[i] ? i : openclose[i];
}

inline uint32_t BPTree::close(uint32_t i) const {
    return structure[i] ? openclose[i] : i;
}

bool BPTree::isleaf(unsigned int idx) const {
    return (structure[idx] && !structure[idx + 1]);
}

uint32_t BPTree::leftchild(uint32_t i) const {
    // aka fchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return i + 1;
}

uint32_t BPTree::rightchild(uint32_t i) const {
    // aka lchild
    if(isleaf(i))
        return 0;  // this is awkward, using 0 which is root, but a root cannot be a child. edge case
    else
        return open(close(i) - 1);
}

uint32_t BPTree::rightsibling(uint32_t i) const {
    // aka nsibling
    uint32_t position = close(i) + 1;
    if(position >= nparens)
        return 0;  // will return 0 if no sibling as root cannot have a sibling
    else if(structure[position])
        return position;
    else 
        return 0;
}

int32_t BPTree::parent(uint32_t i) const {
    return enclose(i);
}

int32_t BPTree::enclose(uint32_t i) const {
    if(structure[i])
        return bwd(i, -2) + 1;
    else
        return bwd(i - 1, -2) + 1; 
}

int32_t BPTree::bwd(uint32_t i, int d) const {
    uint32_t target_excess = excess[i] + d;
    for(int current_idx = i - 1; current_idx >= 0; current_idx--) {
        if(excess[current_idx] == target_excess)
            return current_idx;
    }
    return -1;
}

// The algorithms that this class uses need the tree to be stored in a binary format
// In terms of the Newick format, an opening bracket corresponds to a TRUE,
// a closing bracket to a FALSE, and a tip to a TRUE FALSE
// This function assumes that the tree representation is in cladewise order -
// Ensure this with ape's reorder.phylo() function
// Seems to work whether or not the tree is rooted
void BPTree::rowTree_to_bp(const Rcpp::List & phylo) {
    Rcpp::NumericMatrix edge = phylo["edge"];
    Rcpp::StringVector tips = phylo["tip.label"];
    
    // phylo tips are always numbered from 1 to number of tips;
    uint32_t ntips = tips.size();
    
    // Keeps track of the branch's internal nodes
    std::stack<unsigned int> nodes; 
    
    unsigned int currentNode = 0;
    unsigned int nextNode = 0;
    
    // Goal: Insert true when a branch starts, a false when it closes, and a true-false for each tip.
    
    for (int i = 0; i < edge.nrow(); i++){
        currentNode = edge(i, 0);
        nextNode = edge(i, 1);
        
        if(nodes.size() > 0 && currentNode < nodes.top()) {
            // We've exhausted the branch and moved backwards in the tree
            do {
                nodes.pop();
                structure.push_back(false);
            } while(currentNode != nodes.top());
        }
        
        if(nodes.size() == 0 || currentNode > nodes.top() ) {
            // We are either at the root, or entering a new node
            nodes.push(currentNode);
            structure.push_back(true);
            
        }
        
        if(nextNode <= ntips) {
            // We've found a tip
            structure.push_back(true);
            structure.push_back(false);
        }
        
        if(i == edge.nrow() - 1) {
            // We've reached the end of the tree
            do {
                nodes.pop();
                structure.push_back(false);
            } while(nodes.size() > 0);
        }
    }
    nparens = structure.size();
}

void BPTree::structure_to_openclose() {
    std::stack<unsigned int> oc;
    unsigned int open_idx;
    unsigned int i = 0;

    for(auto it = structure.begin(); it != structure.end(); it++, i++) {
        if(*it) {
            oc.push(i);
        } else {
            open_idx = oc.top();
            oc.pop();
            openclose[i] = open_idx;
            openclose[open_idx] = i;
        }
    }
}

// Add metadata (lengths and names) to the tree representation
// Iterate through the structure, and whenever we hit a true decide if
// it's a leaf or not, and then add the corresponding label/length
// edge.length has (nodes + tips) elements - leaves at the start, nodes at the end
// tip.label has (tips) elements
// root.edge and node.labels are optional, giving the length of the root and
// the internal node (including root) labels, respectively
void BPTree::rowTree_to_metadata(const Rcpp::List & phylo) {
    Rcpp::NumericVector edgelength = phylo["edge.length"];
    Rcpp::NumericMatrix edges = phylo["edge"];
    Rcpp::StringVector tips = phylo["tip.label"];
    
    const uint32_t n_edges = edgelength.size();
    uint32_t ntips = tips.size();

    //Used to find the correct lengths for the nodes - Includes the root
    std::vector<double> edge_v(n_edges + 1, 0.0);
    
    for(unsigned int i = 0; i < n_edges; i++){
        edge_v.at(edges(i,1) - 1) = edgelength[i];
    }
    
    if(phylo.containsElementNamed("root.edge")) {
        edge_v.at(ntips) = phylo["root.edge"];
    }
    
    bool hasNodeLabels = false;
    Rcpp::StringVector nodes;
    
    if(phylo.containsElementNamed("node.labels")) {
        hasNodeLabels = true;
        nodes = phylo["node.labels"];
    }
    
    unsigned int tip_idx = 0; // tip indices run from 0 to ntips-1
    unsigned int node_idx = 0; // node indices run from ntips to ntips + nnodes - 1
    
    for(unsigned int i = 0; i < structure.size(); i++) {
        if(structure[i]){
            std::string label = std::string();
            double length = 0.0;
            
            if(isleaf(i)){
                //Tips can be expected to have both a length and a label
                label =  Rcpp::as<std::string>(tips[tip_idx]);
                length = edge_v[tip_idx];
                tip_idx++;
            }
            
            else{
                //Nodes always have lengths (except the root, which may have it optionally, but defaults to 0.0)
                //Nodes may also optionally have labels (which includes the root label)
                length = edge_v[ntips + node_idx];
                if(hasNodeLabels){
                    label = Rcpp::as<std::string>(nodes[node_idx]);
                }
                node_idx++;
            }
            set_node_metadata(i,label, length);
        }
    }
}

//This takes a label and a length and assigns them to the correct places
void BPTree::set_node_metadata(unsigned int open_idx, std::string name, double length) {
    names[open_idx] = name;
    lengths[open_idx] = length;
}

std::vector<bool> BPTree::get_structure() {
    return structure;
}

std::vector<uint32_t> BPTree::get_openclose() {
    return openclose;
}
