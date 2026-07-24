/*
* BSD 3-Clause License
*
* Copyright (c) 2016-2021, UniFrac development team.
* All rights reserved.
*
* See LICENSE file for more details
*/

#ifndef __FAITH_TREE_H
#define __FAITH_TREE_H 1

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_set>

#include <Rcpp.h>

namespace su {
class BPTree {
    public:
        /* Tracked attributes */
        std::vector<double> lengths;
        std::vector<std::string> names;
        
        /* Total number of parentheses */
        uint32_t nparens;
        
        /* constructor from a defined topology 
        *
        * @param input_structure A boolean vector defining the topology
        * @param input_lengths A vector of double of the branch lengths
        * @param input_names A vector of str of the vertex names
        */
        BPTree(std::vector<bool> input_structure,
                std::vector<double> input_lengths,
                std::vector<std::string> input_names);

        /* Constructor from a TreeSummarizedExperiment 
        *
        * @param treeSE An R treeSE object
        */
        BPTree(const Rcpp::List & rowTree);
        
        ~BPTree();

        /* Postorder tree traversal
        *
        * Get the index position of the ith node in a postorder tree
        * traversal.
        *
        * @param i The ith node in a postorder traversal
        */
        uint32_t postorderselect(uint32_t i)const ;

        /* Preorder tree traversal
        *
        * Get the index position of the ith node in a preorder tree
        * traversal.
        *
        * @param i The ith node in a preorder traversal
        */
        uint32_t preorderselect(uint32_t i) const;

        /* Test if the node at an index position is a leaf
        *
        * @param i The node to evaluate
        */
        bool isleaf(uint32_t i) const;

        /* Get the left child of a node
        *
        * @param i The node to obtain the left child from
        */
        uint32_t leftchild(uint32_t i) const ;

        /* Get the right child of a node
        *
        * @param i The node to obtain the right child from
        */
        uint32_t rightchild(uint32_t i) const;

        /* Get the right sibling of a node
        *
        * @param i The node to obtain the right sibling from
        */
        uint32_t rightsibling(uint32_t i) const;
        
        /* Get the parent of a node
        *
        * @param i The node to obtain the parent of
        */
        int32_t parent(uint32_t i) const;

        /* Get the names at the tips of the tree */
        std::unordered_set<std::string> get_tip_names();

        /* Public getters */
        std::vector<bool> get_structure();
        std::vector<uint32_t> get_openclose();
        
        // Mask self
        BPTree mask(std::vector<bool> topology_mask,
                    std::vector<double> in_lengths);

        BPTree shear(std::unordered_set<std::string> to_keep);

        BPTree collapse();

    private:
        std::vector<bool> structure;          // The topology
        std::vector<uint32_t> openclose;      // Cache'd mapping b/w parentheses
        std::vector<uint32_t> select_0_index; // Cache of select 0
        std::vector<uint32_t> select_1_index; // Cache of select 1
        std::vector<uint32_t> excess;

        /* Construct the select caches */
        void index_and_cache();
        /* Convert ape tree structure to boolean structure */
        void rowTree_to_bp(const Rcpp::List & phylo);
        /* Assign attributes */
        void rowTree_to_metadata(const Rcpp::List & phylo);
        /* Set the cache mapping between parentheses pairs */
        void structure_to_openclose();
        /* Set attributes for a node */
        void set_node_metadata(unsigned int open_idx,
                                std::string label, double length);
        /* Obtain the index of the opening for a given parenthesis */
        inline uint32_t open(uint32_t i) const;
        /* Obtain the index of the closing for a given parenthesis */
        inline uint32_t close(uint32_t i) const; 

        int32_t bwd(uint32_t i, int32_t d) const;
        int32_t enclose(uint32_t i) const;
    };
}

#endif /* __FAITH_TREE_H */
