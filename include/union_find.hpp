#pragma once

#include "Code.hpp"
#include "MaxSATDecoder.hpp"
#include "ldpc/bp.hpp"
#include "ldpc/gf2sparse_linalg.hpp"
#include "ldpc/sparse_matrix_util.hpp"

#include <robin_map.h>
#include <robin_set.h>
#include <vector>

namespace ldpc::uf {

extern const std::vector<double> EMPTY_DOUBLE_VECTOR;
extern tsl::robin_set<int>       EMPTY_INT_ROBIN_SET;

std::vector<int> sort_indices(std::vector<double>& B);

struct Cluster {
    ldpc::bp::BpSparse&      pcm;
    tsl::robin_set<int>&     planar_code_boundary_bits;
    int                      cluster_id;
    bool                     contains_boundary_bits;
    bool                     active;
    bool                     valid;
    tsl::robin_set<int>      bit_nodes;
    tsl::robin_set<int>      check_nodes;
    tsl::robin_set<int>      boundary_check_nodes;
    std::vector<int>         candidate_bit_nodes;
    tsl::robin_set<int>      enclosed_syndromes;
    tsl::robin_map<int, int> spanning_tree_check_roots;
    tsl::robin_set<int>      spanning_tree_bits;
    tsl::robin_set<int>      spanning_tree_leaf_nodes;
    int                      spanning_tree_boundary_bit;
    Cluster**                global_check_membership;
    Cluster**                global_bit_membership;
    tsl::robin_set<Cluster*> merge_list;
    std::vector<int>         cluster_decoding;
    std::vector<int>         matrix_to_cluster_bit_map;
    tsl::robin_map<int, int> cluster_to_matrix_bit_map;
    std::vector<int>         matrix_to_cluster_check_map;
    tsl::robin_map<int, int> cluster_to_matrix_check_map;

    Cluster() = default;

    Cluster(ldpc::bp::BpSparse& parity_check_matrix, int syndrome_index, Cluster** ccm, Cluster** bcm,
            tsl::robin_set<int>& planar_code_boundary_bits = EMPTY_INT_ROBIN_SET);

    ~Cluster();

    int parity();

    void get_candidate_bit_nodes();

    int add_bit_node_to_cluster(int bit_index);

    void merge_with_cluster(Cluster* cl2);

    int grow_cluster(const std::vector<double>& bit_weights = EMPTY_DOUBLE_VECTOR, int bits_per_step = 0);

    int find_spanning_tree_parent(const int check_index);

    void find_spanning_tree();

    std::vector<int> peel_decode(const std::vector<uint8_t>& syndrome);

    ldpc::bp::BpSparse convert_to_matrix(const std::vector<double>& bit_weights = EMPTY_DOUBLE_VECTOR);

    std::vector<int> invert_decode(const std::vector<uint8_t>& syndrome, const std::vector<double>& bit_weights);

    std::vector<int> maxsat_decode(const std::vector<uint8_t>& syndrome, const std::vector<double>& bit_weights);

    void print();
};

class UfDecoder {
private:
    bool                weighted;
    ldpc::bp::BpSparse& pcm;

public:
    std::vector<uint8_t> decoding;
    int                  bit_count;
    int                  check_count;
    tsl::robin_set<int>  planar_code_boundary_bits;
    bool                 pcm_max_bit_degree_2;

    UfDecoder(ldpc::bp::BpSparse& parity_check_matrix);

    std::vector<uint8_t>&
    peel_decode(const std::vector<uint8_t>& syndrome, const std::vector<double>& bit_weights = EMPTY_DOUBLE_VECTOR,
                int bits_per_step = 1);

    std::vector<uint8_t>&
    matrix_decode(const std::vector<uint8_t>& syndrome,
                  const std::vector<double>&  bit_weights   = EMPTY_DOUBLE_VECTOR,
                  int                         bits_per_step = 1);

    std::vector<uint8_t>&
    maxsat_decode(const std::vector<uint8_t>& syndrome,
                  const std::vector<double>&  bit_weights   = EMPTY_DOUBLE_VECTOR,
                  int                         bits_per_step = 1);
};

} // namespace ldpc::uf
