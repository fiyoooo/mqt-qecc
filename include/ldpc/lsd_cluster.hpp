#pragma once

#include "bp.hpp"
#include "gf2dense.hpp"
#include "osd_dense.hpp"

#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <robin_map.h>
#include <robin_set.h>
#include <set>
#include <vector>

namespace ldpc::lsd {

const std::vector<double> NULL_DOUBLE_VECTOR = {};

// helper function to sort indices based on values in a vector
std::vector<int> sort_indices(std::vector<double>& B);

struct ClusterStatistics {
public:
    int                  final_bit_count                     = 0;     // nr of bits in 'final' cluster version, i.e., before solving for solution
    int                  undergone_growth_steps              = 0;     // nr of growth steps the cluster underwent
    int                  nr_merges                           = 0;     // nr of merges the cluster underwent
    std::vector<int>     size_history                        = {};    // history of cluster sizes from 0 to final bit count
    bool                 active                              = false; // if cluster is active, i.e., not merged into another cluster
    int                  got_valid_in_timestep               = -1;    // timestep in which cluster got valid
    int                  got_inactive_in_timestep            = -1;    // timestep in which cluster got inactive, i.e., was absorbed by another
    int                  absorbed_by_cluster                 = -1;    // cluster_id of the cluster the current one was merged into
    int                  nr_of_non_zero_check_matrix_entries = 0;     // nr of non zero entries in the cluster pcm
    double               cluster_pcm_sparsity                = 0;     // nr of non zero entries in the cluster pcm
    std::vector<uint8_t> solution{};                                  // local recovery, solution of cluster
};

struct Statistics {
    std::unordered_map<int, ClusterStatistics>                         individual_cluster_stats;    // clusterid <> stats
    std::unordered_map<int, std::unordered_map<int, std::vector<int>>> global_timestep_bit_history; // timestep <> (clusterid <> added bits)
    long                                                               elapsed_time{};
    osd::OsdMethod                                                     lsd_method;
    int                                                                lsd_order{};
    std::vector<double>                                                bit_llrs;
    std::vector<uint8_t>                                               error;           // the original error
    std::vector<uint8_t>                                               syndrome;        // the syndrome to decode
    std::vector<uint8_t>                                               compare_recover; // a recovery vector to compare against

    void clear();

    [[nodiscard]] std::string toString();
};

class LsdCluster {
public:
    LsdCluster();
    LsdCluster(ldpc::bp::BpSparse&                       parity_check_matrix,
               int                                       syndrome_index,
               std::shared_ptr<std::vector<LsdCluster*>> ccm,
               std::shared_ptr<std::vector<LsdCluster*>> bcm,
               bool                                      on_the_fly = false);

    ~LsdCluster();

    friend class LsdDecoder;

    /**
     * Grows the cluster by adding bit  nodes that are adjacent to boundary check nodes.
     * If bit_weights is provided, the bits are sorted by weight and only a single bit is added per growth step.
     * Otherwise, all bits adjacent to boundary check nodes are added.
     * @param bit_weights
     * @param bits_per_step
     * @return
     */
    void grow_cluster(const std::vector<double>& bit_weights   = NULL_DOUBLE_VECTOR,
                      const int                  bits_per_step = std::numeric_limits<int>::max(),
                      const bool                 is_on_the_fly = false);

    static std::vector<int> randomize_same_weight_indices(const std::vector<int>&    sorted_indices,
                                                          const std::vector<double>& cluster_bit_weights);

    /**
     * Merge this cluster with all clusters that intersect with it.
     * Keeps the larger cluster and merges the smaller cluster into it.
     * That is, the (reduced) parity check matrix of the larger cluster is kept.
     * After merging, the on-the-fly elimination is applied to the larger cluster.
     *
     * If on the fly elimination is applied true is returned if the syndrome is in the cluster.
     */
    void merge_with_intersecting_clusters(const bool is_on_the_fly = false);

    /**
     * Compute a list of candidate bit nodes to add to cluster as neighbours of boundary check nodes.
     * In case there are no new candidate bit nodes for a boundary check node, the check node is removed from the
     * boundary check node list.
     */
    void compute_growth_candidate_bit_nodes();

    /**
     * Adds a bit node to the cluster and updates all lists accordingly.
     * @param bit_index
     * @return true if the bit was added to the cluster, false otherwise.
     */
    bool add_bit_node_to_cluster(const int bit_index, const bool in_merge = false);

    /**
     * Merge this cluster with another cluster.
     * Keeps the larger cluster and merges the smaller cluster into it.
     * That is, the (reduced) parity check matrix of the larger cluster is kept.
     * @param cl2
     */
    static LsdCluster* merge_clusters(LsdCluster* cl1, LsdCluster* cl2);

    /**
     * Adds single check to cluster and updates all lists accordingly
     * @param check_index
     * @param insert_boundary
     */
    int add_check(const int check_index, const bool insert_boundary = false);

    /**
     * Adds single bit to cluster and updates all lists accordingly.
     * @param bit_index
     */
    void add_bit(const int bit_index);

    /**
     * Adds a column to the cluster parity check matrix.
     * The indices of the overall pcm are transferred to the local cluster pcm indices.
     * @param bit_index
     */
    void add_column_to_cluster_pcm(const int bit_index);

    /**
     * Apply on the fly elimination to the cluster.
     * Assumes only a single bit has been added.
     * The previously conducted row operations are applied to the new column, the new column is eliminated.
     * together with the syndrome.
     *
     * @return True if the syndrome is in the image of the cluster parity check matrix.
     */
    bool apply_on_the_fly_elimination();

    /**
     * Reorders the non-pivot columns of the eliminated cluster pcm according to the bit_weights.
     * @param bit_weights
     */
    void sort_non_pivot_cols(const std::vector<double>& bit_weights);

    std::string to_string();

    // TODO fix access, make private again soon
    bool                active{};
    bool                valid{};
    tsl::robin_set<int> bit_nodes;
    tsl::robin_set<int> check_nodes;
    tsl::robin_set<int> boundary_check_nodes;
    tsl::robin_set<int> candidate_bit_nodes;
    tsl::robin_set<int> enclosed_syndromes;

    std::shared_ptr<std::vector<LsdCluster*>> global_check_membership;
    std::shared_ptr<std::vector<LsdCluster*>> global_bit_membership;

    gf2dense::CscMatrix cluster_pcm;

    std::vector<uint8_t> cluster_pcm_syndrome;

    std::vector<int>         cluster_check_idx_to_pcm_check_idx;
    tsl::robin_map<int, int> pcm_check_idx_to_cluster_check_idx;

    std::vector<int>           cluster_bit_idx_to_pcm_bit_idx;
    gf2dense::PluDecomposition pluDecomposition;

private:
    ldpc::bp::BpSparse& pcm;
    int                 cluster_id{};

    tsl::robin_set<LsdCluster*> merge_list;

    int                                                       nr_merges{};
    std::unordered_map<int,
                       std::unordered_map<int,
                                          std::vector<int>>>* global_timestep_bit_history      = nullptr;
    int                                                       curr_timestep                    = 0;
    int                                                       absorbed_into_cluster            = -1;
    int                                                       got_inactive_in_timestep         = -1;
    bool                                                      is_randomize_same_weight_indices = false;
};

} // namespace ldpc::lsd
