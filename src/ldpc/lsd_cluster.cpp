#include "ldpc/lsd_cluster.hpp"

#include <algorithm>
#include <iostream>
#include <random>

namespace ldpc::lsd {

std::vector<int> sort_indices(std::vector<double>& B) {
    std::vector<int> indices(B.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) { return B[i] < B[j]; });
    return indices;
}

void Statistics::clear() {
    this->individual_cluster_stats.clear();
    this->global_timestep_bit_history.clear();
    this->elapsed_time    = 0.0;
    this->lsd_method      = osd::OsdMethod::COMBINATION_SWEEP;
    this->lsd_order       = 0;
    this->bit_llrs        = {};
    this->error           = {};
    this->syndrome        = {};
    this->compare_recover = {};
}

[[nodiscard]] std::string Statistics::toString() {
    // build json like string object from individual cluster stats and global timestep bit history
    std::string result = "{";
    result += "\"elapsed_time_mu\":" + std::to_string(this->elapsed_time) + ",";
    result += "\"lsd_method\":" + std::to_string(static_cast<int>(this->lsd_method)) + ",";
    result += "\"lsd_order\":" + std::to_string(this->lsd_order) + ",";
    // print bit_llrs
    result += "\"bit_llrs\":[";
    for (auto i = 0; i < this->bit_llrs.size(); i++) {
        result += std::to_string(this->bit_llrs.at(i));
        if (i < this->bit_llrs.size() - 1) {
            result += ",";
        }
    }
    result += "],";
    result += "\"error\":[";
    for (auto i = 0; i < this->error.size(); i++) {
        result += std::to_string(this->error.at(i));
        if (i < this->error.size() - 1) {
            result += ",";
        }
    }
    result += "],";
    // print syndrome
    result += "\"syndrome\":[";
    for (auto i = 0; i < this->syndrome.size(); i++) {
        result += std::to_string(this->syndrome.at(i));
        if (i < this->syndrome.size() - 1) {
            result += ",";
        }
    }
    result += "],";
    // print compare_recover
    result += "\"compare_recover\":[";
    for (auto i = 0; i < this->compare_recover.size(); i++) {
        result += std::to_string(this->compare_recover.at(i));
        if (i < this->compare_recover.size() - 1) {
            result += ",";
        }
    }
    result += "],";
    result += "\"individual_cluster_stats\":{";
    for (auto& kv : this->individual_cluster_stats) {
        result += "\"" + std::to_string(kv.first) + "\":{";
        result += "\"active\":" + std::to_string(kv.second.active) + ",";
        result += "\"final_bit_count\":" + std::to_string(kv.second.final_bit_count) + ",";
        result += "\"undergone_growth_steps\":" + std::to_string(kv.second.undergone_growth_steps) + ",";
        result += "\"nr_merges\":" + std::to_string(kv.second.nr_merges) + ",";
        result += "\"got_valid_in_timestep\":" + std::to_string(kv.second.got_valid_in_timestep) + ",";
        result += "\"absorbed_by_cluster\":" + std::to_string(kv.second.absorbed_by_cluster) + ",";
        result += "\"got_inactive_in_timestep\":" + std::to_string(kv.second.got_inactive_in_timestep) + ",";
        result += "\"nr_of_non_zero_check_matrix_entries\":" +
                  std::to_string(kv.second.nr_of_non_zero_check_matrix_entries) + ",";
        result += "\"cluster_pcm_sparsity\":" + std::to_string(kv.second.cluster_pcm_sparsity) + ",";
        // print solution vector
        result += "\"solution\":[";
        for (auto i = 0; i < kv.second.solution.size(); i++) {
            result += std::to_string(kv.second.solution.at(i));
            if (i < kv.second.solution.size() - 1) {
                result += ",";
            }
        }
        result += "],";
        result += "\"size_history\":[";
        for (auto& s : kv.second.size_history) {
            result += std::to_string(s) + ",";
        }
        result.pop_back();
        result += "]},";
    }
    result.pop_back();
    result += "},";
    result += "\"global_timestep_bit_history\":{";
    for (auto& kv : this->global_timestep_bit_history) {
        result += "\"" + std::to_string(kv.first) + "\":{";
        for (auto& kv2 : kv.second) {
            result += "\"" + std::to_string(kv2.first) + "\":[";
            for (auto& b : kv2.second) {
                result += std::to_string(b) + ",";
            }
            // remove last , from result
            result.pop_back();
            result += "],";
        }
        result.pop_back();
        result += "},";
    }
    result.pop_back();
    result += "}";
    result += "}";
    return result;
}

LsdCluster::LsdCluster(ldpc::bp::BpSparse&                       parity_check_matrix,
                       int                                       syndrome_index,
                       std::shared_ptr<std::vector<LsdCluster*>> ccm,
                       std::shared_ptr<std::vector<LsdCluster*>> bcm,
                       bool                                      on_the_fly)
    : pcm(parity_check_matrix), cluster_id(syndrome_index), active(true), valid(false) {
    this->boundary_check_nodes.insert(syndrome_index);
    this->enclosed_syndromes.insert(syndrome_index);
    this->global_check_membership = ccm;
    this->global_bit_membership   = bcm;
    this->check_nodes.insert(syndrome_index);
    this->global_check_membership->at(syndrome_index) = this;
    this->cluster_pcm_syndrome.clear();
    this->pcm_check_idx_to_cluster_check_idx.insert(std::pair<int, int>{syndrome_index, 0});
    this->cluster_check_idx_to_pcm_check_idx.push_back(syndrome_index);
    this->nr_merges        = 0;
    this->pluDecomposition = ldpc::gf2dense::PluDecomposition(this->check_nodes.size(), this->bit_nodes.size(), this->cluster_pcm);
}

LsdCluster::~LsdCluster() {
    this->bit_nodes.clear();
    this->check_nodes.clear();
    this->boundary_check_nodes.clear();
    this->candidate_bit_nodes.clear();
    this->enclosed_syndromes.clear();
    this->merge_list.clear();
    this->cluster_pcm.clear();
    this->cluster_check_idx_to_pcm_check_idx.clear();
    this->pcm_check_idx_to_cluster_check_idx.clear();
    this->cluster_bit_idx_to_pcm_bit_idx.clear();
}

void LsdCluster::grow_cluster(const std::vector<double>& bit_weights, const int bits_per_step, const bool is_on_the_fly) {
    if (!this->active) {
        return;
    }

    this->compute_growth_candidate_bit_nodes();
    this->merge_list.clear();

    if (bit_weights == NULL_DOUBLE_VECTOR) {
        for (auto bit_index : this->candidate_bit_nodes) {
            this->add_bit_node_to_cluster(bit_index);
        }
    } else {
        std::vector<double> cluster_bit_weights;
        cluster_bit_weights.reserve(this->candidate_bit_nodes.size());
        for (auto bit : this->candidate_bit_nodes) {
            cluster_bit_weights.push_back(bit_weights[bit]);
        }

        auto sorted_indices             = sort_indices(cluster_bit_weights);
        auto candidate_bit_nodes_vector = std::vector<int>(this->candidate_bit_nodes.begin(), this->candidate_bit_nodes.end());

        if (this->is_randomize_same_weight_indices) {
            sorted_indices = randomize_same_weight_indices(sorted_indices, cluster_bit_weights);
        }

        int count = 0;
        for (auto i : sorted_indices) {
            if (count == bits_per_step) {
                break;
            }
            int bit_index = candidate_bit_nodes_vector[i];
            this->add_bit_node_to_cluster(bit_index);
            count++;
        }
    }

    this->merge_with_intersecting_clusters(is_on_the_fly);
}

std::vector<int> LsdCluster::randomize_same_weight_indices(const std::vector<int>&    sorted_indices,
                                                           const std::vector<double>& cluster_bit_weights) {
    if (cluster_bit_weights.empty() || sorted_indices.empty()) {
        return {};
    }

    auto reshuffeled_indices = std::vector<int>(sorted_indices.size());
    auto same_weight_indices = std::vector<int>();
    auto other_indices       = std::vector<int>();
    auto least_wt            = cluster_bit_weights[sorted_indices[0]];

    for (auto sorted_index : sorted_indices) {
        if (cluster_bit_weights[sorted_index] == least_wt) {
            same_weight_indices.push_back(sorted_index);
        } else {
            other_indices.push_back(sorted_index);
        }
    }

    if (same_weight_indices.size() > 1) {
        std::shuffle(same_weight_indices.begin(), same_weight_indices.end(), std::mt19937(std::random_device()()));
        for (size_t i = 0; i < same_weight_indices.size(); i++) {
            reshuffeled_indices[i] = same_weight_indices[i];
        }

        for (size_t i = 0; i < other_indices.size(); i++) {
            reshuffeled_indices[i + same_weight_indices.size()] = other_indices[i];
        }
    }

    return reshuffeled_indices;
}

void LsdCluster::merge_with_intersecting_clusters(const bool is_on_the_fly) {
    LsdCluster* larger = this;
    for (auto cl : merge_list) {
        larger = merge_clusters(larger, cl);
    }

    if (is_on_the_fly) {
        larger->valid = larger->apply_on_the_fly_elimination();
    }
}

void LsdCluster::compute_growth_candidate_bit_nodes() {
    std::vector<int> boundary_checks_to_erase;
    this->candidate_bit_nodes.clear();
    // we check for new candidate bit nodes as neighbours of boundary check nodes
    for (auto check_index : boundary_check_nodes) {
        bool erase = true;
        for (auto& e : this->pcm.iterate_row(check_index)) {
            // if bit is not in this cluster, add it to the candidate list.
            if (this->global_bit_membership->at(e.col_index) != this) {
                candidate_bit_nodes.insert(e.col_index);
                erase = false;
            }
        }
        // erase from boundary check nodes if no neighbouring bits are added to candidate list.
        if (erase) {
            boundary_checks_to_erase.push_back(check_index);
        }
    }
    for (auto check_index : boundary_checks_to_erase) {
        this->boundary_check_nodes.erase(check_index);
    }
}

bool LsdCluster::add_bit_node_to_cluster(const int bit_index, const bool in_merge) {
    auto bit_membership = this->global_bit_membership->at(bit_index);
    // if the bit is already in the cluster return.
    if (bit_membership == this) {
        // bit already in current cluster
        return false;
    }
    if (bit_membership == nullptr || in_merge) {
        // if the bit has not yet been assigned to a cluster or we are in merge mode we add it directly.
        this->add_bit(bit_index);
        if (this->global_timestep_bit_history != nullptr) {
            // add bit to timestep history with timestep the map size -1
            (*this->global_timestep_bit_history)[this->curr_timestep][this->cluster_id].push_back(bit_index);
        }
        // add incident checks as well, i.e., whole column to cluster pcm
        this->add_column_to_cluster_pcm(bit_index);
    } else {
        // if bit is in another cluster and we are not in merge mode, we add it to the merge list for later
        this->merge_list.insert(bit_membership);
    }

    return true;
}

LsdCluster* LsdCluster::merge_clusters(LsdCluster* cl1, LsdCluster* cl2) {
    LsdCluster* smaller = nullptr;
    LsdCluster* larger  = nullptr;
    if (cl1->bit_nodes.size() < cl2->bit_nodes.size()) {
        smaller = cl1;
        larger  = cl2;
    } else {
        smaller = cl2;
        larger  = cl1;
    }

    // we merge the smaller into the larger cluster
    for (auto bit_index : smaller->bit_nodes) {
        larger->add_bit_node_to_cluster(bit_index, true);
    }
    // check nodes are added with the bits
    for (auto check_index : smaller->boundary_check_nodes) {
        larger->boundary_check_nodes.insert(check_index);
    }
    for (auto j : smaller->enclosed_syndromes) {
        larger->enclosed_syndromes.insert(j);
    }
    smaller->active                   = false; // smaller, absorbed cluster is deactivated
    smaller->absorbed_into_cluster    = larger->cluster_id;
    smaller->got_inactive_in_timestep = smaller->curr_timestep;
    larger->nr_merges++;
    return larger;
}

int LsdCluster::add_check(const int check_index, const bool insert_boundary) {
    if (insert_boundary) {
        this->boundary_check_nodes.insert(check_index);
    }
    auto inserted = this->check_nodes.insert(check_index);
    if (!inserted.second) {
        this->global_check_membership->at(check_index) = this;
        return this->pcm_check_idx_to_cluster_check_idx[check_index];
    }

    this->global_check_membership->at(check_index) = this;
    this->cluster_check_idx_to_pcm_check_idx.push_back(check_index);
    int local_idx = this->cluster_check_idx_to_pcm_check_idx.size() - 1;
    this->pcm_check_idx_to_cluster_check_idx.insert(
            std::pair<int, int>{check_index, local_idx});
    return local_idx;
}

void LsdCluster::add_bit(const int bit_index) {
    auto inserted = this->bit_nodes.insert(bit_index);
    if (!inserted.second) {
        return;
    }
    this->global_bit_membership->at(bit_index) = this;
    // also add to cluster pcm
    this->cluster_bit_idx_to_pcm_bit_idx.push_back(bit_index);
}

void LsdCluster::add_column_to_cluster_pcm(const int bit_index) {
    std::vector<int> col;
    for (auto& e : this->pcm.iterate_column(bit_index)) {
        int  check_index      = e.row_index;
        auto check_membership = this->global_check_membership->at(check_index);
        if (check_membership == this) {
            // if already in cluster, add to cluster_pcm column of the bit and go to next
            // an index error on the map here indicates an error in the program logic.
            col.push_back(this->pcm_check_idx_to_cluster_check_idx[check_index]);
            continue;
        } else if (check_membership != nullptr) {
            // check is in another cluster
            this->merge_list.insert(check_membership);
        }
        // if check is in another cluster or none, we add it and update cluster_pcm
        auto local_idx = this->add_check(check_index, true);
        col.push_back(local_idx);
    }
    this->cluster_pcm.push_back(col);
}

bool LsdCluster::apply_on_the_fly_elimination() {
    // add columns to existing decomposition matrix
    // new bits are appended to cluster_pcm
    for (auto idx = pluDecomposition.col_count; idx < this->bit_nodes.size(); idx++) {
        this->pluDecomposition.add_column_to_matrix(this->cluster_pcm[idx]);
    }

    // convert cluster syndrome to dense vector fitting the cluster pcm dimensions for solving the system.
    // std::vector<uint8_t> cluster_syndrome;
    this->cluster_pcm_syndrome.resize(this->check_nodes.size(), 0);
    for (auto s : this->enclosed_syndromes) {
        this->cluster_pcm_syndrome[this->pcm_check_idx_to_cluster_check_idx.at(s)] = 1;
    }
    auto syndrome_in_image = this->pluDecomposition.rref_with_y_image_check(this->cluster_pcm_syndrome,
                                                                            pluDecomposition.cols_eliminated);
    return syndrome_in_image;
}

void LsdCluster::sort_non_pivot_cols(const std::vector<double>& bit_weights) {
    if (bit_weights.empty() or this->pluDecomposition.not_pivot_cols.size() < 2) {
        return;
    }
    // global index to cluster index mapping
    tsl::robin_map<int, int> global_to_cluster_colum_map;
    // local weight index to global bit index mapping
    tsl::robin_map<int, int> weight_idx_to_global_idx;
    // local weights to reorder
    std::vector<double> weights;
    for (auto not_pivot_col : this->pluDecomposition.not_pivot_cols) {
        auto global_idx                          = this->cluster_bit_idx_to_pcm_bit_idx[not_pivot_col];
        global_to_cluster_colum_map[global_idx]  = not_pivot_col;
        weight_idx_to_global_idx[weights.size()] = global_idx;
        weights.push_back(bit_weights.at(global_idx));
    }
    auto             sorted_indices = sort_indices(weights);
    std::vector<int> resorted_np_cls;
    resorted_np_cls.reserve(sorted_indices.size());
    // iterate over resorted indices
    // map local weight index to global bit index and then to cluster index
    for (auto idx : sorted_indices) {
        resorted_np_cls.push_back(global_to_cluster_colum_map[weight_idx_to_global_idx[idx]]);
    }
    this->pluDecomposition.not_pivot_cols = resorted_np_cls;
}

std::string LsdCluster::to_string() {
    int               count;
    std::stringstream ss{};
    ss << "........." << std::endl;
    ss << "Cluster ID: " << this->cluster_id << std::endl;
    ss << "Active: " << this->active << std::endl;
    ss << "Enclosed syndromes: ";
    for (auto i : this->enclosed_syndromes)
        ss << i << " ";
    ss << std::endl;
    ss << "Cluster bits: ";
    for (auto i : this->bit_nodes)
        ss << i << " ";
    ss << std::endl;
    ss << "Cluster checks: ";
    for (auto i : this->check_nodes)
        ss << i << " ";
    ss << std::endl;
    ss << "Candidate bits: ";
    for (auto i : this->candidate_bit_nodes)
        ss << i << " ";
    ss << std::endl;
    ss << "Boundary Checks: ";
    for (auto i : this->boundary_check_nodes)
        ss << i << " ";
    ss << std::endl;

    ss << "Cluster bit idx to pcm bit idx: ";
    count = 0;
    for (auto bit_idx : this->cluster_bit_idx_to_pcm_bit_idx) {
        ss << "{" << count << "," << bit_idx << "}";
        count++;
    }
    ss << std::endl;

    ss << "Cluster check idx to pcm check idx: ";
    count = 0;
    for (auto check_idx : this->cluster_check_idx_to_pcm_check_idx) {
        ss << "{" << count << "," << check_idx << "}";
        count++;
    }
    ss << std::endl;

    ss << "Cluster syndrome: ";
    for (auto check_idx : this->cluster_pcm_syndrome) {
        ss << unsigned(check_idx);
    }
    ss << std::endl;
    return ss.str();
}

} // namespace ldpc::lsd