#include "ldpc/lsd.hpp"

#include "ldpc/gf2dense.hpp"
#include "ldpc/sparse_matrix_util.hpp"

namespace ldpc::lsd {

void LsdDecoder::reset_cluster_stats() {
    this->statistics.clear();
}

osd::OsdMethod LsdDecoder::getLsdMethod() const {
    return lsd_method;
}

void LsdDecoder::setLsdMethod(osd::OsdMethod lsdMethod) {
    lsd_method = lsdMethod;
}

void LsdDecoder::set_do_stats(const bool on) {
    this->do_stats = on;
}

bool LsdDecoder::get_do_stats() const {
    return this->do_stats;
}

void LsdDecoder::print_cluster_stats() {
    std::cout << this->statistics.toString() << std::endl;
}

void LsdDecoder::update_growth_stats(const LsdCluster* cl) {
    this->statistics.individual_cluster_stats[cl->cluster_id].undergone_growth_steps++;
    this->statistics.individual_cluster_stats[cl->cluster_id].size_history.push_back(cl->bit_nodes.size());
    this->statistics.individual_cluster_stats[cl->cluster_id].active                   = true;
    this->statistics.individual_cluster_stats[cl->cluster_id].absorbed_by_cluster      = cl->absorbed_into_cluster;
    this->statistics.individual_cluster_stats[cl->cluster_id].got_inactive_in_timestep = cl->got_inactive_in_timestep;
}

void LsdDecoder::update_final_stats(const LsdCluster* cl) {
    this->statistics.individual_cluster_stats[cl->cluster_id].final_bit_count = cl->bit_nodes.size();
    this->statistics.individual_cluster_stats[cl->cluster_id].nr_merges       = cl->nr_merges;
    int nr_nonzero_elems                                                      = gf2dense::count_non_zero_matrix_entries(cl->cluster_pcm);
    this->statistics.individual_cluster_stats[cl->cluster_id].nr_of_non_zero_check_matrix_entries =
            nr_nonzero_elems;
    auto size = cl->pluDecomposition.col_count * cl->pluDecomposition.row_count;
    if (size > 0) {
        this->statistics.individual_cluster_stats[cl->cluster_id].cluster_pcm_sparsity =
                1.0 - ((static_cast<double>(nr_nonzero_elems)) / static_cast<double>(size));
    }
}

std::vector<uint8_t>& LsdDecoder::on_the_fly_decode(std::vector<uint8_t>&      syndrome,
                                                    const std::vector<double>& bit_weights) {
    return this->lsd_decode(syndrome, bit_weights, 1, true);
}

std::vector<uint8_t>&
LsdDecoder::lsd_decode(std::vector<uint8_t>&      syndrome,
                       const std::vector<double>& bit_weights,
                       const int                  bits_per_step,
                       const bool                 is_on_the_fly) {
    auto start_time = std::chrono::high_resolution_clock::now();
    this->statistics.clear();
    this->statistics.syndrome = syndrome;

    fill(this->decoding.begin(), this->decoding.end(), 0);

    std::vector<LsdCluster*> clusters;
    std::vector<LsdCluster*> invalid_clusters;
    auto                     global_bit_membership = std::make_shared<std::vector<LsdCluster*>>(
            std::vector<LsdCluster*>(this->pcm.n));
    auto global_check_membership = std::make_shared<std::vector<LsdCluster*>>(
            std::vector<LsdCluster*>(this->pcm.m));
    // timestep to added bits history for stats
    auto* global_timestep_bits_history = new std::unordered_map<int, std::unordered_map<int, std::vector<int>>>{};
    auto  timestep                     = 0;
    for (auto i = 0; i < this->pcm.m; i++) {
        if (syndrome[i] == 1) {
            auto* cl = new LsdCluster(this->pcm, i, global_check_membership, global_bit_membership);
            clusters.push_back(cl);
            invalid_clusters.push_back(cl);
            if (this->do_stats) {
                this->statistics.individual_cluster_stats[cl->cluster_id] = ClusterStatistics();
                cl->global_timestep_bit_history                           = global_timestep_bits_history;
            }
        }
    }

    while (!invalid_clusters.empty()) {
        std::vector<int> new_bits;
        for (auto cl : invalid_clusters) {
            if (cl->active) {
                cl->curr_timestep = timestep; // for stats
                cl->grow_cluster(bit_weights, bits_per_step, is_on_the_fly);
            }
        }
        invalid_clusters.clear();
        for (auto cl : clusters) {
            if (this->do_stats) {
                this->update_growth_stats(cl);
            }
            if (cl->active && !cl->valid) {
                invalid_clusters.push_back(cl);
            } else if (cl->active && cl->valid && this->do_stats) {
                this->statistics.individual_cluster_stats[cl->cluster_id].got_valid_in_timestep = timestep;
            }
            if (do_stats) {
                this->statistics.individual_cluster_stats[cl->cluster_id].active = cl->active;
            }
        }
        std::sort(invalid_clusters.begin(), invalid_clusters.end(),
                  [](const LsdCluster* lhs, const LsdCluster* rhs) {
                      return lhs->bit_nodes.size() < rhs->bit_nodes.size();
                  });
        timestep++; // for stats
    }

    if (lsd_order == 0) {
        this->statistics.lsd_order  = 0;
        this->statistics.lsd_method = osd::OSD_0;
        for (auto cl : clusters) {
            if (do_stats) {
                this->update_final_stats(cl);
            }
            if (cl->active) {
                auto solution                                                      = cl->pluDecomposition.lu_solve(cl->cluster_pcm_syndrome);
                this->statistics.individual_cluster_stats[cl->cluster_id].solution = solution;
                for (auto i = 0; i < solution.size(); i++) {
                    if (solution[i] == 1) {
                        int bit_idx             = cl->cluster_bit_idx_to_pcm_bit_idx[i];
                        this->decoding[bit_idx] = 1;
                    }
                }
            }
        }
    } else {
        this->statistics.lsd_order  = lsd_order;
        this->statistics.lsd_method = this->lsd_method;
        this->apply_lsdw(clusters, lsd_order, bit_weights);
    }
    auto end_time = std::chrono::high_resolution_clock::now();

    if (do_stats) {
        this->statistics.global_timestep_bit_history = *global_timestep_bits_history;
        this->statistics.bit_llrs                    = bit_weights;
    }
    // always take time
    this->statistics.elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(
                                            end_time - start_time)
                                            .count();
    // cleanup
    for (auto cl : clusters) {
        delete cl;
    }
    global_bit_membership->clear();
    global_check_membership->clear();
    delete global_timestep_bits_history;
    return this->decoding;
}

void LsdDecoder::apply_lsdw(const std::vector<LsdCluster*>& clusters,
                            int                             lsd_order,
                            const std::vector<double>& bit_weights, std::size_t timestep) {
    // cluster growth stage
    for (auto cl : clusters) {
        if (cl->active) {
            // first we measure the dimension of each cluster
            auto cluster_dimension = cl->pluDecomposition.not_pivot_cols.size();
            // if the cluster dimension is smaller than the lsd order, we grow the cluster until it reaches
            //  the lsd order. The number of bits added is limited to be at most the lsd order.
            auto cluster_growth_count = 0;
            auto initial_cluster_size = cl->bit_nodes.size();
            while (cluster_dimension < lsd_order &&
                   cluster_growth_count < lsd_order &&
                   cl->bit_nodes.size() < this->pcm.n &&
                   cluster_growth_count <= initial_cluster_size) {
                cl->curr_timestep = timestep; // for stats
                cl->grow_cluster(bit_weights, 1, true);
                cluster_dimension = cl->pluDecomposition.not_pivot_cols.size();
                cluster_growth_count++;
                if (this->do_stats) {
                    this->update_growth_stats(cl);
                }
                timestep++;
            }
        }
    }
    // apply lsd-w to clusters
    for (auto cl : clusters) {
        if (do_stats) {
            this->update_final_stats(cl);
        }
        if (cl->active) {
            cl->sort_non_pivot_cols(bit_weights);
            auto cl_osd_decoder = osd::DenseOsdDecoder(
                    cl->cluster_pcm,
                    cl->pluDecomposition,
                    this->lsd_method,
                    lsd_order,
                    cl->bit_nodes.size(),
                    cl->check_nodes.size(),
                    bit_weights);
            auto res = cl_osd_decoder.osd_decode(cl->cluster_pcm_syndrome);

            for (auto i = 0; i < res.size(); i++) {
                if (res[i] == 1) {
                    int bit_idx             = cl->cluster_bit_idx_to_pcm_bit_idx[i];
                    this->decoding[bit_idx] = 1;
                }
            }
        }
    }
}

} // namespace ldpc::lsd