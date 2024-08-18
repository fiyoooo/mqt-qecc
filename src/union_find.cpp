//
// Created by Fiona Froehler on 12.08.24.
//

#include "union_find.hpp"

#include "Code.hpp"
#include "MaxSATDecoder.hpp"

std::vector<int> ldpc::uf::Cluster::maxsat_decode(const std::vector<uint8_t>& syndrome, const std::vector<double>& bit_weights) {
    auto   temp = this->convert_to_matrix(bit_weights);
    gf2Mat cluster_pcm;
    for (int i = 0; i < temp.m; ++i) {
        for (const auto& entry : temp.iterate_row(i)) {
            cluster_pcm[i][entry.col_index] = true;
        }
    }

    gf2Vec cluster_syndrome;
    for (auto check_index : check_nodes) {
        cluster_syndrome.push_back(syndrome[check_index]);
    }

    Code          c(cluster_pcm);
    MaxSATDecoder solver(c);
    solver.preconstructZ3Instance();

    try {
        solver.decode(cluster_syndrome);
        this->valid = true;
    } catch (const std::exception& e) {
        // no solution was found
        this->valid = false;
    }

    // translate maxSAT solution back to original bits
    auto cluster_solution = solver.result.estimBoolVector;
    this->cluster_decoding.clear();
    for (auto i = 0; i < cluster_solution.size(); i++) {
        if (cluster_solution[i] == 1) {
            this->cluster_decoding.push_back(this->matrix_to_cluster_bit_map[i]);
        }
    }
    return this->cluster_decoding;
}

std::vector<uint8_t>& ldpc::uf::UfDecoder::maxsat_decode(const std::vector<uint8_t>& syndrome,
                                                         const std::vector<double>&  bit_weights,
                                                         int                         bits_per_step) {
    fill(this->decoding.begin(), this->decoding.end(), 0);

    std::vector<Cluster*> clusters;
    std::vector<Cluster*> invalid_clusters;
    auto**                global_bit_membership   = new Cluster*[pcm.n]();
    auto**                global_check_membership = new Cluster*[pcm.m]();

    // initialize clusters
    for (auto i = 0; i < this->pcm.m; i++) {
        if (syndrome[i] == 1) {
            auto* cl = new ldpc::uf::Cluster(this->pcm, i, global_check_membership, global_bit_membership);
            clusters.push_back(cl);
            invalid_clusters.push_back(cl);
        }
    }

    // grow clusters and decode with maxsat
    while (!invalid_clusters.empty()) {
        for (auto cl : invalid_clusters) {
            if (cl->active) {
                cl->grow_cluster(bit_weights, bits_per_step);
                auto cluster_decoding = cl->maxsat_decode(syndrome, bit_weights);
            }
        }

        // update invalid clusters
        invalid_clusters.clear();
        for (auto cl : clusters) {
            if (cl->active && !cl->valid) {
                invalid_clusters.push_back(cl);
            }
        }
        std::sort(invalid_clusters.begin(), invalid_clusters.end(), [](const Cluster* lhs, const Cluster* rhs) {
            return lhs->bit_nodes.size() < rhs->bit_nodes.size();
        });
    }

    // apply final decoding
    for (auto cl : clusters) {
        if (cl->active) {
            for (int bit : cl->cluster_decoding) {
                this->decoding[bit] = 1;
            }
        }
        delete cl;
    }
    delete[] global_bit_membership;
    delete[] global_check_membership;
    return this->decoding;
}