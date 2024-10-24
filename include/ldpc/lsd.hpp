#pragma once

#include "bp.hpp"
#include "lsd_cluster.hpp"
#include "osd_dense.hpp"

#include <chrono>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ldpc::lsd {

class LsdDecoder {
private:
    ldpc::bp::BpSparse& pcm;
    bool                do_stats;

public:
    std::vector<uint8_t> decoding{};
    Statistics           statistics{};
    int                  bit_count{};
    osd::OsdMethod       lsd_method;
    int                  lsd_order;

    void reset_cluster_stats();

    explicit LsdDecoder(ldpc::bp::BpSparse& parity_check_matrix,
                        osd::OsdMethod      lsdMethod = osd::OsdMethod::COMBINATION_SWEEP,
                        int                 lsd_order = 0) : pcm(parity_check_matrix),
                                             lsd_method(lsdMethod),
                                             lsd_order(lsd_order) {
        this->bit_count = pcm.n;
        this->decoding.resize(this->bit_count);
        this->do_stats = false;
    }

    osd::OsdMethod getLsdMethod() const;

    void setLsdMethod(osd::OsdMethod lsdMethod);

    void set_do_stats(const bool on);

    bool get_do_stats() const;

    void print_cluster_stats();

    void update_growth_stats(const LsdCluster* cl);

    void update_final_stats(const LsdCluster* cl);

    std::vector<uint8_t>& on_the_fly_decode(std::vector<uint8_t>&      syndrome,
                                            const std::vector<double>& bit_weights = NULL_DOUBLE_VECTOR);

    std::vector<uint8_t>&
    lsd_decode(std::vector<uint8_t>&      syndrome,
               const std::vector<double>& bit_weights   = NULL_DOUBLE_VECTOR,
               const int                  bits_per_step = 1,
               const bool                 is_on_the_fly = true);

    void apply_lsdw(const std::vector<LsdCluster*>& clusters,
                    int                             lsd_order,
                    const std::vector<double>& bit_weights, std::size_t timestep = 0);
};

} // namespace ldpc::lsd
