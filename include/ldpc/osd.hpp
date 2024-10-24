#pragma once

#include "bp.hpp"
#include "gf2sparse_linalg.hpp"
#include "sort.hpp"
#include "util.hpp"

#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

namespace ldpc::osd {

enum OsdMethod {
    OSD_OFF,
    OSD_0,
    EXHAUSTIVE,
    COMBINATION_SWEEP
};

class OsdDecoder {
public:
    OsdMethod                                             osd_method;
    int                                                   osd_order;
    int                                                   k, bit_count, check_count;
    ldpc::bp::BpSparse&                                   pcm;
    std::vector<double>&                                  channel_probabilities;
    std::vector<uint8_t>                                  osd0_decoding;
    std::vector<uint8_t>                                  osdw_decoding;
    std::vector<std::vector<uint8_t>>                     osd_candidate_strings;
    std::vector<int>                                      column_ordering;
    ldpc::gf2sparse_linalg::RowReduce<ldpc::bp::BpEntry>* LuDecomposition;

    OsdDecoder(
            ldpc::bp::BpSparse&  parity_check_matrix,
            OsdMethod            osd_method,
            int                  osd_order,
            std::vector<double>& channel_probs);

    int osd_setup();

    ~OsdDecoder();

    std::vector<uint8_t>& decode(std::vector<uint8_t>& syndrome, std::vector<double>& log_prob_ratios);
};

} // namespace ldpc::osd
