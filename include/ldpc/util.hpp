#pragma once

#include <cmath>
#include <iterator>
#include <memory>
#include <vector>

namespace ldpc::util {

std::vector<uint8_t> decimal_to_binary(int decimal_number, int binary_string_length, bool reverse = false);
std::vector<int>     decimal_to_binary_sparse(int decimal_number, int binary_string_length);
std::vector<uint8_t> decimal_to_binary_reverse(int decimal_number, int binary_string_length);
int                  binary_to_decimal(std::vector<uint8_t> binary_number);

} // namespace ldpc::util
