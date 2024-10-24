#pragma once

#include <vector>

namespace ldpc::sort {

struct str {
    double value;
    int    index;
};

unsigned long long ncr(int n, int k);
int                cmp(const void* a, const void* b);
void               soft_decision_col_sort(std::vector<double>& soft_decisions, std::vector<int>& cols, int N);

} // end namespace ldpc::sort
