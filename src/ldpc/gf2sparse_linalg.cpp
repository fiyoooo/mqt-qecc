#include "ldpc/gf2sparse_linalg.hpp"

std::vector<int> ldpc::gf2sparse_linalg::NULL_INT_VECTOR = {};

// cython helper
ldpc::sparse_matrix_base::CsrMatrix
ldpc::gf2sparse_linalg::cy_row_complement_basis(ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry>* mat) {
    return row_complement_basis(*mat).to_csr_matrix();
}

// cython helper
ldpc::sparse_matrix_base::CsrMatrix
ldpc::gf2sparse_linalg::cy_kernel(ldpc::gf2sparse::GF2Sparse<ldpc::gf2sparse::GF2Entry>* mat) {
    return kernel(*mat).to_csr_matrix();
}