#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include <functional>

using MatrixEntry = std::pair<std::pair<int, int>, int>;
using SparseMatrix = std::vector<MatrixEntry>;

// Helper to broadcast a sparse matrix
void broadcast_matrix(SparseMatrix &matrix, int root, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int count = matrix.size();
    MPI_Bcast(&count, 1, MPI_INT, root, comm);

    if (rank != root) {
        matrix.resize(count);
    }
    MPI_Bcast(matrix.data(), count * sizeof(MatrixEntry), MPI_BYTE, root, comm);
}

// Hash B matrix rows by key to accelerate matching
std::unordered_map<int, std::vector<std::pair<int, int>>> hash_matrix_by_row(const SparseMatrix &B) {
    std::unordered_map<int, std::vector<std::pair<int, int>>> hash;
    for (const auto &entry : B) {
        int row = entry.first.first;
        int col = entry.first.second;
        int val = entry.second;
        hash[row].emplace_back(col, val);
    }
    return hash;
}

// Perform local multiplication between A and hashed B
void local_spgemm(const SparseMatrix &A,
                  const std::unordered_map<int, std::vector<std::pair<int, int>>> &B_hash,
                  std::map<std::pair<int, int>, int> &local_C,
                  std::function<int(int, int)> plus,
                  std::function<int(int, int)> times) {
    for (const auto &a : A) {
        int arow = a.first.first;
        int acol = a.first.second;
        int aval = a.second;

        auto it = B_hash.find(acol);
        if (it != B_hash.end()) {
            for (const auto &[bcol, bval] : it->second) {
                std::pair<int, int> pos(arow, bcol);
                int prod = times(aval, bval);
                if (local_C.find(pos) == local_C.end()) {
                    local_C[pos] = prod;
                } else {
                    local_C[pos] = plus(local_C[pos], prod);
                }
            }
        }
    }
}

// Main 2D SpGEMM driver
void spgemm_2d(int m, int p, int n,
               SparseMatrix &A,
               SparseMatrix &B,
               SparseMatrix &C,
               std::function<int(int, int)> plus,
               std::function<int(int, int)> times,
               MPI_Comm row_comm, MPI_Comm col_comm)
{
    int row_rank, row_size, col_rank, col_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);

    std::map<std::pair<int, int>, int> local_C;

    for (int k = 0; k < row_size; ++k) {
        SparseMatrix A_i, B_i;

        if (row_rank == k) A_i = A;
        broadcast_matrix(A_i, k, row_comm);

        if (col_rank == k) B_i = B;
        broadcast_matrix(B_i, k, col_comm);

        auto B_hash = hash_matrix_by_row(B_i);
        local_spgemm(A_i, B_hash, local_C, plus, times);
    }

    C.clear();
    for (const auto &entry : local_C) {
        C.emplace_back(entry.first, entry.second);
    }
}
