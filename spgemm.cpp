#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>
#include <mpi.h>
#include <cassert>
#include "functions.h"
#include <cmath>


void spgemm_2d(int m, int p, int n,
               std::vector<std::pair<std::pair<int,int>, int>> &A,
               std::vector<std::pair<std::pair<int,int>, int>> &B,
               std::vector<std::pair<std::pair<int,int>, int>> &C,
               std::function<int(int, int)> plus,
               std::function<int(int, int)> times,
               MPI_Comm row_comm, MPI_Comm col_comm)
{
    int global_rank, global_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    int grid_size = static_cast<int>(std::sqrt(global_size));

    int row_rank = global_rank / grid_size;
    int col_rank = global_rank % grid_size;

    std::map<std::pair<int, int>, int> local_C;
    std::vector<std::pair<std::pair<int,int>, int>> A_i;
    std::vector<std::pair<std::pair<int,int>, int>> B_i;

    for (int k = 0; k < grid_size; ++k)
    {
        if (col_rank == k)
            A_i = A;

        int A_count = (col_rank == k) ? A_i.size() : 0;
        MPI_Bcast(&A_count, 1, MPI_INT, k, row_comm);

        if (col_rank != k)
            A_i.resize(A_count);

        MPI_Bcast(A_i.data(), A_count * sizeof(A[0]), MPI_BYTE, k, row_comm);

        if (row_rank == k)
            B_i = B;

        int B_count = (row_rank == k) ? B_i.size() : 0;
        MPI_Bcast(&B_count, 1, MPI_INT, k, col_comm);

        if (row_rank != k)
            B_i.resize(B_count);

        MPI_Bcast(B_i.data(), B_count * sizeof(B[0]), MPI_BYTE, k, col_comm);

        // Hash B by row for efficient join
        std::unordered_map<int, std::vector<std::pair<int, int>>> B_hash;
        for (const auto &b : B_i)
            B_hash[b.first.first].emplace_back(b.first.second, b.second);

        // Multiply A × B
        for (const auto &a : A_i) {
            int arow = a.first.first;
            int acol = a.first.second;
            int aval = a.second;

            auto it = B_hash.find(acol);
            if (it != B_hash.end()) {
                for (const auto &[bcol, bval] : it->second) {
                    std::pair<int, int> pos(arow, bcol);
                    int prod = times(aval, bval);
                    if (local_C.count(pos))
                        local_C[pos] = plus(local_C[pos], prod);
                    else
                        local_C[pos] = prod;
                }
            }
        }
    }

    // Convert result to COO format
    C.clear();
    for (const auto &entry : local_C)
        C.push_back({entry.first, entry.second});
}
