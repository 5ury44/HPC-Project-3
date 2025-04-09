#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include "functions.h"

#include <vector>
#include <unordered_map>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include "functions.h"

void distribute_matrix_2d(int m, int n,
                          std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
                          std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
                          int root, MPI_Comm comm_2d) {
    int rank, size;
    MPI_Comm_rank(comm_2d, &rank);
    MPI_Comm_size(comm_2d, &size);

    int coords[2];
    int dims[2], periods[2];
    MPI_Cart_get(comm_2d, 2, dims, periods, coords);

    int dim = static_cast<int>(std::sqrt(size));
    assert(dim * dim == size);  // Ensure square processor grid

    // Calculate block dimensions for load balancing (handle uneven sizes)
    auto get_block_bounds = [](int global_size, int dim, int coord) {
        int base = global_size / dim;
        int rem = global_size % dim;
        int start = coord * base + std::min(coord, rem);
        int end = start + base + (coord < rem ? 1 : 0);
        return std::make_pair(start, end);
    };

    int val, destRank;
    std::vector<int> sendCounts(size, 0);
    std::vector<int> displs(size, 0);
    std::vector<int> sendBuff;

    if (rank == root) {
        std::vector<std::vector<int>> packed(size); // [rank] -> flat int buffer

        for (const auto &entry : full_matrix) {
            int i = entry.first.first;
            int j = entry.first.second;
            int val = entry.second;

            // Determine destination rank by block coords
            int proc_row = std::min(i * dim / m, dim - 1);
            int proc_col = std::min(j * dim / n, dim - 1);

            int destCoords[2] = {proc_row, proc_col};
            MPI_Cart_rank(comm_2d, destCoords, &destRank);

            packed[destRank].push_back(i);
            packed[destRank].push_back(j);
            packed[destRank].push_back(val);
        }

        // Build flat send buffer and metadata
        for (int r = 0; r < size; ++r) {
            sendCounts[r] = packed[r].size();
            displs[r] = (r == 0 ? 0 : displs[r - 1] + sendCounts[r - 1]);
            sendBuff.insert(sendBuff.end(), packed[r].begin(), packed[r].end());
        }
    }

    int recvCount;
    MPI_Scatter(sendCounts.data(), 1, MPI_INT, &recvCount, 1, MPI_INT, root, comm_2d);

    std::vector<int> recvBuf(recvCount);
    MPI_Scatterv(sendBuff.data(), sendCounts.data(), displs.data(), MPI_INT,
                 recvBuf.data(), recvCount, MPI_INT, root, comm_2d);

    // Rebuild COO triplets from flat buffer
    for (int i = 0; i < recvBuf.size(); i += 3) {
        int row = recvBuf[i];
        int col = recvBuf[i + 1];
        int val = recvBuf[i + 2];
        local_matrix.emplace_back(std::make_pair(std::make_pair(row, col), val));
    }


}
