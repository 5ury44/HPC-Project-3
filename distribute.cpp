#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include "functions.h"

void distribute_matrix_2d(int m, int n,
                          std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
                          std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
                          int root, MPI_Comm comm_2d)
{
    int rank, size, coords_rank[2];
    MPI_Comm_rank(comm_2d, &rank);
    MPI_Comm_size(comm_2d, &size);

    int gridsize = std::sqrt(size);

    int block_size_row = m / gridsize;
    int block_size_column = n / gridsize;

    if (rank == root)
    {
        // Step 1: Pre-bucket all entries by destination rank
        std::vector<std::vector<std::pair<std::pair<int, int>, int>>> send_buckets(size);

        for (const auto &elem : full_matrix) {
            int row = elem.first.first;
            int col = elem.first.second;

            int proc_row = std::min(row / block_size_row, gridsize - 1);
            int proc_col = std::min(col / block_size_column, gridsize - 1);

            int dest_coords[2] = {proc_row, proc_col};
            int dest_rank;
            MPI_Cart_rank(comm_2d, dest_coords, &dest_rank);

            send_buckets[dest_rank].push_back(elem);
        }

        // Step 2: Send each processor its block (or assign directly if it's root)
        for (int i = 0; i < size; i++) {
            if (i == root) {
                local_matrix = send_buckets[i];
            } else {
                int buffer_size = send_buckets[i].size();
                MPI_Send(&buffer_size, 1, MPI_INT, i, 0, comm_2d);
                MPI_Send(send_buckets[i].data(), buffer_size * sizeof(send_buckets[i][0]),
                         MPI_BYTE, i, 1, comm_2d);
            }
        }
    }
    else
    {
        // Receive block from root
        int recv_size;
        MPI_Recv(&recv_size, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
        local_matrix.resize(recv_size);
        MPI_Recv(local_matrix.data(), recv_size * sizeof(local_matrix[0]),
                 MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
    }
}
