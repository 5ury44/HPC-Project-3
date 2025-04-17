#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include "functions.h"
#include <mpi.h>

void distribute_matrix_2d(int m, int n, 
                          std::vector<std::pair<std::pair<int, int>, int>> &full_matrix,
                          std::vector<std::pair<std::pair<int, int>, int>> &local_matrix,
                          int root, MPI_Comm comm_2d)
{
    int myRank, szProcesses;
    MPI_Comm_rank(comm_2d, &myRank);
    MPI_Comm_size(comm_2d, &szProcesses);

    // Determine grid dimensions given the number of processors (assuming a square grid)
    int gridDim = static_cast<int>(std::sqrt(szProcesses));
    int blockRows = m / gridDim;
    int blockCols = n / gridDim;

    if (myRank == root) {
        // Create a vector of buckets corresponding to each processor in the grid.
        std::vector<std::vector<std::pair<std::pair<int,int>, int>>> procBuckets(szProcesses);

        // Distribute each entry in full_matrix to the appropriate bucket.
        for (size_t i = 0; i < full_matrix.size(); i++) {
            int curRow = full_matrix[i].first.first;
            int curCol = full_matrix[i].first.second;
        
            // Safer handling of edge coordinates
            int coords[2];
            coords[0] = (curRow < (gridDim - 1) * blockRows) ? curRow / blockRows : gridDim - 1;
            coords[1] = (curCol < (gridDim - 1) * blockCols) ? curCol / blockCols : gridDim - 1;
        
            int targetRank = -1;
            MPI_Cart_rank(comm_2d, coords, &targetRank);
        
            procBuckets[targetRank].push_back(full_matrix[i]);
        }

        // Send data to each processor.
        for (int p = 0; p < szProcesses; p++) {
            if (p == root) {
                // Directly copy the bucket for the root processor.
                local_matrix = procBuckets[p];
            } else {
                // First, send the number of items in this processor's bucket.
                int nItems = procBuckets[p].size();
                MPI_Send(&nItems, 1, MPI_INT, p, 0, comm_2d);
                // Then send the actual bucket data.
                MPI_Send(procBuckets[p].data(), 
                         nItems * sizeof(procBuckets[p][0]), 
                         MPI_BYTE, p, 1, comm_2d);
            }
        }
    } else {
        // Non-root processors receive the block size first and then the bucket.
        int receivedCount = 0;
        MPI_Recv(&receivedCount, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
        local_matrix.resize(receivedCount);
        MPI_Recv(local_matrix.data(), 
                 receivedCount * sizeof(local_matrix[0]), 
                 MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
    }
}
