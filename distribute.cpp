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
    int myRank, numProcesses;
    MPI_Comm_rank(comm_2d, &myRank);
    MPI_Comm_size(comm_2d, &numProcesses);

    // Determine grid dimensions given the number of processors (assuming a square grid)
    int gridDim = static_cast<int>(std::sqrt(numProcesses));
    int blockRows = m / gridDim;
    int blockCols = n / gridDim;

    if (myRank == root) {
        // Create a vector of buckets corresponding to each processor in the grid.
        std::vector<std::vector<std::pair<std::pair<int,int>, int>>> procBuckets(numProcesses);

        // Distribute each entry in full_matrix to the appropriate bucket.
        for (size_t idx = 0; idx < full_matrix.size(); idx++) {
            int curRow = full_matrix[idx].first.first;
            int curCol = full_matrix[idx].first.second;
            
            // Determine the processor row and column using integer division and limit to grid bounds.
            int procRow = std::min(curRow / blockRows, gridDim - 1);
            int procCol = std::min(curCol / blockCols, gridDim - 1);

            int coords[2] = {procRow, procCol};
            int targetRank = -1;
            MPI_Cart_rank(comm_2d, coords, &targetRank);
            
            procBuckets[targetRank].push_back(full_matrix[idx]);
        }

        // Send data to each processor.
        for (int proc = 0; proc < numProcesses; proc++) {
            if (proc == root) {
                // Directly copy the bucket for the root processor.
                local_matrix = procBuckets[proc];
            } else {
                // First, send the number of items in this processor's bucket.
                int numItems = procBuckets[proc].size();
                MPI_Send(&numItems, 1, MPI_INT, proc, 0, comm_2d);
                // Then send the actual bucket data.
                MPI_Send(procBuckets[proc].data(), 
                         numItems * sizeof(procBuckets[proc][0]), 
                         MPI_BYTE, proc, 1, comm_2d);
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
