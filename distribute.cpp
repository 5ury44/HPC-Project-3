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

    int gridDim = static_cast<int>(std::sqrt(szProcesses));
    int blockRows = m / gridDim;
    int blockCols = n / gridDim;

    if (myRank == root) {
        std::vector<std::vector<std::pair<std::pair<int,int>, int>>> procBuckets(szProcesses);

        for (size_t i = 0; i < full_matrix.size(); i++) {
            int curRow = full_matrix[i].first.first;
            int curCol = full_matrix[i].first.second;
        
            int coords[2];
            coords[0] = (curRow < (gridDim - 1) * blockRows) ? curRow / blockRows : gridDim - 1;
            coords[1] = (curCol < (gridDim - 1) * blockCols) ? curCol / blockCols : gridDim - 1;
        
            int targetRank = -1;
            MPI_Cart_rank(comm_2d, coords, &targetRank);
        
            procBuckets[targetRank].push_back(full_matrix[i]);
        }

        for (int p = 0; p < szProcesses; p++) {
            if (p == root) {
                local_matrix = procBuckets[p];
            } else {
                int nItems = procBuckets[p].size();
                MPI_Send(&nItems, 1, MPI_INT, p, 0, comm_2d);
                MPI_Send(procBuckets[p].data(), 
                         nItems * sizeof(procBuckets[p][0]), 
                         MPI_BYTE, p, 1, comm_2d);
            }
        }
    } else {
        int receivedCount = 0;
        MPI_Recv(&receivedCount, 1, MPI_INT, root, 0, comm_2d, MPI_STATUS_IGNORE);
        local_matrix.resize(receivedCount);
        MPI_Recv(local_matrix.data(), 
                 receivedCount * sizeof(local_matrix[0]), 
                 MPI_BYTE, root, 1, comm_2d, MPI_STATUS_IGNORE);
    }
}
